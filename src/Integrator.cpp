
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <utility>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Sequence.h>

#include "ModelFactory.h"

#ifdef DEBUG_BAM_OUTPUT
#include <cram/md5.h>
#include <pbbam/BamWriter.h>
#include <pbbam/BamFile.h>
#include <pbbam/MD5.h>
#include <pacbio/consensus/align/LinearAlignment.h>
#endif

namespace PacBio {
namespace Consensus {

constexpr double NaN = -std::numeric_limits<double>::quiet_NaN();

std::set<std::string> SupportedChemistries() { return ModelFactory::SupportedChemistries(); }
IntegratorConfig::IntegratorConfig(const double minZScore, const double scoreDiff)
    : MinZScore{minZScore}, ScoreDiff{scoreDiff}
{
    if (ScoreDiff < 0) {
        throw std::runtime_error("Score diff must be > 0");
    }
}

AbstractIntegrator::AbstractIntegrator(const IntegratorConfig& cfg) : cfg_{cfg} {}
AbstractIntegrator::AbstractIntegrator(AbstractIntegrator&& ai)
    : cfg_{ai.cfg_}, evals_{std::move(ai.evals_)}
{
}

AbstractIntegrator::~AbstractIntegrator() {}
AddReadResult AbstractIntegrator::AddRead(std::unique_ptr<AbstractTemplate>&& tpl,
                                          const MappedRead& read)
{
    if (read.TemplateEnd <= read.TemplateStart || read.TemplateEnd - read.TemplateStart < 2)
        throw std::invalid_argument("template span < 2!");

    if (read.Length() < 2) throw std::invalid_argument("read span < 2!");

    evals_.emplace_back(Evaluator(std::move(tpl), read, cfg_.MinZScore, cfg_.ScoreDiff));

    const auto status = evals_.back().Status();

    if (status == EvaluatorState::ALPHA_BETA_MISMATCH)
        return AddReadResult::ALPHA_BETA_MISMATCH;
    else if (status == EvaluatorState::POOR_ZSCORE)
        return AddReadResult::POOR_ZSCORE;

    assert(status == EvaluatorState::VALID);

    return AddReadResult::SUCCESS;
}

double AbstractIntegrator::LL(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));
    double ll = 0.0;
    for (auto& eval : evals_) {
        if (eval.Strand() == StrandEnum::FORWARD)
            ll += eval.LL(fwdMut);
        else if (eval.Strand() == StrandEnum::REVERSE)
            ll += eval.LL(revMut);
    }
    return ll;
}

double AbstractIntegrator::LL() const
{
    double ll = 0.0;
    for (const auto& eval : evals_) {
        if (eval) ll += eval.LL();
    }
    return ll;
}

std::vector<double> AbstractIntegrator::LLs(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));
    std::vector<double> lls;
    for (auto& eval : evals_) {
        if (eval.Strand() == StrandEnum::FORWARD)
            lls.push_back(eval.LL(fwdMut));
        else if (eval.Strand() == StrandEnum::REVERSE)
            lls.push_back(eval.LL(revMut));
        else {
            // inactive reads get a strand of "unmapped"
            lls.push_back(NaN);
        }
    }
    return lls;
}

std::vector<double> AbstractIntegrator::LLs() const
{
    std::vector<double> lls;
    for (auto& eval : evals_) {
        if (eval)
            lls.push_back(eval.LL());
        else
            lls.push_back(NaN);
    }
    return lls;
}

std::vector<std::string> AbstractIntegrator::ReadNames() const
{
    std::vector<std::string> readNames;
    for (const Evaluator& eval : evals_) {
        readNames.push_back(eval.ReadName());
    }
    return readNames;
}

double AbstractIntegrator::AvgZScore() const
{
    double mean = 0.0, var = 0.0;
    size_t n = 0;
    for (const auto& eval : evals_) {
        if (eval) {
            double m, v;
            std::tie(m, v) = eval.NormalParameters();
            mean += m;
            var += v;
            ++n;
        }
    }
    return (LL() / n - mean / n) / std::sqrt(var / n);
}

std::vector<double> AbstractIntegrator::ZScores() const
{
    std::vector<double> results;
    results.reserve(evals_.size());
    for (const auto& eval : evals_) {
        if (eval) {
            double mean, var;
            std::tie(mean, var) = eval.NormalParameters();
            results.emplace_back((eval.LL() - mean) / std::sqrt(var));
        } else
            results.emplace_back(std::numeric_limits<double>::quiet_NaN());
    }
    return results;
}
std::vector<std::pair<double, double>> AbstractIntegrator::NormalParameters() const
{
    std::vector<std::pair<double, double>> results;
    results.reserve(evals_.size());
    for (const auto& eval : evals_) {
        if (eval) {
            results.emplace_back(eval.NormalParameters());
        } else {
            auto nan =std::numeric_limits<double>::quiet_NaN();
            results.emplace_back(std::pair<double, double>(nan, nan));
        }
    }
    return results;
}
    
    
Mutation AbstractIntegrator::ReverseComplement(const Mutation& mut) const
{
    return Mutation(mut.Type, TemplateLength() - mut.End(), Complement(mut.Base));
}

MonoMolecularIntegrator::MonoMolecularIntegrator(const std::string& tpl,
                                                 const IntegratorConfig& cfg, const SNR& snr,
                                                 const std::string& model)
    : AbstractIntegrator(cfg)
    , mdl_{model}
    , snr_{snr}
    , fwdTpl_(tpl, ModelFactory::Create(mdl_, snr_))
    , revTpl_(::PacBio::Consensus::ReverseComplement(tpl), ModelFactory::Create(mdl_, snr_))
{
}

MonoMolecularIntegrator::MonoMolecularIntegrator(MonoMolecularIntegrator&& mmi)
    : AbstractIntegrator(std::move(mmi))
    , mdl_{mmi.mdl_}
    , snr_{mmi.snr_}
    , fwdTpl_{std::move(mmi.fwdTpl_)}
    , revTpl_{std::move(mmi.revTpl_)}
{
}

AddReadResult MonoMolecularIntegrator::AddRead(const MappedRead& read)
{
    if (read.Model != mdl_) throw std::invalid_argument("invalid model for integrator!");
    if (read.SignalToNoise != snr_) throw std::invalid_argument("invalid SNR for integrator!");

    if (read.Strand == StrandEnum::FORWARD)
        return AbstractIntegrator::AddRead(
            std::unique_ptr<AbstractTemplate>(new VirtualTemplate(
                fwdTpl_, read.TemplateStart, read.TemplateEnd, read.PinStart, read.PinEnd)),
            read);

    else if (read.Strand == StrandEnum::REVERSE)
        return AbstractIntegrator::AddRead(
            std::unique_ptr<AbstractTemplate>(
                new VirtualTemplate(revTpl_, TemplateLength() - read.TemplateEnd,
                                    TemplateLength() - read.TemplateStart, read.PinEnd, read.PinStart)),
            read);

    throw std::invalid_argument("read is unmapped!");
}

size_t MonoMolecularIntegrator::TemplateLength() const { return fwdTpl_.TrueLength(); }
char MonoMolecularIntegrator::operator[](const size_t i) const { return fwdTpl_[i].Base; }
MonoMolecularIntegrator::operator std::string() const
{
    std::string result;

    result.resize(fwdTpl_.Length());

    for (size_t i = 0; i < fwdTpl_.Length(); ++i)
        result[i] = fwdTpl_[i].Base;

    return result;
}

double MonoMolecularIntegrator::LL(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));
    fwdTpl_.Mutate(fwdMut);
    revTpl_.Mutate(revMut);
    const double ll = AbstractIntegrator::LL(fwdMut);
    fwdTpl_.Reset();
    revTpl_.Reset();
    return ll;
}

void MonoMolecularIntegrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    fwdTpl_.ApplyMutation(fwdMut);
    revTpl_.ApplyMutation(revMut);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandEnum::FORWARD)
            eval.ApplyMutation(fwdMut);
        else if (eval.Strand() == StrandEnum::REVERSE)
            eval.ApplyMutation(revMut);
    }

    assert(fwdTpl_.Length() == revTpl_.Length());

#ifndef NDEBUG
    std::string fwd;
    std::string rev;

    for (size_t i = 0; i < TemplateLength(); ++i) {
        fwd.push_back(fwdTpl_[i].Base);
        rev.push_back(revTpl_[i].Base);
    }

#endif

    assert(fwd == ::PacBio::Consensus::ReverseComplement(rev));
}

void MonoMolecularIntegrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_.ApplyMutations(fwdMuts);
    revTpl_.ApplyMutations(&revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandEnum::FORWARD)
            eval.ApplyMutations(fwdMuts);
        else if (eval.Strand() == StrandEnum::REVERSE)
            eval.ApplyMutations(&revMuts);
    }

    assert(fwdTpl_.Length() == revTpl_.Length());

#ifndef NDEBUG
    std::string fwd;
    std::string rev;

    for (size_t i = 0; i < TemplateLength(); ++i) {
        fwd.push_back(fwdTpl_[i].Base);
        rev.push_back(revTpl_[i].Base);
    }
#endif

    assert(fwd == ::PacBio::Consensus::ReverseComplement(rev));
}

MultiMolecularIntegrator::MultiMolecularIntegrator(const std::string& tpl,
                                                   const IntegratorConfig& cfg)
    : AbstractIntegrator(cfg), fwdTpl_{tpl}, revTpl_{::PacBio::Consensus::ReverseComplement(tpl)}
{
}

AddReadResult MultiMolecularIntegrator::AddRead(const MappedRead& read)
{
    return AbstractIntegrator::AddRead(GetTemplate(read, read.SignalToNoise), read);
}

size_t MultiMolecularIntegrator::TemplateLength() const { return fwdTpl_.length(); }
char MultiMolecularIntegrator::operator[](const size_t i) const { return fwdTpl_[i]; }
MultiMolecularIntegrator::operator std::string() const { return fwdTpl_; }
void MultiMolecularIntegrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    std::vector<Mutation> fwdMuts = {fwdMut};
    std::vector<Mutation> revMuts = {revMut};

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, &fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandEnum::FORWARD)
            eval.ApplyMutation(fwdMut);
        else if (eval.Strand() == StrandEnum::REVERSE)
            eval.ApplyMutation(revMut);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Consensus::ReverseComplement(revTpl_));
}

void MultiMolecularIntegrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandEnum::FORWARD)
            eval.ApplyMutations(fwdMuts);
        else if (eval.Strand() == StrandEnum::REVERSE)
            eval.ApplyMutations(&revMuts);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Consensus::ReverseComplement(revTpl_));
}

std::unique_ptr<AbstractTemplate> MultiMolecularIntegrator::GetTemplate(const MappedRead& read,
                                                                        const SNR& snr)
{
    const size_t len = read.TemplateEnd - read.TemplateStart;

    if (read.Strand == StrandEnum::FORWARD) {
        const size_t start = read.TemplateStart;
        const size_t end = read.TemplateEnd;

        return std::unique_ptr<AbstractTemplate>(
            new Template(fwdTpl_.substr(start, len), ModelFactory::Create(read.Model, snr), start,
                         end, read.PinStart, read.PinEnd));
    } else if (read.Strand == StrandEnum::REVERSE) {
        const size_t start = revTpl_.size() - read.TemplateEnd;
        const size_t end = revTpl_.size() - read.TemplateStart;

        return std::unique_ptr<AbstractTemplate>(
            new Template(revTpl_.substr(start, len), ModelFactory::Create(read.Model, snr), start,
                         end, read.PinEnd, read.PinStart));
    }

    throw std::invalid_argument("read is unmapped!");
}

#ifdef DEBUG_BAM_OUTPUT
std::string AlignmentCigar(const std::string& tpl, const std::string& query)
{
    PairwiseAlignment* p = AlignLinear(tpl, query);
    return p->Cigar();
}

size_t TemplateSpan(const MappedRead& read)
{
    return read.TemplateEnd - read.TemplateStart;
}

int32_t ReadToFlag(const MappedRead& read)
{
    uint32_t flag = 0;
    if (read.Strand == StrandEnum::REVERSE)
        flag |= static_cast<uint32_t>(16);
    return flag;
}

std::string ReadToSequence(const MappedRead& read)
{
    if (read.Strand == StrandEnum::REVERSE) {
        return ::PacBio::Consensus::ReverseComplement(read.Seq);
    } else {
        return read.Seq;
    }
}

std::string MultiMolecularIntegrator::ReadToCigar(const MappedRead& read) const
{
    const size_t tlen = read.TemplateEnd - read.TemplateStart;
    const std::string tpl = std::string(fwdTpl_).substr(read.TemplateStart, tlen);
    if (read.Strand == StrandEnum::REVERSE) {
        const std::string rcSeq = ::PacBio::Consensus::ReverseComplement(read.Seq);
        return AlignmentCigar(tpl, rcSeq);
    } else {
        return AlignmentCigar(tpl, read.Seq);
    }
}

std::vector<PacBio::BAM::BamRecordImpl> MultiMolecularIntegrator::ToBamRecords(const int32_t rid) const
{
    using namespace PacBio::BAM;
    std::vector<BamRecordImpl> records;
    for (auto& eval : evals_) {
        // Create a BamRecord for every (A) valid eval with (B) a valid mapped read
        if (eval) {
            if (const auto& read = eval.Read()) {
                records.emplace_back();
                records.back().Name(read->Name);                                   // QNAME
                records.back().Flag(ReadToFlag(*read));                            // FLAG
                records.back().ReferenceId(rid);                                   // RNAME
                records.back().Position(read->TemplateStart);                      // POS
                records.back().CigarData(ReadToCigar(*read));                      // CIGAR
                records.back().InsertSize(TemplateSpan(*read));                    // TLEN
                records.back().SetSequenceAndQualities(ReadToSequence(*read));     // SEQ
                records.back().AddTag("ip", Tag{read->IPD});                       // TAG "ip" for IPD
                records.back().AddTag("mo", Tag{read->Model});                     // TAG "mo" for sequencing model
                records.back().AddTag("pw", Tag{read->PulseWidth});                // TAG "pw" for PulseWidth
                records.back().AddTag("sn", Tag{read->SignalToNoise.ToVector()});  // TAG "sn" for SNR
            }
        }
    }

    return records;
}

void MultiMolecularIntegrator::WriteBamFile(const std::string& filepath,
                                            const std::string& name,
                                            const bool sorted,
                                            const bool indexed) const
{
    using namespace PacBio::BAM;

    BamHeader header;
    // Add SQ entry with the template information
    SequenceInfo si(name, std::to_string(fwdTpl_.size()));
    si.Checksum(MD5Hash(fwdTpl_));
    header.AddSequence(si);
    // Add CO entry with the raw sequence
    header.AddComment(name + "\t" + std::string(fwdTpl_));

    {   // Blocked to ensure flushing and completion of file I/O
        BamWriter writer(filepath, header);
        const int32_t sid = header.SequenceId(name);
        std::vector<BamRecordImpl> records = ToBamRecords(sid);
    
        if (sorted)
            std::stable_sort(records.begin(), records.end(),
                [](const BamRecordImpl& lhs, const BamRecordImpl& rhs) { 
                    return lhs.Position() < rhs.Position(); 
                });

        for (const auto& record : records)
            writer.Write(record);
    }

    if (indexed) {
        BamFile f(filepath);
        f.CreateStandardIndex();
    }
}

void MultiMolecularIntegrator::WriteReferenceFasta(const std::string& filepath, 
                                                   const std::string& name) const
{
    const size_t COLUMN_SIZE = 60;
    std::string tpl(fwdTpl_);
    std::ofstream fastaFile;
    fastaFile.open(filepath);
    fastaFile << ">" << name << std::endl;
    for (size_t i = 0; i < tpl.size();) {
        fastaFile << tpl.substr(i, COLUMN_SIZE) << std::endl;
        i += COLUMN_SIZE;
    }
    fastaFile.close();
}

void MultiMolecularIntegrator::WriteAlignments(const std::string& filePrefix) const
{
    WriteReferenceFasta(filePrefix + ".fasta");
    WriteBamFile(filePrefix + ".bam");
}

#endif

}  // namespace Consensus
}  // namespace PacBio
