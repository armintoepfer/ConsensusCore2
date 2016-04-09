
#include <cassert>
#include <cmath>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>

#include "../ModelFactory.h"

namespace PacBio {
namespace Consensus {
namespace {

class SP1C1BetaNoCovModel : public ModelConfig
{
    REGISTER_MODEL(SP1C1BetaNoCovModel);

public:
    SP1C1BetaNoCovModel(const SNR& snr);
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double BaseEmissionPr(MoveType move, char from, char to) const;
    double CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const;
    double UndoCounterWeights(size_t nEmissions) const;

    static std::string Name() { return "S/P1-C1/beta"; }
private:
    SNR snr_;
};
    REGISTER_MODEL_IMPL(SP1C1BetaNoCovModel);
    
namespace {

    double matchPmf[16][4]  = {
        
        {    0.987467583,   0.00815966815,   0.00384390417,  0.000367362587},
        {   0.0489892305,     0.950237525,  0.000194606805,  0.000403018444},
        {  0.00113256009,  0.000735859361,      0.98133815,    0.0165934326},
        { 0.000371309534,  0.000329000326,     0.137274406,     0.861852002},
        {    0.992075096,   0.00424731122,   0.00303799373,  0.000461811116},
        {   0.0806231868,     0.919001207,  3.49569275e-05,  0.000178293692},
        { 0.000698558094,  0.000220067471,     0.988339038,    0.0106008344},
        { 0.000442480786,  0.000328643365,     0.138124214,     0.860896452},
        {    0.993153719,    0.0045494866,   0.00201333115,  0.000101441609},
        {   0.0302546227,      0.96913138,    0.0002153392,  0.000257838502},
        {   0.0363827704,  0.000212301565,     0.955128655,   0.00811056487},
        { 0.000774154097,  0.000597866658,    0.0862139813,     0.912223159},
        {    0.989464097,    0.0070582318,    0.0031269763,  0.000173000262},
        {   0.0288085182,     0.970578603,  0.000209353191,    0.0002213711},
        { 0.000921896016,  0.000350856071,     0.988495376,     0.010056583},
        {   0.0660494702,    0.0313132609,     0.116608078,     0.785883684}
    };
    
    double stickPmf[16][4]  = {
        
        { 0.000686185232,     0.145034077,     0.737394958,     0.107964371},
        {    0.323009777,  0.000278091164,     0.569278195,     0.103818752},
        {    0.357565343,     0.193145965,    0.0020635152,      0.42039948},
        {    0.436992797,    0.0750760383,     0.477715765,  0.000729671403},
        { 0.000299731481,     0.316182166,     0.636891797,    0.0427297961},
        {    0.235488828,  0.000169974812,     0.696616791,    0.0656809424},
        {    0.502200864,     0.409471485,  0.000321781565,    0.0838227086},
        {    0.309125231,     0.396910467,      0.29009318,  0.000276508732},
        {  0.00144836803,     0.220893869,     0.472958325,     0.285870653},
        {    0.163607721,  0.000155629485,     0.751968553,    0.0822449133},
        {    0.636274951,    0.0934662256,   0.00134266579,     0.251967978},
        {    0.575088882,    0.0951877758,     0.318483776,  0.000802826171},
        { 0.000353978738,    0.0917892761,     0.529431691,      0.37382333},
        {    0.175702157,  0.000205233362,     0.437185531,     0.384239045},
        {    0.361426238,    0.0236117117,  0.000441153653,     0.608785899},
        {    0.203671379,    0.0207588027,     0.770504977,  0.000449246693}
    };
    
    
    double transProbs[16][4] = {
        
        {    0.903872612,   0.00849349392,    0.0161886892,    0.0714452047}, // CTX - 0
        {    0.935892343,   0.00223659735,    0.0453169902,    0.0165540699}, // CTX - 1
        {    0.904872894,     0.031439782,   0.00653179264,    0.0571555318}, // CTX - 2
        {     0.89524949,   0.00105350139,    0.0162071563,    0.0874898521}, // CTX - 3
        {    0.849642076,    0.0656928784,    0.0383541161,    0.0463109295}, // CTX - 4
        {    0.852505645,     0.039171176,     0.067580177,    0.0407430022}, // CTX - 5
        {     0.90198018,    0.0398023577,    0.0300558285,    0.0281616343}, // CTX - 6
        {    0.856917957,   0.00490832296,    0.0492250567,    0.0889486629}, // CTX - 7
        {    0.901407759,       0.0655271,   0.00833807707,    0.0247270637}, // CTX - 8
        {    0.910359193,    0.0174218349,    0.0630239024,   0.00919506964}, // CTX - 9
        {    0.910379911,    0.0171579767,   0.00872024341,     0.063741869}, // CTX - 10
        {    0.887196572,    0.0146617485,    0.0159119001,    0.0822297793}, // CTX - 11
        {    0.850105939,    0.0745664961,    0.0322468406,    0.0430807245}, // CTX - 12
        {    0.897221167,    0.0330637119,    0.0606995278,   0.00901559332}, // CTX - 13
        {    0.889079776,    0.0520036086,    0.0265717884,    0.0323448268}, // CTX - 14
        {    0.825617274,    0.0688091873,    0.0290846494,    0.0764888889} // CTX - 15
    };

    

    
}

SP1C1BetaNoCovModel::SP1C1BetaNoCovModel(const SNR& snr) : snr_(snr) {}
std::vector<TemplatePosition> SP1C1BetaNoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;
    
    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    
    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t cur = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (cur > 3) throw std::invalid_argument("invalid character in sequence!");
        const auto row = (prev << 2) | cur;
        const auto params = transProbs[row];

        result.emplace_back(TemplatePosition{
            tpl[i - 1],
            params[0],  // match
            params[1],  // branch
            params[2],  // stick
            params[3]   // deletion
        });
        prev = cur;
    }

    result.emplace_back(TemplatePosition{tpl.back(), 1.0, 0.0, 0.0, 0.0});

    return result;
}

double SP1C1BetaNoCovModel::BaseEmissionPr(MoveType move, const char from, const char to) const
{
    // All emissions are now "Covariate Emissions"
    assert(move != MoveType::DELETION);
    return 1.0;
}

double SP1C1BetaNoCovModel::CovEmissionPr(MoveType move, uint8_t nuc, char from, char to) const
{
    auto ito = detail::TranslationTable[static_cast<uint8_t>(to)];
    auto ifrom = detail::TranslationTable[static_cast<uint8_t>(from)];
    ifrom = ifrom > 3 ? 0 : ifrom;
    ito = ito > 3 ? 0 : ito;
    const auto row = (ifrom << 2) | ito;
    
    if (ito > 3 || nuc > 3 || ifrom > 3) throw std::invalid_argument("invalid character in sequence");

    if (move == MoveType::BRANCH)
        return 1.0;
    else if (move == MoveType::STICK)
        return stickPmf[row][nuc];
    else if (move == MoveType::MATCH)
        return matchPmf[row][nuc];

    throw std::invalid_argument("unknown move type!");
}

double SP1C1BetaNoCovModel::UndoCounterWeights(const size_t nEmissions) const { return 0; }
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
