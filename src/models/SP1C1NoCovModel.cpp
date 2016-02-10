
#include <cassert>
#include <cmath>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>

#include "../ModelFactory.h"

namespace PacBio {
namespace Consensus {
namespace {

class SP1C1NoCovModel : public ModelConfig
{
    REGISTER_MODEL(SP1C1NoCovModel);

public:
    SP1C1NoCovModel(const SNR& snr);
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double BaseEmissionPr(MoveType move, char from, char to) const;
    double CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const;
    double UndoCounterWeights(size_t nEmissions) const;

    static std::string Name() { return "S/P1-C1"; }
private:
    SNR snr_;
};

REGISTER_MODEL_IMPL(SP1C1NoCovModel);

    double matchPmf[8][4] = {
        {0.97516262936181, 0.0226803868009865, 0.000904104876493169, 0.00125287896071021},
        {0.00777521848583749, 0.992069805261659, 2.55556083348524e-05, 0.00012942064416881},
        {0.000729957177852375, 0.000960606317664673, 0.955404682823147, 0.0429047536813359},
        {7.40563382904245e-05, 0.000439417167829446, 0.0132299436802719, 0.986256582813608},
        {0.982115284751547, 0.0162249242121629, 0.000406094224056689, 0.00125369681223344},
        {0.00490144805537043, 0.995098025213339, 3.70909634149552e-08, 4.89640326665825e-07},
        {0.00078566308261712, 0.000149098078827161, 0.929347787240117, 0.0697174515984391},
        {5.90604186208778e-05, 0.000590402863495049, 0.0139392468542019, 0.985411289863682},
    };


    double stickPmf[8][4] = {
        {0, 0.187144120055573, 0.146125761154764, 0.666730118789663},
        {0.152613825438472, 0, 0.19898091881288, 0.648405255748648},
        {0.263478215073993, 0.158705661412119, 0, 0.577816123513889},
        {0.339652966590421, 0.245694728009963, 0.414652305399616, 0},
        {0, 0.194445956334637, 0.172385881192111, 0.633168162473252},
        {0.140884762626049, 0, 0.1783771029363, 0.680738134437651},
        {0.173696284650599, 0.0791663237389683, 0, 0.747137391610433},
        {0.629630101703802, 0.251134780702687, 0.119235117593511, 0},
    };



double branchPmf[8][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
                          {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

double SP1C1NoCovParams[8][4] = {
    // Match, Branch, Stick, Delete
    {0.912459935537425,0.0082391897431644,0.0210795390699932,0.05822133564941},
    {0.858034245282967,0.0294528340974306,0.0759756484128592,0.03653727220674},
    {0.913573628851584,0.0120377377203805,0.0199863697578698,0.05440226367016},
    {0.842605304895112,0.100035283428237,0.023547739920972,0.0338116717556791},
    {0.879191827886969,0.0646288774918684,0.0262320142719883,0.02994728034917},
    {0.862217422406686,0.0448665502306592,0.0561215546663366,0.03679447269631},
    {0.911565962931167,0.0291072970021701,0.0235786421895474,0.03574809787711},
    {0.922332894278559,0.0220364924120602,0.0158119417104142,0.03981867159896}


};

SP1C1NoCovModel::SP1C1NoCovModel(const SNR& snr) : snr_(snr) {}
std::vector<TemplatePosition> SP1C1NoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t b = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];

        if (b > 3) throw std::invalid_argument("invalid character in sequence!");

        const bool hpAdd = tpl[i - 1] == tpl[i] ? 0 : 4;
        const auto params = SP1C1NoCovParams[b + hpAdd];

        result.emplace_back(TemplatePosition{
            tpl[i - 1],
            params[0],  // match
            params[1],  // branch
            params[2],  // stick
            params[3]   // deletion
        });
    }

    result.emplace_back(TemplatePosition{tpl.back(), 1.0, 0.0, 0.0, 0.0});

    return result;
}

double SP1C1NoCovModel::BaseEmissionPr(MoveType move, const char from, const char to) const
{
    // All emissions are now "Covariate Emissions"
    assert(move != MoveType::DELETION);
    return 1.0;
}

double SP1C1NoCovModel::CovEmissionPr(MoveType move, uint8_t nuc, char from, char to) const
{
    const uint8_t hpAdd = from == to ? 0 : 4;

    to = detail::TranslationTable[static_cast<uint8_t>(to)];
    if (to > 3 || nuc > 3) throw std::invalid_argument("invalid character in sequence");

    // Which row do we want?
    const uint8_t row = to + hpAdd;

    if (move == MoveType::BRANCH)
        return branchPmf[row][nuc];
    else if (move == MoveType::STICK)
        return stickPmf[row][nuc];
    else if (move == MoveType::MATCH)
        return matchPmf[row][nuc];

    throw std::invalid_argument("unknown move type!");
}

double SP1C1NoCovModel::UndoCounterWeights(const size_t nEmissions) const { return 0; }
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
