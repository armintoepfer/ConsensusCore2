
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
        {0.972986946169289, 0.0233711777162974, 0.000805743669989392, 0.00283613244442378},
        {0.0099294414216902, 0.989745799324021, 3.14907056417152e-05, 0.000293268548647062},
        {0.00137849639827314, 0.001313693526526, 0.953997270203777, 0.0433105398714241},
        {0.000272101671961153, 0.000714773542074764, 0.0122733319882684, 0.986739792797696},
        {0.98132448630424, 0.0168971584108813, 0.000692268865653737, 0.00108608641922515},
        {0.00560128678171135, 0.994397607897536, 3.09723282589416e-07, 7.9559746973039e-07},
        {0.00111157341433363, 0.000502012406163902, 0.926897626647946, 0.0714887875315565},
        {0.000128690228265445, 0.000422843432896165, 0.0126629514075207, 0.986785514931318},
    };


    double stickPmf[8][4] = {
        {0, 0.239146766700572, 0.13213945162201, 0.628713781677418},
        {0.129983237503601, 0, 0.244146935598477, 0.625869826897922},
        {0.263482633722373, 0.156631189102475, 0, 0.579886177175153},
        {0.344480827656803, 0.230468665691304, 0.425050506651892, 0},
        {0, 0.203815845360318, 0.170777571876277, 0.625406582763405},
        {0.144455045732362, 0, 0.18095008111473, 0.674594873152908},
        {0.187399486860146, 0.127010631799708, 0, 0.685589881340147},
        {0.558916082338799, 0.258021588012282, 0.183062329648919, 0}
    };

double branchPmf[8][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
                          {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

double SP1C1NoCovParams[8][4] = {
    // Match, Branch, Stick, Delete
        {0.903873647556401,0.012593793538233,0.0259571817445857,0.0575753771607807},
        {0.85930483737831,0.0243619684256844,0.0794157748051166,0.0369174193908887},
        {0.906826791642645,0.0147882703954811,0.0223633242358615,0.056021613726012},
        {0.854674584716752,0.0858037444633542,0.0252595722665351,0.034262098553358},
        {0.879258662581103,0.0623101287993231,0.0276450467294577,0.030786161890116},
        {0.859419129376976,0.0439876262048763,0.0600861141449959,0.036507130273152},
        {0.913148356418287,0.0249760986728176,0.0250307174761157,0.036844827432779},
        {0.915945207273883,0.0248427165496227,0.0178499098112664,0.041362166365228}
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
