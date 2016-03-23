
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

double matchPmf[8][4] = {
    {0.980417570, 0.011537479, 0.005804964, 0.002239987},  // AA
    {0.026122324, 0.972937583, 0.000367796, 0.000572296},  // CC
    {0.002544283, 0.002239375, 0.962042375, 0.033173967},  // GG
    {0.000509814, 0.001489097, 0.094228328, 0.903772761},  // TT
    {0.979840156, 0.012582917, 0.005185205, 0.002391722},  // NA
    {0.015528755, 0.984439781, 7.91000E-07, 3.07000E-05},  // NC
    {0.002667013, 0.002095727, 0.961571053, 0.033666207},  // NG
    {0.000506358, 0.001057035, 0.116124340, 0.882312267}   // NT
};

double stickPmf[8][4] = {
    {0.000000000, 0.254503401, 0.574809968, 0.170686631},  // AA
    {0.399446202, 0.000000000, 0.510664061, 0.089889737},  // CC
    {0.505214805, 0.188597323, 0.000000000, 0.306187872},  // GG
    {0.361855644, 0.132870306, 0.505274050, 0.000000000},  // TT
    {0.000000000, 0.210676350, 0.615161689, 0.174161960},  // NA
    {0.357451562, 0.000000000, 0.473482915, 0.169065523},  // NC
    {0.577147745, 0.169785817, 0.000000000, 0.253066438},  // NG
    {0.446834358, 0.144605809, 0.408559833, 0.000000000}   // NT
};

double branchPmf[8][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
                          {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

double SP1C1BetaNoCovParams[8][4] = {
    // Match, Branch, Stick, Delete
    {0.888913751, 0.021169653, 0.034937054, 0.054979542},  // AA
    {0.835822697, 0.036126801, 0.091992041, 0.036058461},  // CC
    {0.886427657, 0.022596867, 0.039619893, 0.051355584},  // GG
    {0.821252207, 0.072798639, 0.068161389, 0.037787765},  // TT
    {0.857630366, 0.072058988, 0.036435296, 0.033875351},  // NA
    {0.846000625, 0.032981179, 0.076759732, 0.044258463},  // NC
    {0.881462348, 0.042444137, 0.039293952, 0.036799562},  // NG
    {0.879087800, 0.022178294, 0.057073518, 0.041660389}   // NT
};

SP1C1BetaNoCovModel::SP1C1BetaNoCovModel(const SNR& snr) : snr_(snr) {}
std::vector<TemplatePosition> SP1C1BetaNoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t b = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];

        if (b > 3) throw std::invalid_argument("invalid character in sequence!");

        const bool hpAdd = tpl[i - 1] == tpl[i] ? 0 : 4;
        const auto params = SP1C1BetaNoCovParams[b + hpAdd];

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

double SP1C1BetaNoCovModel::BaseEmissionPr(MoveType move, const char from, const char to) const
{
    // All emissions are now "Covariate Emissions"
    assert(move != MoveType::DELETION);
    return 1.0;
}

double SP1C1BetaNoCovModel::CovEmissionPr(MoveType move, uint8_t nuc, char from, char to) const
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

double SP1C1BetaNoCovModel::UndoCounterWeights(const size_t nEmissions) const { return 0; }
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio