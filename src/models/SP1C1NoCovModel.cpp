
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
        {0.974836771831626, 0.0227950722439149, 0.000792599316271079, 0.00157555660818828},
        {0.00868509499422849, 0.991265182178836, 3.57187589910245e-06, 4.6150951036181e-05},
        {0.00113530059131767, 0.000775090718161884, 0.954964392606382, 0.0431252160841386},
        {0.000256922021690147, 0.000666063636620814, 0.0136937633613173, 0.985383250980372},
        {0.980810763268474, 0.0174954672194116, 0.000669728326881354, 0.00102404118523298},
        {0.00521591634970498, 0.994783109692156, 1.07896887997981e-07, 8.66061250671744e-07},
        {0.00108754853223241, 0.000584025150937109, 0.927240779840482, 0.0710876464763484},
        {0.000115664594236018, 0.000990384634800216, 0.0133056320233962, 0.985588318747567},
    };



    double stickPmf[8][4] = {
        {0, 0.213894084340011, 0.138955823614205, 0.647150092045784},
        {0.142489297188992, 0, 0.216736456276613, 0.640774246534395},
        {0.206410352721341, 0.162804765410477, 0, 0.630784881868183},
        {0.337773459741083, 0.195316463293848, 0.46691007696507, 0},
        {0, 0.22127339467203, 0.144703838261638, 0.634022767066333},
        {0.140923782727102, 0, 0.175722930041557, 0.683353287231341},
        {0.170567249115115, 0.107033493273684, 0, 0.722399257611201},
        {0.585600549923167, 0.25003507340633, 0.164364376670503, 0},
    };


double branchPmf[8][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
                          {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

double SP1C1NoCovParams[8][4] = {
    // Match, Branch, Stick, Delete
    {0.909527410244176,0.00996452515885327,0.0223042816193927,0.0582037829775},
    {0.86058433075905,0.026718785093459,0.0765871901853612,0.0361096939621294},
    {0.898783700185153,0.0243686461784536,0.0226696250130183,0.05417802862337},
    {0.845716679347856,0.0979251820231239,0.0226103875909803,0.03374775103803},
    {0.879785703328608,0.0637053619587933,0.0267899642293115,0.02971897048328},
    {0.860686673972841,0.0463179944820158,0.0557921041242302,0.03720322742091},
    {0.913650705259057,0.0268872628702281,0.0234743423338464,0.03598768953686},
    {0.920297132397726,0.0246397330289372,0.0155097857857348,0.03955334878760}

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
