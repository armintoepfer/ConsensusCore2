
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
        {0.974438503120005, 0.0231154641889822, 0.000693772904813349, 0.00175225978619963},
        {0.0148815232512677, 0.985074652083733, 2.12798129933467e-05, 2.25448520055512e-05},
        {0.00142617067001849, 0.0011709460846, 0.957485390066542, 0.03991749317884},
        {0.000190779168842439, 0.000422847545261742, 0.010793657731581, 0.988592715554315},
        {0.978532977779298, 0.0194094896214458, 0.00050974562200914, 0.0015477869772466},
        {0.00966654187108636, 0.990330387549539, 1.17342392883864e-07, 2.95323698169498e-06},
        {0.000313697473099085, 0.000355906620188189, 0.934285555577205, 0.0650448403295073},
        {7.28271610486195e-05, 0.0003190393805037, 0.0121228940116563, 0.987485239446791},
    };

    double stickPmf[8][4] = {
        {0, 0.13865550950802, 0.128822123792677, 0.732522366699303},
        {0.115173737905619, 0, 0.223763681005227, 0.661062581089154},
        {0.232798551855526, 0.13461704101661, 0, 0.632584407127864},
        {0.33491264127151, 0.195749138118978, 0.469338220609513, 0},
        {0, 0.187261768439734, 0.180246309606338, 0.632491921953929},
        {0.125678953429545, 0, 0.191804013971638, 0.682517032598816},
        {0.171461555636964, 0.0975160790844003, 0, 0.731022365278636},
        {0.670837572326349, 0.201825702245386, 0.127336725428265, 0},
    };



double branchPmf[8][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
                          {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

double SP1C1NoCovParams[8][4] = {
    // Match, Branch, Stick, Delete
    {0.917186930328304,0.00937752453744426,0.0236490481257219,0.04978649700852},
    {0.869775437060557,0.0230955005212403,0.0731701023878155,0.033958960030387},
    {0.91416794226825,0.0131641580599221,0.0188020174975751,0.0538658821742525},
    {0.851460502647672,0.0979579973331125,0.0200978805107598,0.030483619508455},
    {0.880432572806988,0.0629195151878144,0.0280425687131194,0.028605343292078},
    {0.86792462348301,0.0431375425589477,0.0556305091764776,0.0333073247815644},
    {0.920039347845112,0.0268832239875899,0.0210383493712407,0.032039078796057},
    {0.927058749847394,0.023467357278295,0.0133371845336218,0.0361367083406894}


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
