
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
        
        {    0.989796824,   0.00611512367,   0.00361139203,  0.000317110233},
        {   0.0438468606,     0.955861663,  3.60892703e-05,  7.46227591e-05},
        {  0.00210254989,  0.000573693076,      0.98359225,    0.0135291219},
        {   0.0057021146,  0.000480406692,     0.134181337,     0.859468123},
        {    0.989622946,   0.00550951163,   0.00438840662,  0.000297588711},
        {   0.0432620553,     0.956477372,  3.78020659e-05,  5.27299528e-05},
        {  0.00191400896,  0.000767143289,     0.988147124,   0.00902657585},
        {  0.00415625613,   0.00650760837,      0.14358637,     0.845547995},
        {    0.996181635,   0.00343644301,  0.000127966073,  6.58048632e-05},
        {   0.0297878632,      0.97000734,  1.83731345e-05,  3.88397462e-05},
        {  0.00132233598,   0.00062328985,     0.988619216,   0.00926046657},
        {   0.0016622058,   0.00104351365,    0.0903609731,     0.906740841},
        {    0.993556712,   0.00534815673,   0.00087340051,  4.68034345e-05},
        {   0.0217175687,     0.978047742,  2.12689235e-05,  2.47942916e-05},
        { 0.000166772601,  0.000144323764,     0.991363392,   0.00813946371},
        {   0.0414476652,    0.0326367258,     0.116440083,     0.809328411}
    };
    
    double stickPmf[16][4]  = {
        
        {  0.00057849858,     0.168095833,     0.725067537,    0.0987376491},
        {    0.325013518,  0.000279910758,     0.556546155,     0.114521576},
        {     0.41261161,     0.139970701,   0.00200578411,     0.419336711},
        {    0.429212402,    0.0935151135,      0.46817297,   0.00064996537},
        { 0.000328593333,     0.293033832,     0.634411483,    0.0679543788},
        {    0.264359181,  0.000163201626,     0.667418249,    0.0660763117},
        {    0.499794555,     0.412440234,  0.000345204586,    0.0829323469},
        {    0.302876479,     0.361833317,     0.330623106,  0.000333364109},
        {  0.00149158764,     0.224621451,     0.469310556,     0.285185766},
        {    0.217956961,  0.000180642959,     0.677569059,     0.101944978},
        {    0.727547304,    0.0603088806,   0.00100203319,     0.198597329},
        {    0.799264631,     0.105271694,    0.0810926815,    0.0010264995},
        { 0.000370350081,    0.0965439021,     0.564294941,     0.333976256},
        {    0.206735425,  0.000216962947,      0.43235045,     0.357876644},
        {    0.400548963,    0.0428663359,  0.000506914437,     0.549487899},
        {    0.280098325,    0.0152265452,     0.699888104,  0.000379997439}
    };
    
    double branchPmf[16][4]  = {
        
        {    0.981808679,   0.00113695758,   0.00113695758,   0.00113695758},
        {  0.00302928789,     0.951531394,   0.00302928789,   0.00302928789},
        { 0.000410541908,  0.000410541908,     0.993431329,  0.000410541908},
        {  0.00812272058,   0.00812272058,   0.00812272058,     0.870036471},
        {    0.997116005,  0.000180249708,  0.000180249708,  0.000180249708},
        { 0.000347500053,     0.994459945,   0.00034617035,   0.00034617035},
        { 0.000241907509,  0.000241907509,      0.99612948,  0.000241907509},
        {  0.00274484571,   0.00274484571,   0.00274484571,     0.956082469},
        {    0.996924906,  0.000192193357,  0.000192193357,  0.000192193357},
        { 0.000625434729,     0.989993044,  0.000625434729,  0.000625434729},
        {  0.00117473877,   0.00117113559,     0.981258227,   0.00117113559},
        { 0.000839138215,  0.000839138215,  0.000839138215,     0.986573789},
        {    0.997587399,  0.000150787571,  0.000150787571,  0.000150787571},
        { 0.000367034374,      0.99412745,  0.000367034374,  0.000367034374},
        { 0.000212909576,  0.000212909576,     0.996593447,  0.000212909576},
        { 0.000164104279,  0.000160679593,  0.000160679593,     0.997425702}
    };
    
    double transProbs[16][4] = {
        
        {    0.916078871,   0.00948455354,    0.0186910959,    0.0557454795}, // CTX - 0
        {    0.907100404,   0.00424720071,    0.0456251936,    0.0430272018}, // CTX - 1
        {     0.91350938,    0.0341923388,   0.00678834576,    0.0455099352}, // CTX - 2
        {    0.931421805,   0.00149577949,     0.017888726,    0.0491936896}, // CTX - 3
        {    0.874104337,    0.0662682259,    0.0360495689,     0.023577868}, // CTX - 4
        {    0.860029784,    0.0317095319,    0.0721907339,    0.0360699503}, // CTX - 5
        {    0.909532024,    0.0416719746,    0.0287524926,    0.0200435091}, // CTX - 6
        {    0.931046123,    0.0052228796,    0.0415013302,     0.022229667}, // CTX - 7
        {    0.879792849,    0.0669311253,   0.00832333078,    0.0449526951}, // CTX - 8
        {    0.882103577,    0.0165557736,    0.0565802988,    0.0447603503}, // CTX - 9
        {    0.927122924,    0.0100329917,    0.0122308441,    0.0506132405}, // CTX - 10
        {    0.921584598,     0.016066281,    0.0127660113,    0.0495831094}, // CTX - 11
        {    0.822618676,    0.0745595178,    0.0298228614,    0.0729989448}, // CTX - 12
        {    0.837367808,    0.0347371734,    0.0581666692,    0.0697283494}, // CTX - 13
        {    0.851139418,    0.0579469279,    0.0237841934,    0.0671294602}, // CTX - 14
        {    0.847455166,    0.0633698287,    0.0300045522,    0.0591704529} // CTX - 15
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
