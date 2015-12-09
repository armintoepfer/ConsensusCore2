
#pragma once

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace PacBio {
namespace Consensus {

class AlphaBetaMismatch : public std::runtime_error
{
public:
    AlphaBetaMismatch(float v) : std::runtime_error(message(v)) {}
private:
    static std::string message(float v) {
        std::ostringstream ss;
        ss << "Alpha/Beta Mismatch (" << std::setprecision(5) << v << ')';
        return ss.str();
    }
};

class ChemistryNotFound : public std::runtime_error
{
public:
    ChemistryNotFound(const std::string& name)
        : std::runtime_error(std::string("chemistry not found: '") + name + "'")
    {
    }
};

}  // namespace Consensus
}  // namespace PacBio
