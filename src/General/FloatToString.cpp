#include <General/FloatToString.hpp>

std::string toString (float val, uint64_t precision) {
    return std::to_string(val).substr(0, std::to_string(val).find(",") + precision + 1);
}