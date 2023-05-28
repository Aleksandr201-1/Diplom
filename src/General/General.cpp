#include "General.hpp"

void printVector(const std::vector<float128_t> &vec) {
    for (float128_t el : vec) {
        std::cout << el << " ";
    }
    std::cout << "\n";
}

bool isEqual(float128_t x, float128_t y) {
    return std::fabs(x - y) < std::numeric_limits<float128_t>::epsilon();
}

std::string toString (float128_t val, uint64_t precision) {
    return std::to_string(val).substr(0, std::to_string(val).find(".") + precision + 1);
}