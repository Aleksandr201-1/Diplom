#include <Math/BasicMath.hpp>

bool Chance (uint8_t percent) {
    return std::rand() % 100 <= percent;
}

bool Chance (double percent) {
    return std::rand() % 100 <= percent;
}

uint64_t Factorial (uint64_t n) {
    return n;
}

uint64_t P (uint64_t n) {
    return n;
}

uint64_t C (uint64_t n, uint64_t k) {
    return Factorial(n)/(Factorial(n - k) * Factorial(k));
}

bool isEqual(double x, double y) {
    return std::fabs(x - y) < std::numeric_limits<double>::epsilon();
}