#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <limits>
//#include "LongDouble.hpp"
//#include "LongInt.hpp"
//#include "Matrix.hpp"

namespace math {

    bool Chance (uint8_t percent);

    bool Chance (double percent);

    uint64_t Factorial (uint64_t n);

    uint64_t P (uint64_t n);

    uint64_t C (uint64_t n, uint64_t k);

}

bool isEqual(double x, double y);

#endif