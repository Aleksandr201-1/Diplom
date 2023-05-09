#ifndef TOUGH_DET
#define TOUGH_DET

#include <algorithm>
#include <General/General.hpp>
#include <General/Differentiation.hpp>
#include <General/LU.hpp>
#include <General/QR.hpp>

float128_t getToughK (const Matrix<float128_t> &K);

float128_t ToughCoeff (const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &system, const std::vector<float128_t> &Y);

#endif