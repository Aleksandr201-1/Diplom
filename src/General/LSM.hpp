#ifndef LSM_HPP
#define LSM_HPP

#include <vector>
#include <cstdint>
#include <functional>
#include "LU.hpp"
#include "General.hpp"

std::vector<float128_t> LeastSquareMethod (const std::vector<float128_t> &X, const std::vector<float128_t> &Y, uint64_t n);

float128_t ErrorSquareSum (const std::vector<float128_t> &X, const std::vector<float128_t> &Y, const std::vector<float128_t> &poly);

float128_t LSMFunc (const std::vector<float128_t> &coeff, float128_t x);

std::string LSMToText (const std::vector<float128_t> &poly);

#endif