#ifndef LSM_HPP
#define LSM_HPP

#include <vector>
#include <cstdint>
#include <cmath>
#include <functional>
#include <NumericMethods/LU.hpp>
#include <General/FloatToString.hpp>

template <typename T>
std::vector<T> LeastSquareMethod (const std::vector<T> &X, const std::vector<T> &Y, uint64_t n);

template <typename T>
T ErrorSquareSum (const std::vector<T> &X, const std::vector<T> &Y, const std::vector<T> &poly);

template <typename T>
T LSMFunc (const std::vector<T> &coeff, T x);

template <typename T>
std::string LSMToText (const std::vector<T> &poly, uint64_t precision);

#include <NumericMethods/LSM.tpp>

#endif