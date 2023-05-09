#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <cmath>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <functional>
#include <iostream>

using float128_t = long double;

const uint64_t ITERATION_CAP = 200;
const uint64_t PRECISION = 2;

void printVector(const std::vector<float128_t> &vec);

bool isEqual(float128_t x, float128_t y);

std::string toString (float128_t val, uint64_t precision);

std::string readLine ();

#endif