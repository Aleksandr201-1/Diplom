#ifndef DIFFERENTIATION_HPP
#define DIFFERENTIATION_HPP

#include <map>
#include "General.hpp"

enum class DiffConfig {
    POINTS2_ORDER1_WAY1,
    POINTS3_ORDER1_WAY1,
    POINTS4_ORDER1_WAY1,
    POINTS4_ORDER1_WAY2
};

// float128_t derivative1point2 (const std::function<float128_t(std::vector<float128_t> &)> &f, const std::vector<float128_t> &X, float128_t h, uint64_t idx);

// float128_t derivative1point3 (const std::function<float128_t(std::vector<float128_t> &)> &f, const std::vector<float128_t> &X, float128_t h, uint64_t idx);

// float128_t derivative1point4 (const std::function<float128_t(std::vector<float128_t> &)> &f, const std::vector<float128_t> &X, float128_t h, uint64_t idx);

// float128_t derivative1point5 (const std::function<float128_t(std::vector<float128_t> &)> &f, const std::vector<float128_t> &X, float128_t h, uint64_t idx);

float128_t derivative (const std::function<float128_t(std::vector<float128_t> &)> &f,
                       const std::vector<float128_t> &X,
                       float128_t h,
                       uint64_t idx,
                       DiffConfig diff);

float128_t derivative (const std::function<float128_t(float128_t)> &f,
                       float128_t x,
                       float128_t h,
                       DiffConfig diff);

#endif