#ifndef ITERATION_ALGO_HPP
#define ITERATION_ALGO_HPP

#include <vector>
#include <functional>
#include <map>
#include <algorithm>
#include <General/General.hpp>
#include <General/Enum.hpp>
#include <General/Matrix.hpp>
#include <General/LU.hpp>

enum class IterationAlgo {
    NEWTON,
    ZEIDEL,
    SIMPLE_ITERATION,
    EXPLICIT_STEP,
    ERROR
};

std::string IterationAlgoToString (IterationAlgo algo);

IterationAlgo stringToIterationAlgo (const std::string &str);

float128_t norma (const std::vector<float128_t> &a, const std::vector<float128_t> &b);

std::vector<std::vector<float128_t>> ExplicitStep (const std::vector<std::vector<float128_t>> &K,
                                                   const std::vector<std::vector<float128_t>> &Yi,
                                                   uint64_t idx,
                                                   const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f,
                                                   const Matrix<float128_t> &butcher,
                                                   float128_t h);

std::vector<std::vector<float128_t>> SimpleIteration (const std::vector<std::vector<float128_t>> &K,
                                                      const std::vector<std::vector<float128_t>> &Yi,
                                                      uint64_t idx,
                                                      const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f,
                                                      const Matrix<float128_t> &butcher,
                                                      float128_t h,
                                                      float128_t approx,
                                                      IterationAlgo algo);

std::vector<std::vector<float128_t>> NewtonIteration (const std::vector<std::vector<float128_t>> &K,
                                                      const std::vector<std::vector<float128_t>> &Yi,
                                                      uint64_t idx,
                                                      const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f,
                                                      const Matrix<float128_t> &butcher,
                                                      float128_t h,
                                                      float128_t approx);

std::vector<std::vector<float128_t>> IterationStep (const std::vector<std::vector<float128_t>> &K,
                                                    const std::vector<std::vector<float128_t>> &Yi,
                                                    uint64_t idx,
                                                    const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f,
                                                    const Matrix<float128_t> &butcher,
                                                    float128_t h,
                                                    float128_t approx,
                                                    IterationAlgo algo);

#endif