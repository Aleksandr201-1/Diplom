#ifndef ITERATION_ALGO_HPP
#define ITERATION_ALGO_HPP

#include <vector>
#include <functional>
#include <map>
#include <algorithm>
#include <cmath>
#include <General/Enum.hpp>
#include <Math/Matrix.hpp>
#include <NumericMethods/LU.hpp>

enum class IterationAlgo {
    NEWTON,
    ZEIDEL,
    SIMPLE_ITERATION,
    EXPLICIT_STEP,
    ERROR
};

std::string IterationAlgoToString (IterationAlgo algo);

IterationAlgo stringToIterationAlgo (const std::string &str);

double norma (const std::vector<double> &a, const std::vector<double> &b);

std::vector<std::vector<double>> ExplicitStep (const std::vector<std::vector<double>> &K,
                                                   const std::vector<std::vector<double>> &Yi,
                                                   uint64_t idx,
                                                   const std::vector<std::function<double (const std::vector<double> &)>> &f,
                                                   const Matrix<double> &butcher,
                                                   double h);

std::vector<std::vector<double>> SimpleIteration (const std::vector<std::vector<double>> &K,
                                                      const std::vector<std::vector<double>> &Yi,
                                                      uint64_t idx,
                                                      const std::vector<std::function<double (const std::vector<double> &)>> &f,
                                                      const Matrix<double> &butcher,
                                                      double h,
                                                      double approx,
                                                      IterationAlgo algo);

std::vector<std::vector<double>> NewtonIteration (const std::vector<std::vector<double>> &K,
                                                      const std::vector<std::vector<double>> &Yi,
                                                      uint64_t idx,
                                                      const std::vector<std::function<double (const std::vector<double> &)>> &f,
                                                      const Matrix<double> &butcher,
                                                      double h,
                                                      double approx);

std::vector<std::vector<double>> IterationStep (const std::vector<std::vector<double>> &K,
                                                    const std::vector<std::vector<double>> &Yi,
                                                    uint64_t idx,
                                                    const std::vector<std::function<double (const std::vector<double> &)>> &f,
                                                    const Matrix<double> &butcher,
                                                    double h,
                                                    double approx,
                                                    IterationAlgo algo);

#endif