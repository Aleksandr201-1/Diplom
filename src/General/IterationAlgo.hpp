#ifndef ITERATION_ALGO_HPP
#define ITERATION_ALGO_HPP

#include <vector>
#include <functional>
#include "Matrix.hpp"
#include "LU.hpp"

enum class IterationAlgo {
    NEWTON,
    ZEIDEL,
    SIMPLE_ITERATION
};

double norma (const std::vector<double> &a, const std::vector<double> &b);

std::vector<std::vector<double>> ExplicitStep (const std::vector<std::vector<double>> &K, const std::vector<double> &X, const std::vector<std::vector<double>> &Yi, uint64_t i, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx);

std::vector<std::vector<double>> SimpleIteration (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx);

std::vector<std::vector<double>> Zeidel (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx);

std::vector<std::vector<double>> Newton (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx);

std::vector<std::vector<double>> IterationStep (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx, IterationAlgo algo);

#endif