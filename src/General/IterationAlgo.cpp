#include "IterationAlgo.hpp"
#include <iostream>
double norma (const std::vector<double> &a, const std::vector<double> &b) {
    uint64_t size = std::min(a.size(), b.size());
    double ans = 0;
    for (uint64_t i = 0; i < size; ++i) {
        ans = std::max(ans, std::abs(a[i] - b[i]));
    }
    return ans;
}

double derivative (const std::function<double (const std::vector<double> &, const std::vector<std::vector<double>> &, const Matrix<double> &)> &f, const std::vector<double> &X, const std::vector<std::vector<double>> &K, const Matrix<double> &butcher, double h, uint64_t idx, uint64_t idy) {
    //std::vector<double> args(X);
    //Matrix<double> tmp(butcher);
    std::vector<std::vector<double>> Ktmp(K);
    double a, b, c, d;
    //args[idx] += 2 * h;
    Ktmp[idx][idy] += 2 * h;
    a = f(X, Ktmp, butcher);
    Ktmp[idx][idy] -= h;
    b = f(X, Ktmp, butcher);
    Ktmp[idx][idy] -= 2 * h;
    c = f(X, Ktmp, butcher);
    Ktmp[idx][idy] -= h;
    d = f(X, Ktmp, butcher);
    return (-a + 8 * b - 8 * c + d) / 12 * h;
}

std::vector<std::vector<double>> SimpleIteration (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, IterationAlgo algo, double approx) {
    std::vector<std::vector<double>> prev = X, next = X;
    uint64_t orderOfTask = f.size(), orderOfApprox = butcher.size().m - 1;
    std::vector<double> new_args(args);
    double diff;
    uint64_t iter = 0;
    do {
        diff = 0;
        prev = next;
        for (uint64_t j = 0; j < orderOfApprox; ++j) {
            new_args = args;
            //аргументы для функций
            new_args[0] += butcher(j, 0) * h;
            for (uint64_t k = 1; k < args.size() - 1; ++k) {
                for (uint64_t l = 0; l < orderOfApprox; ++l) {
                    if (algo == IterationAlgo::SIMPLE_ITERATION) {
                        new_args[k] += butcher(j, l + 1) * prev[k - 1][l] * h;
                    } else {
                        new_args[k] += butcher(j, l + 1) * next[k - 1][l] * h;
                    }
                }
            }

            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask; ++k) {
                next[k][j] = f[k](new_args);
            }
        }
        for (uint64_t i = 0; i < prev.size(); ++i) {
            diff = std::max(diff, norma(prev[i], next[i]));
        }
        std::cout << "diff: " << diff << "\n";
        ++iter;
    } while (diff > approx);
    return next;
}

std::vector<std::vector<double>> Newton (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx) {
    std::vector<std::vector<double>> prev = X, next = X;
    uint64_t orderOfTask = f.size(), orderOfApprox = butcher.size().m - 1;
    std::vector<double> new_args(args);
    double diff;
    uint64_t iter = 0;
    do {
        diff = 0;
        prev = next;
        for (uint64_t i = 0; i < orderOfApprox; ++i) {
            //аргументы для функций
            new_args = args;
            new_args[0] += butcher(i, 0) * h;
            for (uint64_t k = 1; k < args.size() - 1; ++k) {
                for (uint64_t l = 0; l < orderOfApprox; ++l) {
                    new_args[k] += butcher(i, l + 1) * next[k - 1][l] * h;
                }
            }
            //F
            std::vector<std::function<double (const std::vector<double> &, const std::vector<std::vector<double>> &, const Matrix<double> &)>> F(X.size());
            for (uint64_t j = 0; j < X.size(); ++j) {
                //args.size = orderOfTask * orderOfApprox + orderOfTask + 2
                auto lambda = [=] (const std::vector<double> &args, const std::vector<std::vector<double>> &K, const Matrix<double> &butcher) -> double {
                    double ans = butcher(j, i);
                    return ans;
                };
                F[j] = lambda;
            }
            Matrix<double> Yakobi(orderOfApprox);
            for (uint64_t j = 0; j < orderOfApprox; ++j) {
                for (uint64_t k = 0; k < orderOfApprox; ++k) {
                    Yakobi(j, k) = derivative(F[j], new_args, butcher, h, j, k);
                }
            }
            Yakobi = LUReverseMatrix(Yakobi);
            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask; ++k) {
                next[k][j] = f[k](new_args);
            }
        }
        for (uint64_t i = 0; i < prev.size(); ++i) {
            diff = std::max(diff, norma(prev[i], next[i]));
        }
        ++iter;
    } while (diff > approx);
    return next;
}

std::vector<std::vector<double>> Iteration (const std::vector<std::vector<double>> &X, const std::vector<double> &args, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx, IterationAlgo algo) {
    switch (algo) {
        case IterationAlgo::ZEIDEL:
        case IterationAlgo::SIMPLE_ITERATION:
            return SimpleIteration(X, args, f, butcher, h, algo, approx);
        case IterationAlgo::NEWTON:
            return Newton(X, args, f, butcher, h, approx);
    }
    return {{}};
}