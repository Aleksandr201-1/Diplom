#include "LSM.hpp"

std::vector<double> LeastSquareMethod (const std::vector<double> &X, const std::vector<double> &Y, uint64_t n) {
    Matrix<double> M(n + 1, n + 1);
    std::vector<double> y(n + 1);
    for (uint64_t i = 0; i < n + 1; ++i) {
        double tmp = 0;
        for (uint64_t j = 0; j < X.size(); ++j) {
            tmp += std::pow(X[j], i);
        }
        M(0, i) = tmp;
    }
    for (uint64_t i = 1; i < n + 1; ++i) {
        double tmp = 0;
        for (uint64_t j = 0; j < n; ++j) {
            M(i, j) = M(i - 1, j + 1);
        }
        for (uint64_t j = 0; j < X.size(); ++j) {
            tmp += std::pow(X[j], n + i);
        }
        M(i, n) = tmp;
    }
    for (uint64_t i = 0; i < y.size(); ++i) {
        double tmp = 0;
        for (uint64_t j = 0; j < Y.size(); ++j) {
            tmp += Y[j] * std::pow(X[j], i);
        }
        y[i] = tmp;
    }
    return LUsolveSLAE(M, y);
}

double ErrorSquareSum (const std::vector<double> &X, const std::vector<double> &Y, const std::vector<double> &poly) {
    double ans = 0;
    for (uint64_t i = 0; i < X.size(); ++i) {
        double tmp = 0;
        for (uint64_t j = 0; j < poly.size(); ++j) {
            tmp += poly[j] * std::pow(X[i], j);
        }
        ans += (tmp - Y[i]) * (tmp - Y[i]);
    }
    return ans;
}

double LSMFunc (const std::vector<double> &coeff, double x) {
    double ans = 0;
    for (uint64_t i = 0; i < coeff.size(); ++i) {
        ans += coeff[i] * std::pow(x, i);
    }
    return ans;
}

std::string LSMToText (const std::vector<double> &poly) {
    std::string func = std::to_string(poly[0]);
    for (uint64_t i = 1; i < poly.size(); ++i) {
        if (poly[i] < 0) {
            func += " - ";
        } else {
            func += " + ";
        }
        func += toString(std::abs(poly[i]), PRECISION) + "x^" + std::to_string(i);
    }
    return func;
}