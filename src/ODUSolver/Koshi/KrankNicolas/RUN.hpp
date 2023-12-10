#ifndef RUN_HPP
#define RUN_HPP

#include "General.hpp"
#include "Matrix.hpp"

//проверка матрицы на 3х-диагональность
template <class T>
bool checkIf3Diagonal (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;
    T null = T(0);
    bool result = true;
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t left, right;
        left = i - 1 > n ? 0 : i - 1;
        right = i + 2 > n ? n : i + 2;
        for (uint64_t j = 0; j < left; ++j) {
            if (!isEqual(matrix(i, j), null)) {
                result = false;
                break;
            }
        }
        for (uint64_t j = left; j < right; ++j) {
            if (isEqual(matrix(i, j), null)) {
                result = false;
                break;
            }
        }
        for (uint64_t j = right; j < n; ++j) {
            if (!isEqual(matrix(i, j), null)) {
                result = false;
                break;
            }
        }
    }
    return result;
}

//метод прогонки для решения СЛАУ
template <class T>
std::vector<T> RUNsolveSLAE (const Matrix<T> &matrix, const std::vector<T> &ans) {
    uint64_t n = matrix.size().n;
    std::vector<T> P(n), Q(n), x(n);

    if (!matrix.isSquare() || n != ans.size()) {
        return {};
    }
    if (n != ans.size()) {
        return {};
    }
    if (!checkIf3Diagonal(matrix)) {
        return {};
    }

    T a = 0, b = matrix(0, 0), c = matrix(0, 1), d = ans[0];
    P[0] = -c / b;
    Q[0] = d / b;
    for (uint64_t i = 1; i < n; ++i) {
        a = matrix(i, i - 1);
        b = matrix(i, i);
        c = i + 1 < n ? matrix(i, i + 1) : 0;
        d = ans[i];
        P[i] = -c / (b + a * P[i - 1]);
        Q[i] = (d - a * Q[i - 1]) / (b + a * P[i - 1]);
    }
    x[n - 1] = Q[n - 1];
    for (uint64_t i = n - 2; i < n; --i) {
        x[i] = P[i] * x[i + 1] + Q[i];
    }
    return x;
}

#endif