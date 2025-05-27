#ifndef LU_HPP
#define LU_HPP

#include <iostream>
#include <tuple>
#include <Math/Matrix.hpp>

template <class T>
Matrix<T> createP (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;
    Matrix<T> P(n), A(matrix);
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t toSwap = i;
        T max = std::abs(matrix(i, i));
        for (uint64_t j = i + 1; j < n; ++j) {
            if (std::abs(matrix(j, i)) > max) {
                max = std::abs(matrix(j, i));
                toSwap = j;
            }
        }
        if (max == T(0)) {
            exit(-1);
        }
        if (toSwap != i) {
            P.swapRows(i, toSwap);
            A.swapRows(i, toSwap);
        }
        for (uint64_t j = i + 1; j < n; ++j) {
            T tmp = P(j, i) / P(i, i);
            for (uint64_t k = i + 1; k < n; ++k) {
                A(j, k) -= tmp * A(i, k);
            }
        }
    }
    return P;
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> LU (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;
    Matrix<T> L(n), U(matrix), P(n);
    P = createP(matrix);
    U = P * U;
    for(uint64_t k = 1; k < n; ++k) {
        for(uint64_t i = k - 1; i < n; ++i) {
            for(uint64_t j = i; j < n; ++j) {
                L(j, i) = U(j, i) / U(i, i);
            }
        }
        for(uint64_t i = k; i < n; ++i) {
            for(uint64_t j = k - 1; j < n; ++j) {
                U(i, j) = U(i, j) - L(i, k - 1) * U(k - 1, j);
            }
        }
    }
    return std::make_tuple(L, U, P);
}

template <class T>
std::vector<T> LUsolveSLAE (const Matrix<T> &matrix, const std::vector<T> &ans) {
    uint64_t n = matrix.size().n;

    if (!matrix.isSquare() || n != ans.size()) {
        std::cerr << "Matrix is not square. Stop working.\n";
        exit(-1);
    }
    if (n != ans.size()) {
        std::cerr << "Matrix and vector have different sizes. Stop working.\n";
        exit(-1);
    }

    Matrix<T> A(matrix);
    std::vector<T> x(n, T(0)), y(n, T(0));
    Matrix<T> L, U, P;
    std::tie(L, U, P) = LU(A);
    std::vector<T> b = (P * Matrix<T>(n, 1, ans)).toVector();

    y[0] = b[0];
    for (uint64_t i = 1; i < n; ++i) {
        y[i] = b[i];
        for (uint64_t j = 0; j < i; ++j) {
            y[i] -= y[j] * L(i, j);
        }
    }

    x[n - 1] = y[n - 1] / U(n - 1, n - 1);
    for (uint64_t i = n - 2; i < n - 1; --i) {
        x[i] = y[i];
        for (uint64_t j = i + 1; j < n; ++j) {
            x[i] -= x[j] * U(i, j);
        }
        x[i] /= U(i, i);
    }
    return x;
}

template <class T>
Matrix<T> LUReverseMatrix (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;

    if (!matrix.isSquare()) {
        std::cerr << "Matrix is not square. Stop working.\n";
        exit(-1);
    }

    Matrix<T> A(matrix);
    Matrix<T> L, U, P;
    std::tie(L, U, P) = LU(A);

    Matrix<T> Ar(A), Lr(n), Ur(n);
    
    for (uint64_t i = 0; i < n; ++i) {
        Lr(i, i) = 1.0 / L(i, i);
        for (uint64_t j = i + 1; j < n; ++j) {
            T tmp = 0;
            for (uint64_t k = 0; k < j; ++k) {
                tmp -= L(j, k) * Lr(k, i);
            }
            Lr(j, i) = tmp / L(j, j);
        }
    }

    for (uint64_t i = n - 1; i < n; --i) {
        Ur(i, i) = 1.0 / U(i, i);
        for (uint64_t j = i - 1; j < n; --j) {
            T tmp = 0;
            for (uint64_t k = i; k > j; --k) {
                tmp -= U(j, k) * Ur(k, i);
            }
            Ur(j, i) = tmp / U(j, j);
        }
    }

    Ar = Ur * Lr * P;
    return Ar;
}

template <class T>
T LUGetDet (const Matrix<T> &matrix) {
    T ans = 1;
    uint64_t n = matrix.size().n;

    if (!matrix.isSquare()) {
        std::cerr << "Matrix is not square. Stop working.\n";
        exit(-1);
    }

    Matrix<T> A(matrix);
    Matrix<T> L, U, P;
    std::tie(L, U, P) = LU(A);
    for (uint64_t i = 0; i < U.size().n; ++i) {
        ans *= U(i, i);
    }

    return ans;
}

#endif