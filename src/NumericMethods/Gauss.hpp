#ifndef GAUSS_HPP
#define GAUSS_HPP

#include <iostream>
#include <Math/Matrix.hpp>

template<class T>
std::vector<T> GaussSolveSLAE (const Matrix<T> &matrix, const std::vector<T> &ans, std::string ur = "Q") {
    Matrix<T> A = matrix;
    std::vector<T> b(ans), x(ans), used(ans.size(), 0);
    uint64_t n = ans.size();
    for (uint64_t i = 0; i < n; ++i) {
        T maxAbsEl = 0;
        uint64_t idx = -1;
        for (uint64_t j = 0; j < n; ++j) {
            if (std::abs(maxAbsEl) < std::abs(A(j, i)) && used[j] == 0) {
                maxAbsEl = A(j, i);
                idx = j;
            }
            //maxAbsEl = std::max(maxAbsEl, std::abs(A(j, i)));
        }
        if (idx == -1) {
            std::cout << "solvin " << ur << "\n";
            std::cout << "matrix M:\n" << matrix << "\n";
            std::cout << "vector ans:\n";
            for (auto el : ans) {
                std::cout << el << " ";
            }
            std::cout << "\nmatrix A:\n" << A << "\n";
            std::cout << "Undefined x\n";
            exit(-1);
        }
        used[idx] = 1;
        for  (uint64_t j = 0; j < n; ++j) {
            if (A(j, i) == 0) {
                continue;
            }
            if (j != idx) {
                T div = A(j, i) / maxAbsEl;
                b[j] -= div * b[idx];
                for (uint64_t k = i; k < n; ++k) {
                    A(j, k) -= div * A(idx, k);
                }
                A(j, i) = 0.0;
            }
        }
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            if (A(i, j) != 0.0) {
                x[j] = b[i] / A(i, j);
                break;
            }
        }
    }
    return x;
}

#endif