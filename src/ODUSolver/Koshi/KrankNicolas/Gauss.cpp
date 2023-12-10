#include <iostream>
#include <cmath>
#include <limits>
#include "Matrix.hpp"

bool isEqual(double x, double y) {
    return std::fabs(x - y) < std::numeric_limits<double>::epsilon();
}

std::vector<double> GaussSolveSLAE (const Matrix<double> &matrix, const std::vector<double> &ans) {
    Matrix<double> A = matrix;
    std::vector<double> b(ans), x(ans), used(ans.size(), 0);
    uint64_t n = ans.size();
    for (uint64_t i = 0; i < n; ++i) {
        double maxAbsEl = 0;
        uint64_t idx = -1;
        for (uint64_t j = 0; j < n; ++j) {
            if (std::abs(maxAbsEl) < std::abs(A(j, i)) && used[j] == 0) {
                maxAbsEl = A(j, i);
                idx = j;
            }
            //maxAbsEl = std::max(maxAbsEl, std::abs(A(j, i)));
        }
        if (idx == -1) {
            std::cout << "Undefined x\n";
            exit(-1);
        }
        used[idx] = 1;
        for  (uint64_t j = 0; j < n; ++j) {
            if (A(j, i) == 0) {
                continue;
            }
            if (j != idx) {
                double div = A(j, i) / maxAbsEl;
                b[j] -= div * b[idx];
                for (uint64_t k = i; k < n; ++k) {
                    A(j, k) -= div * A(idx, k);
                }
                A(j, i) = 0.0;
                // for (uint64_t k = 0; k < n; ++k) {
                //     if (k != idx) {
                //         A(k, i) = 0;
                //     }
                // }
            }
        }
        //std::cout << A << '\n';
        // for  (uint64_t j = idx + 1; j < n; ++j) {
        //     double div = A(j, i) / maxAbsEl;
        //     x[j] -= div * maxAbsEl;
        //     for (uint64_t k = i; k < n; ++k) {
        //         A(j, k) -= div * maxAbsEl;
        //     }
        // }
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            if (A(i, j) != 0.0) {
                x[j] = b[i] / A(i, j);
                break;
            }
        }
    }
    // for (uint64_t i = 0; i < b.size(); ++i) {
    //     b[i] /= A(i, i);
    // }
    return x;
}

int main() {
    uint64_t n;
    std::cin >> n;
    Matrix<double> matrix(n);
    std::vector<double> ans(n);
    std::cin >> matrix;
    for (uint64_t i = 0; i < ans.size(); ++i) {
        std::cin >> ans[i];
    }
    ans = GaussSolveSLAE(matrix, ans);
    for (auto el : ans) {
        std::cout << el << ' ';
    }
    std::cout << '\n';
    return 0;
}
