#ifndef QR_HPP
#define QR_HPP

#include <Math/Matrix.hpp>

const uint64_t ITERATION_CAP = 100;

template <class T>
T getMainEl (const Matrix<T> &matrix, uint64_t i) {
    T ans = matrix(i, i), tmp = 0;
    T sign = ans >= T(0) ? T(1) : T(-1);
    for (uint64_t j = i; j < matrix.size().n; ++j) {
        tmp += matrix(j, i) * matrix(j, i);
    }
    return ans + sign * std::sqrt(tmp);
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>> QR (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;
    Matrix<T> Q(n), R(matrix), v(n, 1);
    auto pred = [] (T el) -> bool {
        return el == T(0);
    };
    for (uint64_t i = 0; i < n - 1; ++i) {
        for (uint64_t j = 0; j < i; ++j) {
            v(j, 0) = 0;
        }
        v(i, 0) = getMainEl(R, i);
        for (uint64_t j = i + 1; j < n; ++j) {
            v(j, 0) = R(j, i);
        }

        if (std::count_if(v.toVector().cbegin(), v.toVector().cend(), pred) == n) {
            ++i;
            continue;
        }

        Matrix<T> H = Matrix<T>(n) - 2 * ((v * v.transp()) / (v.transp() * v)(0, 0));
        R = H * R;
        Q = Q * H;
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < i; ++j) {
            R(i, j) = T(0);
        }
    }
    return std::make_tuple(Q, R);
}

template <class T>
bool QRCheckForEnd (const Matrix<T> &matrix, T approx) {
    T sum = 0;
    for (uint64_t j = 1; j < matrix.size().n; ++j) {
        sum += matrix(j, 0) * matrix(j, 0);
    }
    return std::sqrt(sum) <= approx;
}

template <class T>
std::vector<T> QRFindLambda (const Matrix<T> &matrix, T approx) {
    std::vector<T> lambda;
    uint64_t iteration = 0, n = matrix.size().n;
    Matrix<T> Q, R, A(matrix);
    while (1) {
        ++iteration;
        std::tie(Q, R) = QR(A);
        A = R * Q;
        if (QRCheckForEnd(A, approx) || iteration > ITERATION_CAP) {
            break;
        }
    }
    for (uint64_t i = 0; i < n; ++i) {
        T sum = 0;
        for (uint64_t j = i + 1; j < A.size().n; ++j) {
            sum += A(j, i) * A(j, i);
        }
        sum = std::sqrt(sum);
        uint64_t j = i + 1;
        if (sum <= approx) {
            lambda.push_back(A(i, i));
        } else {
            T b = -(A(i, i) + A(j, j)), c = A(i, i) * A(j, j) - A(i, j) * A(j, i), a = T(1);
            T D = b * b - 4 * a * c;
            ++i;
            lambda.push_back(-b / 2);
            lambda.push_back(-b / 2);
        }
    }
    return lambda;
}

#endif