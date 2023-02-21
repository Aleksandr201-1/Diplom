#ifndef QR_HPP
#define QR_HPP

#include "General.hpp"
#include "Matrix.hpp"

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
    //std::cout << "====CREATING QR DECOMPOSITION===\n";
    //std::cout << "Matrix A:\n" << matrix << "\n";
    uint64_t n = matrix.size().n;
    Matrix<T> Q(n), R(matrix), v(n, 1);
    //std::cout << "------------\n";
    //std::cout << "Step 0:\nA(0) = A\n";
    auto pred = [] (T el) -> bool {
        return el == T(0);
    };
    for (uint64_t i = 0; i < n - 1; ++i) {
        //std::cout << "------------\n";
        //std::cout << "Step " << i + 1 << ":\n";

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

        //std::cout << "H = E - 2 * (v * vT / vT * v)\n";
        Matrix<T> H = Matrix<T>(n) - 2 * ((v * v.transp()) / (v.transp() * v)(0, 0));
        //std::cout << "Matrix H(" << i << "):\n" << H << "\n";
        R = H * R;
        //std::cout << "Matrix A(" << i << "):\n" << R << "\n";
        Q = Q * H;
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < i; ++j) {
            R(i, j) = T(0);
        }
    }
    //R(1, 0) = R(2, 0) = R(2, 1) = T(0);
    //std::cout << "==============DONE==============\n";
    //std::cout << "QR: " << Q * R << "\nQ: " << Q << "\nR: " << R << "\n";
    //exit(1);
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
std::vector<double> QRFindLambda (const Matrix<T> &matrix, T approx) {
    std::vector<double> lambda;
    uint64_t iteration = 0, n = matrix.size().n;
    //std::cout << "Matrix A:\n" << matrix << "\nEpsilon = " << approx << "\n";
    Matrix<T> Q, R, A(matrix);
    //std::cout << "A(k) = Q(k) * R(k)\nA(k + 1) = R(k) * Q(k)\n";
    while (1) {
        ++iteration;
        //std::cout << "------------\n";
        //std::cout << "Iteration " << iteration << ":\n";
        std::tie(Q, R) = QR(A);
        A = R * Q;
        //std::cout << "Matrix R(" << iteration - 1 << "):\n" << R << "\nMatrix Q(" << iteration - 1 << "):\n" << Q << "\n";
        //std::cout << "Matrix A(" << iteration << "):\n" << A << "\n";
        if (QRCheckForEnd(A, approx) || iteration > ITERATION_CAP) {
            break;
        }
    }
    //std::cout << "Final A:\n" << A << "\n";
    for (uint64_t i = 0; i < n; ++i) {
        T sum = 0;
        for (uint64_t j = i + 1; j < A.size().n; ++j) {
            sum += A(j, i) * A(j, i);
        }
        sum = std::sqrt(sum);
        uint64_t j = i + 1;
        if (sum <= approx) {
            //std::cout << "x[" << i << "] = " << A(i, i) << "\n";
            lambda.push_back(A(i, i));
            //lambda.push_back(1);
        } else {
            T b = -(A(i, i) + A(j, j)), c = A(i, i) * A(j, j) - A(i, j) * A(j, i), a = T(1);
            T D = b * b - 4 * a * c;
            //std::cout << "x[" << i << "] = " << -b / 2 << " + " << std::sqrt(std::abs(D)) / 2 << "i\n";
            ++i;
            //std::cout << "x[" << i << "] = " << -b / 2<< " - " << std::sqrt(std::abs(D)) / 2 << "i\n";
            lambda.push_back(-b / 2);
            lambda.push_back(-b / 2);
        }
    }
    return lambda;
}

#endif