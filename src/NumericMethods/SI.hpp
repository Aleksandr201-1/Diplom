#ifndef SIMPLE_ITERATION_HPP
#define SIMPLE_ITERATION_HPP

//Баев, 
//Турбулентные течения реагирубющих газов, Либби, Вильямс
//Сполдинг, Горение и массообмен

#include <Math/Matrix.hpp>

template <class T>
T SINormal (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;
    T max = 0;
    for (uint64_t i = 0; i < n; ++i) {
        T tmp = 0;
        for (uint64_t j = 0; j < n; ++j) {
            tmp += std::abs(matrix(i, j));
        }
        if (tmp > max) {
            max = tmp;
        }
    }
    return max;
}

template <class T>
T FindEpsilon (const std::vector<T> &oldX, const std::vector<T> &newX, T a) {
    T eps = 0;
    for (uint64_t i = 0; i < oldX.size(); ++i) {
        T tmp = std::abs(oldX[i] - newX[i]);
        if (tmp > eps) {
            eps = tmp;
        }
    }
    return a < T(1) ? eps * a / (T(1) - a) : eps;
}

template <class T>
std::vector<T> SISolveSLAE (const Matrix<T> &matrix, const std::vector<T> &b, T approx) {
    uint64_t n = matrix.size().n;
    std::vector<T> beta(n), x(n, 0);
    Matrix<T> alpha(n, n);

    if (!matrix.isSquare() || n != b.size()) {
        return x;
    }
    if (n != b.size()) {
        return x;
    }

    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < i; ++j) {
            alpha(i, j) = -matrix(i, j) / matrix(i, i);
        }
        for (uint64_t j = i + 1; j < n; ++j) {
            alpha(i, j) = -matrix(i, j) / matrix(i, i);
        }
        beta[i] = b[i] / matrix(i, i);
        x[i] = beta[i];
    }
    T a = SINormal(alpha);
    uint64_t iteration = 0;
    T epsilon = T(0);
    if (a > T(1)) {
        //std::cout << "Norma greater than 1. Using second end criteria (||x(k) - x(k - 1)|| < e)\n";
    } else {
        //std::cout << "Norma below 1. Using first end criteria (||x(k) - x*|| < e)\n";
    }

    //std::cout << "\n------------\n";
    //std::cout << "Iteration 0:\n";
    //printVector("x", x);
    while (1) {
        ++iteration;
        //std::cout << "------------\n";
        //std::cout << "Iteration " << iteration << ":\n";
        std::vector<T> newX(x);
        for (uint64_t i = 0; i < n; ++i) {
            newX[i] = beta[i];
            //std::cout << "x(" << iteration << ")[" << i << "] = beta[" << i << "]";
            //if (method == SI_ZEIDEL_METHOD) {
                for (uint64_t j = 0; j < i; ++j) {
                    newX[i] += alpha(i, j) * newX[j];
                    //std::cout << " + alpha(" << i << ", " << j << ") * x(" << iteration << ")[" << j << "]";
                }
                for (uint64_t j = i; j < n; ++j) {
                    newX[i] += alpha(i, j) * x[j];
                    //std::cout << " + alpha(" << i << ", " << j << ") * x(" << iteration - 1 << ")[" << j << "]";
                }
            //} else {
            //     for (uint64_t j = 0; j < n; ++j) {
            //         newX[i] += alpha(i, j) * x[j];
            //         //std::cout << " + alpha(" << i << ", " << j << ") * x(" << iteration - 1 << ")[" << j << "]";
            //     }
            // }
            //std::cout << "\n";
        }
        epsilon = FindEpsilon(x, newX, a);
        x = newX;
        if (epsilon < approx || iteration > 200) {
            break;
        }
    }
    return x;
}

#endif