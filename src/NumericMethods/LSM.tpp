#ifndef LSM_TPP

template <typename T>
std::vector<T> LeastSquareMethod (const std::vector<T> &X, const std::vector<T> &Y, uint64_t n) {
    Matrix<T> M(n + 1, n + 1);
    std::vector<T> y(n + 1);
    for (uint64_t i = 0; i < n + 1; ++i) {
        T tmp = 0;
        for (uint64_t j = 0; j < X.size(); ++j) {
            tmp += std::pow(X[j], i);
        }
        M(0, i) = tmp;
    }
    for (uint64_t i = 1; i < n + 1; ++i) {
        T tmp = 0;
        for (uint64_t j = 0; j < n; ++j) {
            M(i, j) = M(i - 1, j + 1);
        }
        for (uint64_t j = 0; j < X.size(); ++j) {
            tmp += std::pow(X[j], n + i);
        }
        M(i, n) = tmp;
    }
    for (uint64_t i = 0; i < y.size(); ++i) {
        T tmp = 0;
        for (uint64_t j = 0; j < Y.size(); ++j) {
            tmp += Y[j] * std::pow(X[j], i);
        }
        y[i] = tmp;
    }
    return LUsolveSLAE(M, y);
}

template <typename T>
T ErrorSquareSum (const std::vector<T> &X, const std::vector<T> &Y, const std::vector<T> &poly) {
    T ans = 0;
    for (uint64_t i = 0; i < X.size(); ++i) {
        T tmp = 0;
        for (uint64_t j = 0; j < poly.size(); ++j) {
            tmp += poly[j] * std::pow(X[i], j);
        }
        ans += (tmp - Y[i]) * (tmp - Y[i]);
    }
    return ans;
}

template <typename T>
T LSMFunc (const std::vector<T> &coeff, T x) {
    T ans = 0;
    for (uint64_t i = 0; i < coeff.size(); ++i) {
        ans += coeff[i] * std::pow(x, i);
    }
    return ans;
}

template <typename T>
std::string LSMToText (const std::vector<T> &poly, uint64_t precision) {
    std::string func = std::to_string(poly[0]);
    for (uint64_t i = 1; i < poly.size(); ++i) {
        if (poly[i] < 0) {
            func += " - ";
        } else {
            func += " + ";
        }
        func += toString(std::abs(poly[i]), precision) + "x^" + std::to_string(i);
    }
    return func;
}

#endif