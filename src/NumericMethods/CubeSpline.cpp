#include <NumericMethods/CubeSpline.hpp>

Matrix<double> CubeSpline (const std::vector<double> &X, const std::vector<double> &Y) {
    uint64_t n = X.size() - 1;
    Matrix<double> M(n, 4);
    for (uint64_t i = 0; i < Y.size() - 1; ++i) {
        M(i, 0) = Y[i];
    }
    --n;
    Matrix<double> S(n, n);
    std::vector<double> c(n);
    S(0, 0) = 2 * (X[2] - X[0]);
    S(0, 1) = X[2] - X[1];
    c[0] = 3 * ((Y[2] - Y[1]) / (X[2] - X[1]) - (Y[1] - Y[0]) / (X[1] - X[0]));
    for (uint64_t i = 1; i < c.size() - 1; ++i) {
        S(i, i - 1) = X[i + 1] - X[i];
        S(i, i) = 2 * (X[i + 2] - X[i]);
        S(i, i + 1) = X[i + 2] - X[i + 1];
        c[i] = 3 * ((Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]));
    }
    S(n - 1, n - 2) = X[n - 1] - X[n - 2];
    S(n - 1, n - 1) = 2 * (X[n] - X[n - 2]);
    c[n - 1] = 3 * ((Y[n + 1] - Y[n]) / (X[n + 1] - X[n]) - (Y[n] - Y[n - 1]) / (X[n] - X[n - 1]));
    c = LUsolveSLAE(S, c);
    for (uint64_t i = 0; i < c.size(); ++i) {
        M(i + 1, 2) = c[i];
    }
    for (uint64_t i = 0; i < n; ++i) {
        M(i, 1) = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) - 1.0 / 3.0 * (X[i + 1] - X[i]) * (M(i + 1, 2) + 2 * M(i, 2));
        M(i, 3) = (M(i + 1, 2) - M(i, 2)) / (3 * (X[i + 1] - X[i]));
    }
    M(n, 1) = (Y[n + 1] - Y[n]) / (X[n + 1] - X[n]) - 2.0 / 3.0 * (X[n + 1] - X[n]) * M(n, 2);
    M(n, 3) = -M(n, 2) / (3 * (X[n + 1] - X[n]));
    return M;
}

double CubeSplineFunc (const std::vector<double> &X, const Matrix<double> &M, double x) {
    double ans = 0;
    uint64_t i = 0;
    while (x > X[i + 1]) {
        ++i;
    }
    for (uint64_t j = 0; j < M.size().m; ++j) {
        ans += M(i, j) * std::pow(x - X[i], j);
    }
    return ans;
}

std::string CubeSplineToText (const std::vector<double> &X, const Matrix<double> &M, double x) {
    std::string str = "0";
    uint64_t i = 0;
    while (x > X[i + 1]) {
        ++i;
    }
    for (uint64_t j = 0; j < M.size().m; ++j) {
        str += "+" + toString(M(i, j), 2) + "* (x - " + toString(X[i], 2) + ") ** " + toString(j, 2);
        //str += "+" + std::to_string(M(i, j)) + "* (x - " + std::to_string(X[i]) + ") ** " + std::to_string(j);
    }
    return str;
}