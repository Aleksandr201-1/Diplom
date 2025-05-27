#include <NumericMethods/Integral.hpp>

double IntegralRectangle (const std::function<double(double)> &f, double X1, double X2, double step) {
    double ans = 0;
    uint64_t size = (X2 - X1) / step;
    for (uint64_t i = 0; i < size; ++i) {
        ans += f((X1 * 2 + step) / 2) * step;
        X1 += step;
    }
    return ans;
}

double IntegralRectangle (const std::vector<double> &X, const std::vector<double> &Y) {
    double ans = 0;
    double X1 = X[0];
    auto f = [&] (double x) -> double {
        return LinearInterpolation(X, Y, x);
    };
    for (uint64_t i = 0; i < X.size() - 1; ++i) {
        double step = X[i] - X[i + 1];
        ans += f((X1 * 2 + step) / 2) * step;
        X1 += step;
    }
    return ans;
}

double IntegralRectangle (const std::function<double(uint64_t)> &fx, const std::function<double(uint64_t)> &fy, uint64_t size) {
    double ans = 0;
    double X1 = fx(0);
    for (uint64_t i = 0; i < size - 1; ++i) {
        double step = fx(i) - fx(i + 1);
        ans += fy(i) * step;
        X1 += step;
    }
    return ans;
}

double IntegralTrapeze (const std::function<double(double)> &f, double X1, double X2, double step) {
    double ans = 0;
    uint64_t size = (X2 - X1) / step;
    ans = f(X1) + f(X2);
    ans /= 2;
    for (uint64_t i = 1; i < size; ++i) {
        ans += f(X1 + i * step);
    }
    ans *= step;
    return ans;
}

double IntegralTrapeze (const std::vector<double> &X, const std::vector<double> &Y) {
    double ans = 0;
    uint64_t size = Y.size();
    ans = Y[0] + Y[size - 1];
    ans /= 2;
    for (uint64_t i = 1; i < size; ++i) {
        double step = X[i] - X[i - 1];
        ans += Y[i] * step;
    }
    return ans;
}

double IntegralTrapeze (const std::function<double(uint64_t)> &fx, const std::function<double(uint64_t)> &fy, uint64_t size) {
    double ans = 0;
    ans = fy(0) + fy(size - 1);
    ans /= 2;
    for (uint64_t i = 1; i < size; ++i) {
        double step = fx(i) - fx(i - 1);
        ans += fy(i) * step;
    }
    return ans;
}

double IntegralSimpson (const std::function<double(double)> &f, double X1, double X2, double step) {
    double ans = 0;
    uint64_t size = (X2 - X1) / step;
    ans += f(X1) + f(X2);
    for (uint64_t i = 1; i < size; ++i) {
        ans += f(X1 + step * i) * (i & 1 ? 4 : 2);
    }
    ans *= step / 3;
    return ans;
}

double IntegralSimpson (const std::vector<double> &X, const std::vector<double> &Y) {
    double ans = 0;
    ans += Y.front() + Y.back();
    for (uint64_t i = 1; i < X.size(); ++i) {
        ans += Y[i] * (i & 1 ? 4 : 2) * (X[i - 1] - X[i]) / 3;
    }
    return ans;
}

double IntegralSimpson (const std::function<double(uint64_t)> &fx, const std::function<double(uint64_t)> &fy, uint64_t size) {
    double ans = 0;
    ans += fy(0) + fy(size - 1);
    for (uint64_t i = 1; i < size; ++i) {
        ans += fy(i) * (i & 1 ? 4 : 2) * (fx(i - 1) - fx(i)) / 3;
    }
    return ans;
}

double IntegralRunge (double ans1, double ans2, double k, double p) {
    return std::abs((ans1 - ans2) / (std::pow(k, p) - 1));
}