#include "Differentiation.hpp"

const std::map<DiffConfig, std::vector<float128_t>> diff_coeffs = {
    {DiffConfig::POINTS2_ORDER1_WAY1, {-1, 0, 1, 2, 1}},               // 2 points, 1 diff order, 1 way of realisation
    {DiffConfig::POINTS3_ORDER1_WAY1, {0, 0, -3, 4, -1, 2, 1}},        // 3 points, 1 diff order, 1 way of realisation
    {DiffConfig::POINTS4_ORDER1_WAY1, {0, 1, -8, 0, 8, -1, 0, 12, 1}}, // 4 points, 1 diff order, 1 way of realisation
    {DiffConfig::POINTS4_ORDER1_WAY2, {0, 0, 0, -11, 18, -9, 2, 6, 1}} // 4 points, 1 diff order, 2 way of realisation
};

// const std::vector<std::vector<float128_t>> diff_coeffs = {
//     {-1, 0, 1, 2, 1},               // 2 points, 1 diff order
//     {0, 0, -3, 4, -1, 2, 1},        // 3 points, 1 diff order
//     {0, 1, -8, 0, 8, -1, 0, 12, 1}, // 4 points, 1 diff order
//     {0, 0, 0, -11, 18, -9, 2, 6, 1} // 4 points, 1 diff order
// };

// float128_t derivative1point4 (const std::function<float128_t(std::vector<float128_t> &)> &f, const std::vector<float128_t> &X, float128_t h, uint64_t idx) {
//     std::vector<float128_t> args(X);
//     float128_t a, b, c, d;
//     args[idx] += 2 * h;
//     a = f(args);
//     args[idx] -= h;
//     b = f(args);
//     args[idx] -= 2 * h;
//     c = f(args);
//     args[idx] -= h;
//     d = f(args);
//     return (-a + 8 * b - 8 * c + d) / 12 * h;
// }

float128_t derivative (const std::function<float128_t(std::vector<float128_t> &)> &f,
                       const std::vector<float128_t> &X,
                       float128_t h,
                       uint64_t idx,
                       DiffConfig diff) {
    std::vector<float128_t> args(X);
    float128_t ans = 0;
    const auto &coeffs = diff_coeffs.find(diff)->second;
    uint64_t size = coeffs.size();
    uint64_t points = (size - 1) / 2;
    args[idx] -= points * h;
    for (uint64_t i = 0; i < size - 2; ++i) {
        args[idx] += h;
        ans += coeffs[i] * f(args);
    }
    ans /= std::pow(h, coeffs[size - 1]) * coeffs[size - 2];
    return ans;
}

float128_t derivative (const std::function<float128_t(float128_t)> &f,
                       float128_t x,
                       float128_t h,
                       DiffConfig diff) {
    float128_t ans = 0;
    const auto &coeffs = diff_coeffs.find(diff)->second;
    uint64_t size = coeffs.size();
    uint64_t points = (size - 1) / 2;
    x -= points * h;
    for (uint64_t i = 0; i < size - 2; ++i) {
        x += h;
        ans += coeffs[i] * f(x);
    }
    ans /= std::pow(h, coeffs[size - 1]) * coeffs[size - 2];
    return ans;
}