#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <SFML/System/Vector2.hpp>
#include <cstdint>
#include <vector>
#include <functional>
//#include <iostream>

double LinearInterpolation (const std::vector<std::pair<double, double>> &points, double x);
double LinearInterpolation (const std::vector<sf::Vector2f> &points, double x);

template <typename T>
T LinearInterpolation (const std::vector<T> &X, const std::vector<T> &Y, T x) {
    T ans = 0;
    uint64_t i = 0;
    while (x > X[i + 1]) {
        ++i;
        if (i + 1 == X.size()) {
            --i;
            break;
        }
    }
    T coeff = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
    return (Y[i] - X[i] * coeff) + coeff * x;
}

#endif