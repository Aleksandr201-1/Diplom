#ifndef INTEGRAL_HPP
#define INTEGRAL_HPP

#include <cstdint>
#include <cmath>
#include <functional>
#include <vector>
#include <NumericMethods/Interpolation.hpp>

double IntegralRectangle (const std::function<double(double)> &f, double X1, double X2, double step);
double IntegralRectangle (const std::vector<double> &X, const std::vector<double> &Y);
double IntegralRectangle (const std::function<double(uint64_t)> &fx, const std::function<double(uint64_t)> &fy, uint64_t size);

double IntegralTrapeze   (const std::function<double(double)> &f, double X1, double X2, double step);
double IntegralTrapeze (const std::vector<double> &X, const std::vector<double> &Y);
double IntegralTrapeze (const std::function<double(uint64_t)> &fx, const std::function<double(uint64_t)> &fy, uint64_t size);

double IntegralSimpson   (const std::function<double(double)> &f, double X1, double X2, double step);
double IntegralSimpson (const std::vector<double> &X, const std::vector<double> &Y);
double IntegralSimpson (const std::function<double(uint64_t)> &fx, const std::function<double(uint64_t)> &fy, uint64_t size);

double IntegralRunge     (double ans1, double ans2, double k, double p);

#endif