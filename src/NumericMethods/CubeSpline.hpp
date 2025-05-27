#ifndef CUBE_SPLINE_HPP
#define CUBE_SPLINE_HPP

#include <vector>
#include <cstdint>
#include <cmath>
#include <Math/Matrix.hpp>
#include <NumericMethods/LU.hpp>
#include <General/FloatToString.hpp>

Matrix<double> CubeSpline (const std::vector<double> &X, const std::vector<double> &Y);
double CubeSplineFunc (const std::vector<double> &X, const Matrix<double> &M, double x);
std::string CubeSplineToText (const std::vector<double> &X, const Matrix<double> &M, double x);

#endif