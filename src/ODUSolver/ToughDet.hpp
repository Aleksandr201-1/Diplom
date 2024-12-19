#ifndef TOUGH_DET
#define TOUGH_DET

#include <algorithm>
#include <NumericMethods/Differentiation.hpp>
#include <NumericMethods/LU.hpp>
#include <NumericMethods/QR.hpp>

double getToughK (const Matrix<double> &K);

double ToughCoeff (const std::vector<std::function<double (const std::vector<double> &)>> &system, const std::vector<double> &Y, const std::vector<double> &Yprev);

#endif