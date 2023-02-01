#ifndef TOUGH_DET
#define TOUGH_DET

#include "General.hpp"
#include "QR.hpp"

double ToughCoeff (const Task &task);

double ToughCoeff (const std::vector<std::function<double(const std::vector<double> &)>> &system, const Task &task);

#endif