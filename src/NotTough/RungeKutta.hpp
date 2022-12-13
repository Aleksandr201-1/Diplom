#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <vector>
#include "../General/General.hpp"
#include "../General/FuncMaker.hpp"

std::pair<std::vector<double>, std::vector<double>> RungeKutta4 (const Task &task, double h);
std::pair<std::vector<double>, std::vector<double>> RungeKutta6 (const Task &task, const Matrix<double> butcher, double h);

#endif