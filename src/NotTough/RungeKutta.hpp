#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <vector>
#include "../General/General.hpp"
#include "../General/FuncMaker.hpp"
#include "../General/ButcherTable.hpp"
#include "../General/ToughDet.hpp"
//#include "../General/LSM.hpp"

std::pair<std::vector<double>, std::vector<double>> RungeKutta4 (const Task &task, double h);
std::pair<std::vector<double>, std::vector<double>> RungeKutta6 (const Task &task, const Matrix<double> butcher, double h);

std::pair<std::vector<double>, std::vector<double>> Falberg (const Task &task, const Matrix<double> butcher, double h);
std::pair<std::vector<double>, std::vector<double>> NonExpl (const Task &task, const Matrix<double> butcher, double h);

#endif