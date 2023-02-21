#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <vector>
#include "../General/General.hpp"
#include "../General/FuncMaker.hpp"
#include "../General/ButcherTable.hpp"
#include "../General/ToughDet.hpp"
//#include "../General/LSM.hpp"

const double min_h = 0.05;

enum class IterationAlgo {
    NEWTON,
    ZEIDEL
};

//std::pair<std::vector<double>, std::vector<double>> RungeKutta4 (const Task &task, double h);
std::pair<std::vector<double>, std::vector<double>> RungeKutta6 (const Task &task, const Matrix<double> butcher, double h);
std::pair<std::vector<double>, std::vector<double>> Falberg (const Task &task, const Matrix<double> butcher, double h);
std::pair<std::vector<double>, std::vector<double>> NonExplZeidel (const Task &task, const Matrix<double> butcher, double h);

std::pair<std::vector<double>, std::vector<double>> ODUSolve (SolveMethod method, const Task &task, const Matrix<double> butcher, double h, IterationAlgo iter_alg = IterationAlgo::ZEIDEL, double approx = 0.01);

#endif