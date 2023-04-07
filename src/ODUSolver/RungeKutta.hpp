#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <vector>
#include <General/General.hpp>
#include <General/FuncMaker.hpp>
#include <General/ButcherTable.hpp>
#include <General/ToughDet.hpp>
#include <General/IterationAlgo.hpp>
// #include "../General/General.hpp"
// #include "../General/FuncMaker.hpp"
// #include "../General/ButcherTable.hpp"
// #include "../General/ToughDet.hpp"
// #include "../General/IterationAlgo.hpp"
//#include "../General/LSM.hpp"

//const double min_h = 0.05;
const uint64_t MIN_H_SIZE = 100;
const uint64_t MAX_H_SIZE = 10;

// enum class IterationAlgo {
//     NEWTON,
//     ZEIDEL,
//     SIMPLE_ITERATION
// };

//std::pair<std::vector<double>, std::vector<double>> RungeKutta4 (const Task &task, double h);
std::vector<std::vector<double>> RungeKutta (const Task &task, const Matrix<double> butcher, double h);
std::vector<std::vector<double>> Falberg (const Task &task, const Matrix<double> butcher, double h);
std::vector<std::vector<double>> NonExpl (const Task &task, const Matrix<double> butcher, double h, IterationAlgo iter_alg, double approx);

std::vector<std::vector<double>> ODESolve (SolveMethod method, const Task &task, const Matrix<double> butcher, double h, IterationAlgo iter_alg = IterationAlgo::SIMPLE_ITERATION, double approx = 0.01);

#endif