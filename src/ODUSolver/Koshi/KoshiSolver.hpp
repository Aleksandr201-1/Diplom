#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <vector>
#include <Butcher/ButcherTable.hpp>
#include <ODUSolver/ToughDet.hpp>
#include <ODUSolver/IterationAlgo.hpp>
#include <ODUSolver/Koshi/KoshiTask.hpp>

std::vector<std::vector<double>> KoshiSolver (SolveMethod method, const KoshiTask &task, const Matrix<double> butcher, double h, IterationAlgo iter_alg = IterationAlgo::SIMPLE_ITERATION, double approx = 0.01);

#endif