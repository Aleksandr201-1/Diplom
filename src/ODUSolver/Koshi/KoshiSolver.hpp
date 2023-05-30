#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <vector>
#include <General/General.hpp>
#include <General/ButcherTable.hpp>
#include <ODUSolver/ToughDet.hpp>
#include <ODUSolver/IterationAlgo.hpp>
#include <ODUSolver/Koshi/KoshiTask.hpp>

std::vector<std::vector<float128_t>> KoshiSolver (SolveMethod method, const KoshiTask &task, const Matrix<float128_t> butcher, float128_t h, IterationAlgo iter_alg = IterationAlgo::SIMPLE_ITERATION, float128_t approx = 0.01);

#endif