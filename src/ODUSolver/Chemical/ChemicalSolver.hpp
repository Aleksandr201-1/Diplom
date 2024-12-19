#ifndef CHEMICAL_SOLVER_HPP
#define CHEMICAL_SOLVER_HPP

#include <vector>
//#include <General/FuncMaker.hpp>
#include <Butcher/ButcherTable.hpp>
#include <ODUSolver/ToughDet.hpp>
#include <ODUSolver/IterationAlgo.hpp>
#include <ODUSolver/Chemical/ChemicalTask.hpp>

std::vector<std::vector<double>> drop (const std::vector<std::vector<double>> &vv, uint64_t newSize);

double NewtonFindT (const std::function<double (double)> &f, double x, double approx, DiffConfig conf = DiffConfig::POINTS2_ORDER1_WAY3);

std::vector<std::vector<double>> ChemicalSolver (SolveMethod method,
                                                     const ChemicalSystem &task,
                                                     const Matrix<double> butcher,
                                                     double h_min,
                                                     double h_max,
                                                     double h_last,
                                                     IterationAlgo iter_alg,
                                                     double approx,
                                                     ReactionType type);

#endif