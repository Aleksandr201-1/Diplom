#ifndef CHEMICAL_SOLVER_HPP
#define CHEMICAL_SOLVER_HPP

#include <vector>
#include <General/General.hpp>
//#include <General/FuncMaker.hpp>
#include <General/ButcherTable.hpp>
#include <ODUSolver/ToughDet.hpp>
#include <ODUSolver/IterationAlgo.hpp>
#include <ODUSolver/Chemical/ChemicalTask.hpp>

std::vector<std::vector<float128_t>> ChemicalSolver (SolveMethod method,
                                                     const ChemicalSystem &task,
                                                     const Matrix<float128_t> butcher,
                                                     float128_t h,
                                                     IterationAlgo iter_alg,
                                                     float128_t approx,
                                                     bool Tconst);

#endif