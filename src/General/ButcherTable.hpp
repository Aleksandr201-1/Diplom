#ifndef BUTCHER_TABLE_HPP
#define BUTCHER_TABLE_HPP

#include <algorithm>
#include <map>
#include "Matrix.hpp"

enum class SolveMethod {
    RUNGE_KUTTA,
    FALBERG,
    CHESKINO,
    MERSON,
    RADO,
    GAUSS,
    LOBATTO,
    L_STABLE_DIAGONAL,
    DORMAN_PRINCE,
    NOT_A_METHOD
};

std::string solveMethodToString (SolveMethod method);

SolveMethod stringToSolveMethod (const std::string &str);

Matrix<double> createButcherTable (SolveMethod method, uint64_t order, uint64_t way);

#endif