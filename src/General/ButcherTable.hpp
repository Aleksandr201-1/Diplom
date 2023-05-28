#ifndef BUTCHER_TABLE_HPP
#define BUTCHER_TABLE_HPP

#include <algorithm>
#include <map>
#include "General.hpp"
#include "Matrix.hpp"
#include "Enum.hpp"

enum class SolveMethod {
    //using enum UsableList;
    RUNGE_KUTTA,
    FALBERG,
    CHESKINO,
    MERSON,
    RADO,
    GAUSS,
    LOBATTO,
    L_STABLE_DIAGONAL,
    DORMAN_PRINCE,
    ERROR
};

std::string solveMethodToString (SolveMethod method);

SolveMethod stringToSolveMethod (const std::string &str);

Matrix<float128_t> createButcherTable (SolveMethod method, uint64_t order, uint64_t way);

#endif