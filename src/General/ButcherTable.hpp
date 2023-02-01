#ifndef BUTCHER_TABLE_HPP
#define BUTCHER_TABLE_HPP

#include "Matrix.hpp"

enum class SolveMethod {
    RUNGE_KUTTA,
    FALBERG,
    CHESKINO,
    MERSON,
    RADO,
    GAUSS,
    LOBATTO
};

Matrix<double> createButcherTable (SolveMethod method, uint64_t order, uint64_t way);

#endif