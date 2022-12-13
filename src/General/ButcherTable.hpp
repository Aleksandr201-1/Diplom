#ifndef BUTCHER_TABLE_HPP
#define BUTCHER_TABLE_HPP

#include "Matrix.hpp"

Matrix<double> createButcherTable (uint64_t order, uint64_t way);

#endif