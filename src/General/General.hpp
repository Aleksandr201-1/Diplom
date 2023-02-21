#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <cmath>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <functional>
#include <iostream>
#include "FuncMaker.hpp"
// #include "ButcherTable.hpp"
// //#include "ToughDet.hpp"
// #include "LSM.hpp"

struct Task {
    std::vector<FunctionalTree> trees;
    std::vector<std::function<double(const std::vector<double> &)>> odu_system;
    std::vector<double> Y;
    double X0, Xn;
    uint64_t order;
};

const uint64_t ITERATION_CAP = 200;
const uint64_t PRECISION = 2;

uint64_t getOrder (const std::string &task);

void printVector(const std::vector<double> &vec);

double findEpsillon ();

bool isEqual(double x, double y);

double derivative (const std::function<double(double)> &f, double x, uint64_t degree = 1);

std::string toString (double val, uint64_t precision);

std::string readLine ();

Task getTaskInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn);

#endif