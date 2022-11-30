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
#include "LSM.hpp"

struct Task {
    std::vector<FunctionalTree> trees;
    std::vector<double> Y;
    double X0, Xn, a, b, h;
    uint64_t n;
};

const uint64_t ITERATION_CAP = 200;
const uint64_t PRECISION = 2;

uint64_t getOrder (const std::string &task);

void printVector(const std::vector<double> &vec);

double findEpsillon ();

bool isEqual(double x, double y);

double derivative (const std::function<double(double)> &f, double x, uint64_t degree = 1);

std::string toString (double val, uint64_t precision);

Task getTaskInfo(const std::vector<std::string> &system);

#endif