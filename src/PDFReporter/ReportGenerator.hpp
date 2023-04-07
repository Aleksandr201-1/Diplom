#ifndef REPORT_GENERATOR_HPP
#define REPORT_GENERATOR_HPP

#include <iostream>
#include <functional>
#include <iomanip>
#include <chrono>
#include <ctime>
//#include <regex>
#include <General/General.hpp>
#include <General/IterationAlgo.hpp>
#include <General/LSM.hpp>
#include <General/ButcherTable.hpp>
// #include "../General/General.hpp"
// #include "../General/IterationAlgo.hpp"
// #include "../General/LSM.hpp"
// #include "../General/ButcherTable.hpp"

enum class ReportType {
    TXT,
    TEX,
    PDF
};

struct ReportInfo {
    SolveMethod method;
    IterationAlgo algo;
    uint64_t order, way;
    Matrix<double> butcher;

    std::vector<std::string> input_task;
    Task task;
    double h, tough_coeff, approx;

    //std::vector<double> X, Ynum;
    bool multigraph;
    std::vector<std::vector<double>> solution;
    //std::pair<std::vector<double>, std::vector<double>> solution;
    //std::function<double (double)> analitic;
    std::vector<FunctionalTree> analitic;

    ReportInfo ();
    ~ReportInfo ();
    //std::function<double (const std::vector<double> &)> analitic;
};

void generateReport (const ReportInfo &info, ReportType type, std::ostream &out = std::cout);

#endif