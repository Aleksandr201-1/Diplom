#ifndef REPORT_GENERATOR_HPP
#define REPORT_GENERATOR_HPP

#include <iostream>
#include <functional>
#include <iomanip>
#include <chrono>
#include <ctime>
//#include <regex>
#include "../General/General.hpp"
#include "../General/LSM.hpp"
#include "../General/ButcherTable.hpp"

enum class ReportType {
    TXT,
    TEX,
    PDF
};

struct ReportInfo {
    SolveMethod method;
    uint64_t order, way;
    Matrix<double> butcher;

    std::vector<std::string> input_task;
    Task task;
    double h, tough_coeff;

    //std::vector<double> X, Ynum;
    std::pair<std::vector<double>, std::vector<double>> solution;
    //std::function<double (double)> analitic;
    FunctionalTree analitic;

    ReportInfo ();
    ~ReportInfo ();
    //std::function<double (const std::vector<double> &)> analitic;
};

void generateReport (const ReportInfo &info, ReportType type, std::ostream &out = std::cout);

#endif