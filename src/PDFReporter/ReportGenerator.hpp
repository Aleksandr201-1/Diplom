#ifndef REPORT_GENERATOR_HPP
#define REPORT_GENERATOR_HPP

#include <iostream>
#include <functional>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <ODUSolver/Koshi/KoshiTask.hpp>
#include <ODUSolver/Chemical/ChemicalTask.hpp>
#include <ODUSolver/IterationAlgo.hpp>
#include <NumericMethods/LSM.hpp>
#include <Butcher/ButcherTable.hpp>

std::vector<std::vector<double>> getAnaliticSolution (const std::vector<double> &X, const std::vector<FuncMaker> &func);

std::tuple<double, double> getAnaliticCompare (const std::vector<double> &Yn, const std::vector<double> &Ya);

enum class ReportType {
    TXT,
    TEX,
    PDF
};

struct ReportInfo {
    TaskType type;
    SolveMethod method;
    ReactionType react;
    IterationAlgo algo;
    uint64_t order, way;
    Matrix<double> butcher;

    std::vector<std::string> input_task;
    Task *task;
    double h_min, h_max, h_last, tough_coeff, approx;

    bool multigraph;
    std::vector<std::pair<std::string, std::string>> graph_info;
    std::vector<std::vector<double>> solution;
    std::vector<FuncMaker> analitic;

    std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> table;
    uint64_t workTime;

    std::string fileInput;

    ReportInfo ();
    ~ReportInfo ();
};

void generateReport (const ReportInfo &info, ReportType type, std::ostream &out = std::cout);

#endif