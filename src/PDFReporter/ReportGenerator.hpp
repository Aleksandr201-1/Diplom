#ifndef REPORT_GENERATOR_HPP
#define REPORT_GENERATOR_HPP

#include <iostream>
#include <functional>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <General/General.hpp>
#include <ODUSolver/Koshi/KoshiTask.hpp>
#include <ODUSolver/IterationAlgo.hpp>
#include <General/LSM.hpp>
#include <General/ButcherTable.hpp>

enum class ReportType {
    TXT,
    TEX,
    PDF
};

struct ReportInfo {
    TaskType type;
    SolveMethod method;
    IterationAlgo algo;
    uint64_t order, way;
    Matrix<float128_t> butcher;

    std::vector<std::string> input_task;
    Task *task;
    float128_t h, tough_coeff, approx;

    bool multigraph;
    std::vector<std::pair<std::string, std::string>> graph_info;
    std::vector<std::vector<float128_t>> solution;
    std::vector<FunctionalTree> analitic;

    std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> table;
    uint64_t workTime;

    std::string fileInput;

    ReportInfo ();
    ~ReportInfo ();
};

void generateReport (const ReportInfo &info, ReportType type, std::ostream &out = std::cout);

#endif