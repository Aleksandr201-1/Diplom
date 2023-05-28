#ifndef INPUT_HANDLE_HPP
#define INPUT_HANDLE_HPP

#include <iostream>
#include <vector>
#include <PDFReporter/ReportGenerator.hpp>

void help (const std::string &name);

std::string readLine (std::istream &input = std::cin);

void argsHandler (const std::vector<std::string> &args, ReportInfo &info);

#endif