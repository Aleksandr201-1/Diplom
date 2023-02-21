#include "ButcherTable.hpp"
//#include <iostream>

std::string solveMethodToString (SolveMethod method) {
    switch (method) {
        case SolveMethod::RUNGE_KUTTA:
            return "RungeKutta";
        case SolveMethod::FALBERG:
            return "Falberg";
        case SolveMethod::CHESKINO:
            return "Cheskino";
        case SolveMethod::MERSON:
            return "Merson";
        case SolveMethod::RADO:
            return "Rado";
        case SolveMethod::GAUSS:
            return "Gauss";
        case SolveMethod::LOBATTO:
            return "Lobatto";
        default:
            return "";
    }
}

Matrix<double> createButcherTable (SolveMethod method, uint64_t order, uint64_t way) {
    Matrix<double> butcher;
    std::string filename = "./src/General/Butcher/BT-" + solveMethodToString(method) + "-" + std::to_string(order) + "-" + std::to_string(way) + ".bin";
    //std::cout << filename << "\n";
    //system("ls -l");
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::logic_error("createButcherTable: cant open file with name \"" + filename + "\"");
    }
    file >> butcher;
    file.close();
    return butcher;
}