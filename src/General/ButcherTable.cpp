#include "ButcherTable.hpp"
//#include <iostream>

const std::map<SolveMethod, std::string> solve_methods = {
    {SolveMethod::RUNGE_KUTTA, "RungeKutta"},
    {SolveMethod::FALBERG, "Falberg"},
    {SolveMethod::CHESKINO, "Cheskino"},
    {SolveMethod::MERSON, "Merson"},
    {SolveMethod::RADO, "Rado"},
    {SolveMethod::GAUSS, "Gauss"},
    {SolveMethod::LOBATTO, "Lobatto"},
    {SolveMethod::L_STABLE_DIAGONAL, "LStableDiagonal"},
    {SolveMethod::DORMAN_PRINCE, "DormanPrince"}
};

std::string solveMethodToString (SolveMethod method) {
    return enumToString(method, solve_methods);
}

SolveMethod stringToSolveMethod (const std::string &str) {
    return stringToEnum(str, solve_methods);
    //return SolveMethod::ERROR;
}

Matrix<float128_t> createButcherTable (SolveMethod method, uint64_t order, uint64_t way) {
    Matrix<float128_t> butcher;
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