#include "ButcherTable.hpp"
//#include <iostream>

const std::map<std::string, SolveMethod> methods = {
    {"RungeKutta", SolveMethod::RUNGE_KUTTA},
    {"Falberg", SolveMethod::FALBERG},
    {"Cheskino", SolveMethod::CHESKINO},
    {"Merson", SolveMethod::MERSON},
    {"Rado", SolveMethod::RADO},
    {"Gauss", SolveMethod::GAUSS},
    {"Lobatto", SolveMethod::LOBATTO},
    {"LStableDiagonal", SolveMethod::L_STABLE_DIAGONAL},
    {"DormanPrince", SolveMethod::DORMAN_PRINCE}//,
    //{"NotAMethod", SolveMethod::NOT_A_METHOD}
};

std::string solveMethodToString (SolveMethod method) {
    auto check = [method] (const auto &pair) -> bool {
        return pair.second == method;
    };
    auto result = std::find_if(methods.begin(), methods.end(), check);
    if (result != methods.end()) {
        return result->first;
    }
    return "NotAMethod";
}

SolveMethod stringToSolveMethod (const std::string &str) {
    auto it = methods.find(str);
    if (it != methods.end()) {
        return it->second;
    }
    return SolveMethod::NOT_A_METHOD;
    // if (str == "RungeKutta") {
    //     return SolveMethod::RUNGE_KUTTA;
    // }
    // if (str == "Falberg") {
    //     return SolveMethod::FALBERG;
    // }
    // if (str == "Cheskino") {
    //     return SolveMethod::CHESKINO;
    // }
    // if (str == "Merson") {
    //     return SolveMethod::MERSON;
    // }
    // if (str == "Rado") {
    //     return SolveMethod::RADO;
    // }
    // if (str == "Gauss") {
    //     return SolveMethod::GAUSS;
    // }
    // if (str == "Lobatto") {
    //     return SolveMethod::LOBATTO;
    // }
    // if (str == "LStableDiagonal") {
    //     return SolveMethod::L_STABLE_DIAGONAL;
    // }
    // if (str == "DormanPrince") {
    //     return SolveMethod::DORMAN_PRINCE;
    // }
    // return SolveMethod::NOT_A_METHOD;
    //throw std::logic_error("stringToSolveMethod: method \"" + str + "\" doesnt exist");
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