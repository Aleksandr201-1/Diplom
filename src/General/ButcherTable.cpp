#include "ButcherTable.hpp"
//#include <iostream>
Matrix<double> createButcherTable (uint64_t order, uint64_t way) {
    Matrix<double> butcher;
    std::string filename = "./src/General/Butcher/BT-RK-" + std::to_string(order) + "-" + std::to_string(way) + ".bin";
    //std::cout << filename << "\n";
    //system("ls -l");
    std::ifstream file(filename);
    file >> butcher;
    file.close();
    return butcher;
}