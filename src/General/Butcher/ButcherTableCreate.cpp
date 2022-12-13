#include <iostream>
#include "../Matrix.hpp"

int main () {
    //BT-RK-4-1.bin
    Matrix<double> matrix(5, 5, {
        0, 0, 0, 0, 0,
        1.0 / 2.0, 1.0 / 2.0, 0, 0, 0,
        1.0 / 2.0, 0, 1.0 / 2.0, 0, 0,
        1, 0, 0, 1, 0,
        0, 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0
    });
    
    std::string filename = "BT-RK-4-1.bin";

    std::ofstream file(filename);
    file << matrix;
    file.close();

    return 0;
}