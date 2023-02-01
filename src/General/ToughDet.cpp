#include "ToughDet.hpp"

double ToughCoeff (const Task &task) {
    double coeff = 0;
    uint64_t order = task.order;
    Matrix<double> A(order, order);
    std::vector<double> args(order, 0);
    for (uint64_t i = 0; i < order; ++i) {
        for (uint64_t j = 0; j < order; ++j) {

        }
    }
    return coeff;
}

double ToughCoeff (const std::vector<std::function<double(const std::vector<double> &)>> &system, const Task &task) {
    double coeff = 0;
    uint64_t order = task.order;
    Matrix<double> A(order, order);
    std::vector<double> args(order + 2, 0);
    args[0] = task.Xn;
    for (uint64_t i = 0; i < order; ++i) {
        args[i + 1] = 1;
        for (uint64_t j = 0; j < order; ++j) {
            A(j, i) = system[j](args);
        }
        args[i + 1] = 0;
    }
    std::cout << "\ncurrent A:\n" << A << "\n";
    std::vector<double> lambdas = QRFindLambda(A, 0.01);
    std::cout << "lambdas: ";
    printVector(lambdas);
    double max = *std::max_element(lambdas.begin(), lambdas.end());
    double min = *std::min_element(lambdas.begin(), lambdas.end());
    if (max >= 0 && min >= 0) {
        return 0;
    } else {
        max = std::abs(max);
        min = std::abs(min);
        return std::max(max / min, min / max);
    }
    // if (min > max) {
    //     std::swap(max, min);
    // }
    // //auto Matr = Matrix<double>(3, 3, {1, 2, 3, 4, 4, 0, 7, 0, 9});
    // auto Matr = Matrix<double>(2, 2, {0, 1, 1, -2});
    // std::cout << "\n" << Matr << "\n";
    // printVector(QRFindLambda(Matr, 0.01));
    // return max / min;
}