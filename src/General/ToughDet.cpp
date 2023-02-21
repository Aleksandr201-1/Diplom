#include "ToughDet.hpp"

double ToughCoeff (const Task &task) {
    double coeff = 0, curr_coeff;
    uint64_t order = task.order;
    auto &system = task.odu_system;
    Matrix<double> A(order, order);
    std::vector<double> args(order + 2, 0);
    double step = (task.Xn - task.X0) / 10;
    for (double x = task.X0; x <= task.Xn; x += step) {
        args[0] = x;
        for (uint64_t i = 0; i < order; ++i) {
            args[i + 1] = 1;
            for (uint64_t j = 0; j < order; ++j) {
                A(j, i) = system[j](args);
            }
            args[i + 1] = 0;
        }
        //std::cout << std::scientific << A(3, 0) << "\n";
        // if (A(3, 0) == 0.0) {
        //     std::cout << "good\n";
        // } else {
        //     std::cout << "bad\n";
        // }
        std::cout << "\ncurrent A:\n" << A << "\n";
        std::vector<double> lambdas = QRFindLambda(A, 0.01);
        std::cout << "lambdas: ";
        printVector(lambdas);
        double max = *std::max_element(lambdas.begin(), lambdas.end());
        double min = *std::min_element(lambdas.begin(), lambdas.end());
        if (max >= 0 || min >= 0) {
            curr_coeff = 0;
        } else {
            if (max == min) {
                curr_coeff = std::abs(max);
            } else {
                curr_coeff = std::max(max / min, min / max);
            }
        }
        coeff = std::max(coeff, curr_coeff);
        std::cout << "curr coeff: " << coeff << "\n";
    }
    return coeff;
}