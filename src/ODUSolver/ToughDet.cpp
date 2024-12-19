#include "ToughDet.hpp"

const uint64_t COUNT_OF_STEPS = 1000;
//std::vector<std::function<double (const std::vector<double> &)>>

double getToughK (const Matrix<double> &K) {
    double max = 0;
    
    for (uint64_t i = 0; i < K.size().n; ++i) {
        double tmp = 0;
        for (uint64_t j = 0; j < K.size().m; ++j) {
            tmp += std::abs(K(i, j));
        }
        max = std::max(max, tmp);
    }

    for (uint64_t i = 0; i < K.size().m; ++i) {
        double tmp = 0;
        for (uint64_t j = 0; j < K.size().n; ++j) {
            tmp += std::abs(K(j, i));
        }
        max = std::max(max, tmp);
    }

    return max;
}

double ToughCoeff (const std::vector<std::function<double (const std::vector<double> &)>> &system, const std::vector<double> &Y, const std::vector<double> &Yprev) {
    double coeff = 0, curr_coeff;
    uint64_t order = system.size();
    std::cout << "order: " << order << "\n";
    std::cout << "Y size: " << Y.size() << "\n";
    Matrix<double> A(order - 2, order - 2);
    double det = 0;
    for (uint64_t i = 0; i < order - 2; ++i) {
        for (uint64_t j = 0; j < order - 2; ++j) {
            //A(j, i) = derivative(system[j], Y, 0.00001, i + 1, DiffConfig::POINTS4_ORDER1_WAY2);
            A(j, i) = (system[j](Y) - system[j](Yprev)) / (Y[i + 1] - Yprev[i + 1]);
        }
    }
    //std::cout << "matrix Yakobi:\n" << A << "\n";
    //std::cout << "Yakobi:\n" << A << "\n";
    // Matrix<double> L, U, P;
    // double new_det = 1;
    // std::tie(L, U, P) = LU(A);
    // for (uint64_t i = 0; i < U.size().n; ++i) {
    //     new_det *= U(i, i);
    // }
    std::vector<double> lambdas = QRFindLambda(A, 0.001);
    double max = *std::max_element(lambdas.begin(), lambdas.end());
    double min = *std::min_element(lambdas.begin(), lambdas.end());
            //if (max >= 0 || min >= 0) {
            //    curr_coeff = 0;
            //} else {
                if (max == min) {
                    curr_coeff = std::abs(max);
                } else {
                    curr_coeff = std::max(std::abs(max / min), std::abs(min / max));
                }
            //}
            //coeff = std::max(coeff, curr_coeff);

    //det = std::max(det, -new_det);
//     A = Matrix<double>(6, 6, {-15.4673900864288,       -2696504.45348263,       -5.15579669547631,       -37.4167682121373,       -5.15579669547620,       -2164792.68116993,     
//    64.5219430333220,       -2696499.29768593,        914781.940065800,        64.5219430333220,        914781.934546957,        2164787.52537324,     
//   5.392998589853018E-014,   2696499.29768593,      -5.518842554199565E-003, -5.518842554194221E-023, -5.130852536653115E-015, -7.309394978856807E-015,
//   -32.2609715166610,       7.899968673056870E-014,  0.000000000000000E+000,  -32.2609715166610,       -914781.934546957,       7.899562598102051E-014,
//   -914761.311360175,        2696509.60927932,      -914771.617434723,        10.3115933909526,       -914771.622953566,        2164797.83696663,
//    -4.329575050746476E-014,  6.008593542900178E-014, -7.309394923668382E-015,  1.829563869093914E-014,   914781.934546957,       -2164787.52537324});
//     lambdas = QRFindLambda(A, double(0.001));
    //std::cout << "true matrix:\n" << A << "\n";
    //std::cout << "mult: " << LUGetDet(A) * LUGetDet(LUReverseMatrix(A)) << "\n";
    //std::cout << "mult: " << getToughK(A) * getToughK(LUReverseMatrix(A)) << "\n";
    std::cout << "lambdas: ";
    for (auto el : lambdas) {
        std::cout << el << " ";
    }
    std::cout << "\n";
    //printVector(lambdas);
    //std::cout << new_det << " " << det << "\n";
    //std::cout << A.det() << "\n";
    //coeff = std::max(coeff, curr_coeff);
    return curr_coeff;
    //return det;
}


// double ToughCoeff (const Task &task, const std::vector<double> &Y, double x) {
//     double coeff = 0, curr_coeff;
//     uint64_t order = task.order;
//     auto &system = task.odu_system;
//     Matrix<double> A(order, order);
//     std::vector<double> args(order + 1, 0);
//     for (uint64_t i = 1; i < args.size() - 1; ++i) {
//         args[i] = task.Y[i];
//     }
//     double step = (task.Xn - task.X0) / COUNT_OF_STEPS;
//     double det = 0;
//     for (double x = task.X0; x <= task.Xn; x += step) {
//         args[0] = x;
//         for (uint64_t i = 0; i < order; ++i) {
//             //args[i + 1] = 1;
//             for (uint64_t j = 0; j < order; ++j) {
//                 //A(j, i) = system[j](args);
//                 A(j, i) = derivative(system[j], args, 0.01, i + 1);
//             }
//             //args[i + 1] = 0;
//         }
//         Matrix<double> L, U, P;
//         double new_det = 1;
//         std::tie(L, U, P) = LU(A);
//         for (uint64_t i = 0; i < U.size().n; ++i) {
//             new_det *= U(i, i);
//         }
//         det = std::max(det, -new_det);
//         if (x == task.X0) {
//             std::cout << "matrix:\n" << A << "\n";
//             std::cout << new_det << " " << det << "\n";
//             std::cout << A.det() << "\n";
//         }
//         //det = std::max(det, std::abs(A.det()));
//         //std::cout << std::scientific << A(3, 0) << "\n";
//         // if (A(3, 0) == 0.0) {
//         //     std::cout << "good\n";
//         // } else {
//         //     std::cout << "bad\n";
//         // }
//         //std::cout << "\ncurrent A:\n" << A << "\n";
//         //std::vector<double> lambdas = QRFindLambda(A, 0.01);
//         //std::cout << "lambdas: ";
//         //printVector(lambdas);
//         //double max = *std::max_element(lambdas.begin(), lambdas.end());
//         //double min = *std::min_element(lambdas.begin(), lambdas.end());
//         //if (max >= 0 || min >= 0) {
//         //    curr_coeff = 0;
//         //} else {
//             // if (max == min) {
//             //     curr_coeff = std::abs(max);
//             // } else {
//             //     curr_coeff = std::max(max / min, min / max);
//             // }
//         //}
//         //coeff = std::max(coeff, curr_coeff);
//         //std::cout << "curr coeff: " << coeff << "\n";
//     }
//     //return coeff;
//     return det;
// }