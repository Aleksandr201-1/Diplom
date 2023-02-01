#include "RungeKutta.hpp"

//модуль определения жёсткости
//сохранение функционального дерева в файл
//добавить в метод рунге-кутты возможность вставлять критерий выбора шага
//итерационные методы (неявные, то есть с полной таблицей бутчера)
//печать отчёта
//графический интерфейс

std::pair<std::vector<double>, std::vector<double>> RungeKutta4 (const Task &task, double h) {
    double X0 = task.X0, Xn = task.Xn;
    const std::vector<FunctionalTree> &trees = task.trees;
    uint64_t order = task.order;
    //std::cout << task.Y.size() << " " << X0 << " " << Xn << "\n";

    std::vector<double> X, Y, Z;
    for (double i = X0; i <= Xn; i += h) {
        X.push_back(i);
    }
    Y.push_back(task.Y[0]);
    Z.push_back(task.Y[1]);
    auto g = [&] (double x, double y, double z) -> double {
        static auto c2 = trees[0].getCoeff(3);
        return -trees[0]({x, y, z, 0}) / c2(x);
    };
    //std::cout << Y[0] << " " << Z[0] << "\n";
    //std::vector<std::vector<double>> K(order, std::vector<double>(4));
    // std::vector<double> delta(order);
    for (uint64_t i = 1; i < X.size(); ++i) {
        // double delta;
        // for (uint64_t j = 0; j < order; ++j) {
        //     K[j][0] = h * Yi[j][i - 1];
        //     delta = 1.0 / 6.0 * (K[j][0] + 2 * K[j][1] + 2 * K[j][2] + K[j][3]);
        //     Yi[j].push_back(Yi[j][i - 1] + delta);
        // }
        double K[2][4];
        double delta_y, delta_z;
        double K1 = h * Z[i - 1];
        //double J1 = h * Z'[i - 1];
        double L1 = h * g(X[i - 1], Y[i - 1], Z[i - 1]);
        double K2 = h * (Z[i - 1] + 1.0 / 2.0 * L1);
        double L2 = h * g(X[i - 1] + 1.0 / 2.0 * h, Y[i - 1] + 1.0 / 2.0 * K1, Z[i - 1] + 1.0 / 2.0 * L1);
        double K3 = h * (Z[i - 1] + 1.0 / 2.0 * L2);
        double L3 = h * g(X[i - 1] + 1.0 / 2.0 * h, Y[i - 1] + 1.0 / 2.0 * K2, Z[i - 1] + 1.0 / 2.0 * L2);
        double K4 = h * (Z[i - 1] + L3);
        double L4 = h * g(X[i - 1] + h, Y[i - 1] + K3, Z[i - 1] + L3);
        //std::cout << K4 << " " << L4 << "\n";
        delta_y = 1.0 / 6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
        delta_z = 1.0 / 6.0 * (L1 + 2 * L2 + 2 * L3 + L4);
        Y.push_back(Y[i - 1] + delta_y);
        Z.push_back(Z[i - 1] + delta_z);
        K[0][0] = K1;
        K[1][0] = L1;
        K[0][1] = K2;
        K[1][1] = L2;
        K[0][2] = K3;
        K[1][2] = L3;
        K[0][3] = K4;
        K[1][3] = L4;
        // for (uint64_t j = 0; j < 4; ++j) {
        //     for (uint64_t k = 0; k < 2; ++k) {
        //         std::cout << "K[" << k << "][" << j << "] = " << K[k][j] << "\n";
        //     }
        // }
    }
    printVector(Z);
    return std::make_pair(X, Y);
}

double norma (const std::vector<double> &a, const std::vector<double> &b) {
    uint64_t size = std::min(a.size(), b.size());
    // if (a.size() != b.size()) {
    //     return 0;
    // }
    double ans = 0;
    for (uint64_t i = 0; i < a.size(); ++i) {
        ans = std::max(ans, std::abs(a[i] - b[i]));
    }
    return ans;
}

std::pair<std::vector<double>, std::vector<double>> RungeKutta6 (const Task &task, const Matrix<double> butcher, double h) {
    double X0 = task.X0, Xn = task.Xn;
    const std::vector<FunctionalTree> &trees = task.trees;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    // for (double i = X0; i <= Xn; i += h) {
    //     Yi[0].push_back(i);
    // }
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }
    auto g = [&] (const std::vector<double> &args) -> double {
        static auto c2 = trees[0].getCoeff(orderOfTask + 1);
        return -trees[0](args) / c2(args[0]);
    };
    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    uint64_t i = 1;
    //for (uint64_t i = 1; i < Yi[0].size(); ++i) {
    while (Yi[0].back() + h < Xn) {
        //std::cout << i << "\n";

        //начальная инициализация
        // for (uint64_t j = 0; j < args.size() - 1; ++j) {
        //     args[j] = Yi[j][i - 1];
        // }
        // //args[0] = Yi[0][i - 1];
        // for (uint64_t j = 0; j < orderOfTask - 1; ++j) {
        //     //args[j + 1] = Yi[j + 1][i - 1];
        //     K[j][0] = h * Yi[orderOfTask - j][i - 1];
        // }
        // //args[orderOfTask] = Yi[orderOfTask][i - 1];
        // K[orderOfTask - 1][0] = h * g(args);

        for (uint64_t j = 0; j < orderOfApprox; ++j) {

            //аргументы для функций
            args[0] = Yi[0][i - 1] + butcher(j, 0) * h;
            for (uint64_t k = 1; k < args.size() - 1; ++k) {
                args[k] = Yi[k][i - 1];
                for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                    args[k] += butcher(j, l + 1) * K[k - 1][l];
                }
            }

            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask - 1; ++k) {
                //args[k + 1] = Yi[k][i - 1];
                K[k][j] = Yi[orderOfTask - k][i - 1];// + butcher(j, k + 1) * K[orderOfTask - k - 1][j - 1];
                for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                    // K[k][j] += butcher(j, l + 1) * K[orderOfTask - 1 - k][j - 1];
                    K[k][j] += butcher(j, l + 1) * K[orderOfTask - 1 - k][l]; //???
                }
                K[k][j] *= h;
            }

            //аргументы
            // args[0] = Yi[0][i - 1] + butcher(j, 0) * h;
            // for (uint64_t k = 1; k < args.size() - 1; ++k) {
            //     args[k] = Yi[k][i - 1];
            //     for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
            //         args[k] += butcher(j, l + 1) * K[k - 1][l];
            //     }
            // }
            K[orderOfTask - 1][j] = h * g(args);
        }
        // for (uint64_t j = 0; j < orderOfApprox; ++j) {
        //     for (uint64_t k = 0; k < orderOfTask; ++k) {
        //         std::cout << "K[" << k << "][" << j << "] = " << K[k][j] << "\n";
        //     }
        // }

        // K[0][0] = h * Yi[2][i - 1]; //k1
        // K[1][0] = h * g({Yi[0][i - 1], Yi[1][i - 1], Yi[2][i - 1], 0}); //l1
        // K[0][1] = h * (Yi[2][i - 1] + 1.0 / 2.0 * K[1][0]); //k2
        // K[1][1] = h * g({Yi[0][i - 1] + 1.0 / 2.0 * h, Yi[1][i - 1] + 1.0 / 2.0 * K[0][0], Yi[2][i - 1] + 1.0 / 2.0 * K[1][0], 0}); //l2
        // K[0][2] = h * (Yi[2][i - 1] + 1.0 / 2.0 * K[1][1]); //k3
        // K[1][2] = h * g({Yi[0][i - 1] + 1.0 / 2.0 * h, Yi[1][i - 1] + 1.0 / 2.0 * K[0][1], Yi[2][i - 1] + 1.0 / 2.0 * K[1][1], 0}); //l3
        // K[0][3] = h * (Yi[2][i - 1] + K[1][2]); //k4
        // K[1][3] = h * g({Yi[0][i - 1] + h, Yi[1][i - 1] + K[0][2], Yi[2][i - 1] + K[1][2], 0}); //l4
        //std::cout << K[0][3] << " " << K[1][3] << "\n";

        //добавление нового значения
        Yi[0].push_back(Yi[0][i - 1] + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1][i - 1] + delta);
        }

        //изменение шага
        //std::cout << K[0][1] << " " << K[0][2] << " " << K[0][0] << " " << K[0][1] << "\n";
        if (!isEqual(K[0][1], K[0][2]) && !isEqual(K[0][0], K[0][1])) {
            double Q = std::abs(K[0][1] - K[0][2]) / std::abs(K[0][0] - K[0][1]);
            //std::cout << K[0][1] << " " << K[0][2] << " " << K[0][0] << " " << K[0][1] << "\n";
            // if (i == 4) {
            //     exit(-1);
            // }
            if (Q > 0.6) { //0.1
                h /= 2;
            } else if (Q < 0.01) {
                h *= 2;
            }
            std::cout << "Q = " << Q << "\nh = " << h << "\n";
        }
        std::cout << "X = " << Yi[0].back() << "\n";
        ++i;
        //printVector(Yi[1]);
    }
    //printVector(Yi[2]);
    return std::make_pair(Yi[0], Yi[1]);
}

std::pair<std::vector<double>, std::vector<double>> Falberg (const Task &task, const Matrix<double> butcher, double h) {
    double X0 = task.X0, Xn = task.Xn;
    const std::vector<FunctionalTree> &trees = task.trees;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    std::vector<double> control;
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }
    control.push_back(task.Y[0]);
    auto g = [&] (const std::vector<double> &args) -> double {
        static auto c2 = trees[0].getCoeff(orderOfTask + 1);
        return -trees[0](args) / c2(args[0]);
    };
    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    uint64_t i = 1;
    while (Yi[0].back() + h < Xn) {

        for (uint64_t j = 0; j < orderOfApprox; ++j) {

            //аргументы для функций
            args[0] = Yi[0][i - 1] + butcher(j, 0) * h;
            for (uint64_t k = 1; k < args.size() - 1; ++k) {
                args[k] = Yi[k][i - 1];
                for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                    args[k] += butcher(j, l + 1) * K[k - 1][l];
                }
            }

            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask - 1; ++k) {
                K[k][j] = Yi[orderOfTask - k][i - 1];
                for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                    K[k][j] += butcher(j, l + 1) * K[orderOfTask - 1 - k][l]; //???
                }
                K[k][j] *= h;
            }

            K[orderOfTask - 1][j] = h * g(args);
        }

        //добавление нового значения
        Yi[0].push_back(Yi[0][i - 1] + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1][i - 1] + delta);
        }
        
        double delta_control = 0;
        for (uint64_t k = 1; k <= orderOfApprox; ++k) {
            delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
        }
        control.push_back(control[i - 1] + delta_control);

        //изменение шага
        double R = norma(Yi[1], control);
        double eps = 0.1;
        if (R > eps) {
            h /= 2;
        } else if (R < eps / 64) {
            h *= 2;
        }
        // if (!isEqual(K[0][1], K[0][2]) && !isEqual(K[0][0], K[0][1])) {
        //     double Q = std::abs(K[0][1] - K[0][2]) / std::abs(K[0][0] - K[0][1]);
        //     if (Q > 0.1) {
        //         h /= 2;
        //     } else if (Q < 0.01) {
        //         h *= 2;
        //     }
        // }
        ++i;
        std::cout << R << "\n";
    }
    return std::make_pair(Yi[0], Yi[1]);
}

// void Zeidel (std::vector<std::vector<double>> &K, const std::vector<std::vector<double>> &Yi, const std::vector<std::function<double(const std::vector<double> &)>> &func, const Matrix<double> butcher, double h) {
//     std::vector<double> args(Yi.size() + 1, 0);
//     for (uint64_t i = 0; i < K.size(); ++i) {
//         for (uint64_t j = 0; j < K[0].size(); ++j) {
//             for (uint64_t k = 0; k < args.size(); ++k) {
//                 args[k] = Yi[k][0];
//             }
//         }
//     }
// }

// void Newton (std::vector<std::vector<double>> &K, const std::vector<std::vector<double>> &Yi, const std::vector<std::function<double(const std::vector<double> &)>> &func, const Matrix<double> butcher, double h) {
    
// }

std::pair<std::vector<double>, std::vector<double>> NonExplZeidel (const Task &task, const Matrix<double> butcher, double h) {
    double X0 = task.X0, Xn = task.Xn;
    const std::vector<FunctionalTree> &trees = task.trees;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    //std::vector<double> control;
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }
    //control.push_back(task.Y[0]);
    auto g = [&] (const std::vector<double> &args) -> double {
        static auto c2 = trees[0].getCoeff(orderOfTask + 1);
        return -trees[0](args) / c2(args[0]);
    };
    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    uint64_t i = 1;
    while (Yi[0].back() + h < Xn) {

        std::vector<double> Ktmp;
        do {
            for (uint64_t j = 0; j < orderOfApprox; ++j) {
                Ktmp = K[orderOfTask - 1];

                //аргументы для функций
                args[0] = Yi[0][i - 1] + butcher(j, 0) * h;
                for (uint64_t k = 1; k < args.size() - 1; ++k) {
                    args[k] = Yi[k][i - 1];
                    for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                        args[k] += butcher(j, l + 1) * K[k - 1][l];
                    }
                }

                //коэффициенты
                for (uint64_t k = 0; k < orderOfTask - 1; ++k) {
                    K[k][j] = Yi[orderOfTask - k][i - 1];
                    for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                        K[k][j] += butcher(j, l + 1) * K[orderOfTask - 1 - k][l]; //???
                    }
                    K[k][j] *= h;
                }

                K[orderOfTask - 1][j] = h * g(args);
                std::cout << "norma: " << norma(Ktmp, K[orderOfTask - 1]) << "\n";
            }
        } while (norma(Ktmp, K[orderOfTask - 1]) > 0.01);

        //добавление нового значения
        Yi[0].push_back(Yi[0][i - 1] + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1][i - 1] + delta);
        }

        // //изменение шага
        // double R = norma(Yi[1], control);
        // double eps = 0.1;
        // if (R > eps) {
        //     h /= 2;
        // } else if (R < eps / 64) {
        //     h *= 2;
        // }
        ++i;
        // std::cout << R << "\n";
    }
    return std::make_pair(Yi[0], Yi[1]);
}

std::pair<std::vector<double>, std::vector<double>> NonExplNewton (const Task &task, const Matrix<double> butcher, double h) {
    double X0 = task.X0, Xn = task.Xn;
    const std::vector<FunctionalTree> &trees = task.trees;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    //std::vector<double> control;
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }
    //control.push_back(task.Y[0]);
    auto g = [&] (const std::vector<double> &args) -> double {
        static auto c2 = trees[0].getCoeff(orderOfTask + 1);
        return -trees[0](args) / c2(args[0]);
    };
    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    uint64_t i = 1;
    while (Yi[0].back() + h < Xn) {

        std::vector<double> Ktmp;
        do {
            for (uint64_t j = 0; j < orderOfApprox; ++j) {
                Ktmp = K[orderOfTask - 1];

                //аргументы для функций
                args[0] = Yi[0][i - 1] + butcher(j, 0) * h;
                for (uint64_t k = 1; k < args.size() - 1; ++k) {
                    args[k] = Yi[k][i - 1];
                    for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                        args[k] += butcher(j, l + 1) * K[k - 1][l];
                    }
                }

                //коэффициенты
                for (uint64_t k = 0; k < orderOfTask - 1; ++k) {
                    K[k][j] = Yi[orderOfTask - k][i - 1];
                    for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                        K[k][j] += butcher(j, l + 1) * K[orderOfTask - 1 - k][l]; //???
                    }
                    K[k][j] *= h;
                }

                K[orderOfTask - 1][j] = h * g(args);
                std::cout << "norma: " << norma(Ktmp, K[orderOfTask - 1]) << "\n";
            }
        } while (norma(Ktmp, K[orderOfTask - 1]) > 0.01);

        //добавление нового значения
        Yi[0].push_back(Yi[0][i - 1] + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1][i - 1] + delta);
        }

        // //изменение шага
        // double R = norma(Yi[1], control);
        // double eps = 0.1;
        // if (R > eps) {
        //     h /= 2;
        // } else if (R < eps / 64) {
        //     h *= 2;
        // }
        ++i;
        // std::cout << R << "\n";
    }
    return std::make_pair(Yi[0], Yi[1]);
}