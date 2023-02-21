#include "RungeKutta.hpp"

//модуль определения жёсткости
//сохранение функционального дерева в файл
//добавить в метод рунге-кутты возможность вставлять критерий выбора шага
//итерационные методы (неявные, то есть с полной таблицей бутчера)
//печать отчёта
//графический интерфейс

double norma (const std::vector<double> &a, const std::vector<double> &b) {
    uint64_t size = std::min(a.size(), b.size());
    double ans = 0;
    for (uint64_t i = 0; i < a.size(); ++i) {
        ans = std::max(ans, std::abs(a[i] - b[i]));
    }
    return ans;
}

std::pair<std::vector<double>, std::vector<double>> RungeKutta6 (const Task &task, const Matrix<double> butcher, double h) {
    const double &X0 = task.X0, &Xn = task.Xn;
    auto &func = task.odu_system;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }

    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    //uint64_t i = 1;

    //for (uint64_t i = 1; i < Yi[0].size(); ++i) {
    while (true) {
        for (uint64_t j = 0; j < orderOfApprox; ++j) {

            //аргументы для функций
            args[0] = Yi[0].back() + butcher(j, 0) * h;
            for (uint64_t k = 1; k < args.size() - 1; ++k) {
                args[k] = Yi[k].back();
                for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                    args[k] += butcher(j, l + 1) * K[k - 1][l];
                }
            }

            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask; ++k) {
                K[k][j] = h * func[k](args);
            }
        }

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

        //изменение шага
        if (!isEqual(K[0][1], K[0][2]) && !isEqual(K[0][0], K[0][1])) {
            double Q = 0;
            for (uint64_t j = 0; j < K.size(); ++j) {
                if (!isEqual(K[j][1], K[j][2]) && !isEqual(K[j][0], K[j][1])) {
                    double tmp = std::abs(K[j][1] - K[j][2]) / std::abs(K[j][0] - K[j][1]);
                    Q = std::max(Q, tmp);
                }
            }
            //std::cout << K[0][1] << " " << K[0][2] << " " << K[0][0] << " " << K[0][1] << "\n";
            if (Q > 0.1 && h > min_h) { //0.1
                h /= 2;
            } else if (Q < 0.01) {
                h *= 2;
            }
            // if (Yi[0].back() + h > Xn && Yi[0].back() != Xn) {
            //     h = Xn - Yi[0].back();
            // }
            //std::cout << "Q = " << Q << "\nh = " << h << "\n";
        }
        //изменение шага
        if (Yi[0].back() + h > Xn) {
            if (Xn - Yi[0].back() > min_h) {
                h = Xn - Yi[0].back();
            } else {
                break;
            }
        }
        std::cout << "X = " << Yi[0].back() << "\n";
        std::cout << "h = " << h << "\n";
        //++i;
    }
    return std::make_pair(Yi[0], Yi[1]);
}

std::pair<std::vector<double>, std::vector<double>> Falberg (const Task &task, const Matrix<double> butcher, double h) {
    const double &X0 = task.X0, &Xn = task.Xn;
    auto &func = task.odu_system;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::cout << "R\n" << orderOfApprox << "\n";
    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    std::vector<double> control;
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }
    control.push_back(task.Y[0]);

    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    //uint64_t i = 1;
    std::cout << "R\n";

    while (true) {

        for (uint64_t j = 0; j < orderOfApprox; ++j) {

            //аргументы для функций
            args[0] = Yi[0].back() + butcher(j, 0) * h;
            for (uint64_t k = 1; k < args.size() - 1; ++k) {
                args[k] = Yi[k].back();
                for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                    args[k] += butcher(j, l + 1) * K[k - 1][l];
                }
            }

            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask; ++k) {
                K[k][j] = h * func[k](args);
            }
        }

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }
        
        //изменение шага
        double delta_control = 0;
        for (uint64_t k = 1; k <= orderOfApprox; ++k) {
            delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
        }
        control.push_back(control.back() + delta_control);
        double R = norma(Yi[1], control);
        double eps = 0.1;
        if (R > eps && h > min_h) {
            h /= 2;
        } else if (R < eps / 64) {
            h *= 2;
        }
        if (Yi[0].back() + h > Xn) {
            if (Xn - Yi[0].back() > min_h) {
                h = Xn - Yi[0].back();
            } else {
                break;
            }
        }
        // if (Yi[0].back() + h > Xn && Yi[0].back() != Xn) {
        //     h = Xn - Yi[0].back();
        // }
        //++i;
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
    auto &func = task.odu_system;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y'' ...
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }

    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);

    for (uint64_t i = 0; i < args.size() - 1; ++i) {
        args[i] = Yi[i][0];
    }
    for (uint64_t i = 0; i < K[0].size(); ++i) {
        for (uint64_t j = 0; j < K.size(); ++j) {
            //K[i][j] = butcher(0, j + 1);// * Yi[i + 1][j];
            K[j][i] = func[j](args);
        }
    }

    while (true) {
        std::vector<double> Ktmp;
        auto KK = K;
        double nor = 0;
        std::cout << "func: " << func[0]({0, 1, 1}) << "\n";
        do {
            Ktmp = K[orderOfTask - 1];
            KK = K;
            nor = 0;
            for (uint64_t j = 0; j < orderOfApprox; ++j) {

                //аргументы для функций
                args[0] = Yi[0].back() + butcher(j, 0) * h;
                for (uint64_t k = 1; k < args.size() - 1; ++k) {
                    args[k] = Yi[k].back();
                    for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                        args[k] += butcher(j, l + 1) * K[k - 1][l] * h;
                    }
                }

                //коэффициенты
                for (uint64_t k = 0; k < orderOfTask; ++k) {
                    K[k][j] = func[k](args);
                }
                std::cout << "norma: " << norma(Ktmp, K[orderOfTask - 1]) << "\n";
                // std::cout << "X: " << Yi[0].back() << "\n";
                // std::cout << "h: " << h << "\n";
            }
            for (uint64_t i = 0; i < K.size(); ++i) {
                nor = std::max(nor, norma(KK[i], K[i]));
            }
        //} while (norma(Ktmp, K[orderOfTask - 1]) > 0.001);
        } while (nor > 0.01);

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1] * h;
            }
            std::cout << "delta: " << delta << "\n";
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

        //изменение шага
        if (Yi[0].back() + h > Xn) {
            if (Xn - Yi[0].back() > min_h) {
                h = Xn - Yi[0].back();
            } else {
                break;
            }
        }
    }
    return std::make_pair(Yi[0], Yi[1]);
}

std::pair<std::vector<double>, std::vector<double>> NonExplNewton (const Task &task, const Matrix<double> butcher, double h) {
    double X0 = task.X0, Xn = task.Xn;
    auto &func = task.odu_system;
    uint64_t orderOfTask = task.order, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y''...
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(task.Y[i]);
    }

    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));
    std::vector<double> args(orderOfTask + 2, 0);
    //uint64_t i = 1;
    for (uint64_t i = 0; i < K.size(); ++i) {
        for (uint64_t j = 0; j < K[0].size(); ++j) {
            K[i][j] = butcher(0, j + 1) * Yi[i + 1][j];
        }
    }

    while (Yi[0].back() + h <= Xn) {

        std::vector<double> Ktmp;
        do {
            Ktmp = K[orderOfTask - 1];
            for (uint64_t j = 0; j < orderOfApprox; ++j) {

                //аргументы для функций
                args[0] = Yi[0].back() + butcher(j, 0) * h;
                for (uint64_t k = 1; k < args.size() - 1; ++k) {
                    args[k] = Yi[k].back();
                    for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
                        args[k] += butcher(j, l + 1) * K[k - 1][l] * h;
                    }
                }

                //коэффициенты
                for (uint64_t k = 0; k < orderOfTask; ++k) {
                    K[k][j] = func[k](args);
                }
                std::cout << "norma: " << norma(Ktmp, K[orderOfTask - 1]) << "\n";
            }
        } while (norma(Ktmp, K[orderOfTask - 1]) > 0.01);

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1] * h;
            }
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

        //изменение шага
        if (Yi[0].back() + h > Xn && Yi[0].back() != Xn) {
            h = Xn - Yi[0].back();
        }
        //++i;
    }
    return std::make_pair(Yi[0], Yi[1]);
}

std::pair<std::vector<double>, std::vector<double>> ODUSolve (SolveMethod method, const Task &task, const Matrix<double> butcher, double h, IterationAlgo iter_alg, double approx) {
    switch (method) {
        case SolveMethod::RUNGE_KUTTA:
            return RungeKutta6(task, butcher, h);
        case SolveMethod::FALBERG:
        case SolveMethod::CHESKINO:
        case SolveMethod::MERSON:
            return Falberg(task, butcher, h);
        default:
            return NonExplZeidel(task, butcher, h);
    }
}