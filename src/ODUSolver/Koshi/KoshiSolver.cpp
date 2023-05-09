#include "KoshiSolver.hpp"

//модуль определения жёсткости
//сохранение функционального дерева в файл
//добавить в метод рунге-кутты возможность вставлять критерий выбора шага
//итерационные методы (неявные, то есть с полной таблицей бутчера)
//печать отчёта
//графический интерфейс

std::vector<std::vector<float128_t>> RungeKutta (const KoshiTask &task, const Matrix<float128_t> butcher, float128_t h) {
    float128_t X0, Xn;
    std::tie(X0, Xn) = task.getBorders();
    const auto &func = task.getODE();
    const auto &Y0 = task.getY0();
    uint64_t orderOfTask = func.size(), orderOfApprox = butcher.size().m - 1;
    float128_t min_h = (Xn - X0) / MIN_H_SIZE;
    float128_t max_h = (Xn - X0) / MAX_H_SIZE;
    float128_t tough = 0;

    std::vector<std::vector<float128_t>> Yi(orderOfTask + 1); //X Y Y' Y''...

    //инициализация начальных условий
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(Y0[i]);
    }

    std::vector<std::vector<float128_t>> K(orderOfTask, std::vector<float128_t>(orderOfApprox, 0));
    //std::vector<float128_t> args(orderOfTask + 2, 0);
    //uint64_t i = 1;
    uint64_t stepNum = 1;
    //std::cout << Yi.size() << " " << Yi[0].size() << " " << orderOfTask << " " << orderOfApprox << "\n";
    //exit(0);
    //for (uint64_t i = 1; i < Yi[0].size(); ++i) {
    while (true) {
        // for (uint64_t j = 0; j < orderOfApprox; ++j) {

        //     //аргументы для функций
        //     args[0] = Yi[0].back() + butcher(j, 0) * h;
        //     for (uint64_t k = 1; k < args.size() - 1; ++k) {
        //         args[k] = Yi[k].back();
        //         for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
        //             args[k] += butcher(j, l + 1) * K[k - 1][l];
        //         }
        //     }

        //     //коэффициенты
        //     for (uint64_t k = 0; k < orderOfTask; ++k) {
        //         K[k][j] = h * func[k](args);
        //     }
        // }
        K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
        //tough = std::max(tough, getToughK(K));

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            float128_t delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

        //изменение и корректировка шага (Рунге-Ромберг)
        // if (Yi[0].size() > 1) {
        //     std::vector<float128_t> tmpY2h(Yi.size());
        //     for (uint64_t i = 0; i < Yi.size(); ++i) {
        //         tmpY2h[i] = Yi[i][stepNum];
        //     }
        //     float128_t R = 0, tmpR;
        //     //float128_t h2Y;
        //     float128_t half_step = h / 2;
        //     for (uint64_t i = 0; i < 2; ++i) {}
        //     K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, half_step);
        //     //K = ExplicitStep(K, args, Yi, stepNum - 1, func, butcher, half_step);
        //     //tmpY1h[0] = Yi[0][stepNum - 1] + half_step;
        //     Yi[0][stepNum] = Yi[0][stepNum - 1] + half_step;
        //     for (uint64_t j = 0; j < orderOfTask; ++j) {
        //         float128_t delta = 0;
        //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
        //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
        //         }
        //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        //         Yi[j + 1][stepNum] = Yi[j + 1][stepNum - 1] + delta;
        //     }

        //     K = ExplicitStep(K, Yi, stepNum, func, butcher, half_step);
        //     Yi[0][stepNum] = Yi[0][stepNum] + half_step;
        //     for (uint64_t j = 0; j < orderOfTask; ++j) {
        //         float128_t delta = 0;
        //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
        //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
        //         }
        //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        //         Yi[j + 1][stepNum] = Yi[j + 1][stepNum] + delta;
        //     }
        //     // for (uint64_t i = 0; i < Yi.size(); ++i) {
        //     //     Yi[i][stepNum] = tmpY1h[i];
        //     // }
        //     for (uint64_t j = 0; j < Yi.size(); ++j) {
        //         tmpR = (Yi[j][stepNum] - tmpY2h[j]) / (std::pow(2, orderOfApprox) - 1);
        //         Yi[j][stepNum] += tmpR;
        //         R = std::max(R, std::abs(tmpR));
        //     }
        //     //std::cout << "R = " << R << "\n";
        //     if (R > 0.01 && h / 2 >= min_h) {
        //         h /= 2;
        //     } else if (R < 0.001 && h * 2 <= max_h) {
        //         //h *= 2;
        //     }
        // }
        //
        // if (Yi[0].size() > 2) {
        //     std::vector<float128_t> tmpY(Yi.size());
        //     float128_t R = 0, tmpR;
        //     //float128_t h2Y;
        //     K = ExplicitStep(K, Yi, stepNum - 2, func, butcher, h*2);
        //     tmpY[0] = Yi[0][stepNum - 2] + 2 * h;
        //     for (uint64_t j = 0; j < orderOfTask; ++j) {
        //         float128_t delta = 0;
        //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
        //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
        //         }
        //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        //         tmpY[j + 1] = Yi[j + 1][stepNum - 2] + delta;
        //     }

        //     for (uint64_t j = 0; j < Yi.size(); ++j) {
        //         tmpR = (Yi[j][stepNum - 1] - tmpY[j]) / (std::pow(2, orderOfApprox) - 1);
        //         Yi[j][stepNum - 1] += tmpR;
        //         R = std::max(R, std::abs(tmpR));
        //     }

        //     if (R > 0.01 && h / 2 > min_h) {
        //         h /= 2;
        //     } else if (R < 0.001) {
        //         //h *= 2;
        //     }
        // }

        //изменение шага
        // if (!isEqual(K[0][1], K[0][2]) && !isEqual(K[0][0], K[0][1])) {
        //     float128_t Q = 0;
        //     for (uint64_t j = 0; j < K.size(); ++j) {
        //         if (!isEqual(K[j][1], K[j][2]) && !isEqual(K[j][0], K[j][1])) {
        //             float128_t tmp = std::abs(K[j][1] - K[j][2]) / std::abs(K[j][0] - K[j][1]);
        //             Q = std::max(Q, tmp);
        //         }
        //     }
        //     //std::cout << K[0][1] << " " << K[0][2] << " " << K[0][0] << " " << K[0][1] << "\n";
        //     if (Q > 0.1 && h > min_h) { //0.1
        //         h /= 2;
        //     } else if (Q < 0.01) {
        //         h *= 2;
        //     }
        //     // if (Yi[0].back() + h > Xn && Yi[0].back() != Xn) {
        //     //     h = Xn - Yi[0].back();
        //     // }
        //     //std::cout << "Q = " << Q << "\nh = " << h << "\n";
        // }

        //изменение шага
        if (Yi[0].back() + h > Xn) {
            if (Xn - Yi[0].back() >= min_h) {
                h = Xn - Yi[0].back();
            } else {
                break;
            }
        }
        //std::cout << "X = " << Yi[0].back() << "\n";
        //std::cout << "h = " << h << "\n";
        //++i;
        ++stepNum;
    }
    //return std::make_pair(Yi[0], Yi[1]);
    return Yi;
}

// std::vector<std::vector<float128_t>> Falberg (const KoshiTask &task, const Matrix<float128_t> butcher, float128_t h) {
//     const float128_t &X0 = task.X0, &Xn = task.Xn;
//     auto &func = task.odu_system;
//     uint64_t orderOfTask = func.size(), orderOfApprox = butcher.size().m - 1;
//     float128_t min_h = (Xn - X0) / MIN_H_SIZE;
//     float128_t tough = 0;

//     std::cout << "R\n" << orderOfApprox << "\n";
//     std::vector<std::vector<float128_t>> Yi(orderOfTask + 1); //X Y Y' Y''...
//     std::vector<float128_t> control;

//     //инициализация начальных условий
//     Yi[0].push_back(X0);
//     for (uint64_t i = 0; i < orderOfTask; ++i) {
//         Yi[i + 1].push_back(task.Y[i]);
//     }
//     control.push_back(task.Y[0]);

//     std::vector<std::vector<float128_t>> K(orderOfTask, std::vector<float128_t>(orderOfApprox, 0));
//     std::vector<float128_t> args(orderOfTask + 2, 0);
//     //uint64_t i = 1;
//     std::cout << "R\n";

//     while (true) {

//         for (uint64_t j = 0; j < orderOfApprox; ++j) {

//             //аргументы для функций
//             args[0] = Yi[0].back() + butcher(j, 0) * h;
//             for (uint64_t k = 1; k < args.size() - 1; ++k) {
//                 args[k] = Yi[k].back();
//                 for (uint64_t l = 0; l < orderOfApprox; ++l) { //change orderOfApprox to j?
//                     args[k] += butcher(j, l + 1) * K[k - 1][l];
//                 }
//             }

//             //коэффициенты
//             for (uint64_t k = 0; k < orderOfTask; ++k) {
//                 K[k][j] = h * func[k](args);
//             }
//         }
//         tough = std::max(tough, getToughK(K));

//         //добавление нового значения
//         Yi[0].push_back(Yi[0].back() + h);
//         for (uint64_t j = 0; j < orderOfTask; ++j) {
//             float128_t delta = 0;
//             for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//                 delta += butcher(orderOfApprox, k) * K[j][k - 1];
//             }
//             Yi[j + 1].push_back(Yi[j + 1].back() + delta);
//         }
        
//         //изменение шага
//         float128_t delta_control = 0;
//         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//             delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
//         }
//         control.push_back(control.back() + delta_control);
//         float128_t R = std::abs(Yi[1].back() - control.back()); //norma(Yi[1], control);
//         float128_t eps = 0.1;
//         if (R > eps && h > min_h) {
//             h /= 2;
//         } else if (R < eps / 64) {
//             h *= 2;
//         }
//         if (Yi[0].back() + h > Xn) {
//             if (Xn - Yi[0].back() > min_h) {
//                 h = Xn - Yi[0].back();
//             } else {
//                 break;
//             }
//         }
//         // if (Yi[0].back() + h > Xn && Yi[0].back() != Xn) {
//         //     h = Xn - Yi[0].back();
//         // }
//         //++i;
//         std::cout << R << "\n";
//     }
//     //return std::make_pair(Yi[0], Yi[1]);
//     return Yi;
// }

// void Zeidel (std::vector<std::vector<float128_t>> &K, const std::vector<std::vector<float128_t>> &Yi, const std::vector<std::function<float128_t(const std::vector<float128_t> &)>> &func, const Matrix<float128_t> butcher, float128_t h) {
//     std::vector<float128_t> args(Yi.size() + 1, 0);
//     for (uint64_t i = 0; i < K.size(); ++i) {
//         for (uint64_t j = 0; j < K[0].size(); ++j) {
//             for (uint64_t k = 0; k < args.size(); ++k) {
//                 args[k] = Yi[k][0];
//             }
//         }
//     }
// }

// void Newton (std::vector<std::vector<float128_t>> &K, const std::vector<std::vector<float128_t>> &Yi, const std::vector<std::function<float128_t(const std::vector<float128_t> &)>> &func, const Matrix<float128_t> butcher, float128_t h) {
    
// }

// std::vector<std::vector<float128_t>> NonExpl (const ChemicalSystem &task, const Matrix<float128_t> butcher, float128_t h, IterationAlgo iter_alg, float128_t approx) {
//     float128_t X0 = 0, Xn = h * 200;
//     //std::tie(X0, Xn) = task.getBorders();
//     auto func = task.getODE();
//     auto Y = task.getY0();
//     uint64_t orderOfTask = func.size(), orderOfApprox = butcher.size().m - 1;
//     float128_t min_h = (Xn - X0) / MIN_H_SIZE;
//     float128_t max_h = (Xn - X0) / MAX_H_SIZE;
//     float128_t tough = 0;

//     std::vector<std::vector<float128_t>> Yi(orderOfTask + 1); //X Y Y' Y'' ...
//     std::vector<float128_t> control;

//     //инициализация начальных условий
//     Yi[0].push_back(X0);
//     for (uint64_t i = 0; i < orderOfTask; ++i) {
//         Yi[i + 1].push_back(Y[i]);
//     }
//     control.push_back(Y[0]);

//     //order of approx
//     //----->
//     //K1 K2 K3 K4 | order of task
//     //L1 L2 L3 L4 |
//     //M1 M2 M3 M4 V
//     std::vector<std::vector<float128_t>> K(orderOfTask, std::vector<float128_t>(orderOfApprox, 0));
//     //std::vector<float128_t> args(orderOfTask + 2, 0);

//     uint64_t stepNum = 1;

//     while (true) {
//         //auto KK = K;
//         //float128_t nor = 0;
//         //uint64_t iter = 0;
//         //std::cout << "func: " << func[0]({0, 1, 1}) << "\n";
//         //if (Yi.size() < 2) {
//             //for (uint64_t i = 0; i < args.size() - 1; ++i) {
//             //    args[i] = Yi[i].back();
//             //}
//             //for (uint64_t i = 0; i < K[0].size(); ++i) {
//                 for (uint64_t j = 0; j < K.size(); ++j) {
//                     //K[i][j] = butcher(0, j + 1) * Yi[i + 1][j];
//                     //K[j][0] = func[j](args);
//                     //K[j][i] = 0;
//                     //K[j][i] = butcher(0, j + 1);
//                 }
//             //}
//         //}
//         if (1) {
//             if (stepNum == 1) {
//                 K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
//             }
//             K = IterationStep(K, Yi, stepNum - 1, func, butcher, h, approx, iter_alg);
//         } else {
//             K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
//         }
//         //tough = std::max(tough, getToughK(K));

//         //добавление нового значения
//         Yi[0].push_back(Yi[0].back() + h);
//         for (uint64_t j = 0; j < orderOfTask; ++j) {
//             float128_t delta = 0;
//             for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//                 delta += butcher(orderOfApprox, k) * K[j][k - 1] * h;
//             }
//             //std::cout << "delta: " << delta << "\n";
//             Yi[j + 1].push_back(Yi[j + 1].back() + delta);
//         }
//         // std::vector<float128_t> Y(Yi.size());
//         // for (uint64_t j = 0; j < Yi.size(); ++j) {
//         //     Y[j] = Yi[j].back();
//         // }
//         //tough = std::max(tough, ToughCoeff(func, Y));
//         //std::cout << "Tough " << stepNum << ": " << tough << "\n";

//         // if (Yi[0].size() > 1) {
//         //     std::vector<float128_t> tmpY2h(Yi.size());
//         //     for (uint64_t i = 0; i < Yi.size(); ++i) {
//         //         tmpY2h[i] = Yi[i][stepNum];
//         //     }
//         //     float128_t R = 0, tmpR;
//         //     //float128_t h2Y;
//         //     float128_t half_step = h / 2;
//         //     for (uint64_t i = 0; i < 2; ++i) {}

//         //     K = IterationStep(K, Yi, stepNum - 1, func, butcher, half_step, approx, iter_alg);
//         //     Yi[0][stepNum] = Yi[0][stepNum - 1] + half_step;
//         //     for (uint64_t j = 0; j < orderOfTask; ++j) {
//         //         float128_t delta = 0;
//         //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//         //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
//         //         }
//         //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
//         //         Yi[j + 1][stepNum] = Yi[j + 1][stepNum - 1] + delta;
//         //     }

//         //     //K = ExplicitStep(K, Yi, stepNum, func, butcher, half_step);
//         //     K = IterationStep(K, Yi, stepNum, func, butcher, half_step, approx, iter_alg);
//         //     Yi[0][stepNum] = Yi[0][stepNum] + half_step;
//         //     for (uint64_t j = 0; j < orderOfTask; ++j) {
//         //         float128_t delta = 0;
//         //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//         //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
//         //         }
//         //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
//         //         Yi[j + 1][stepNum] = Yi[j + 1][stepNum] + delta;
//         //     }
//         //     // for (uint64_t i = 0; i < Yi.size(); ++i) {
//         //     //     Yi[i][stepNum] = tmpY1h[i];
//         //     // }
//         //     for (uint64_t j = 0; j < Yi.size(); ++j) {
//         //         tmpR = (Yi[j][stepNum] - tmpY2h[j]) / (std::pow(2, orderOfApprox) - 1);
//         //         Yi[j][stepNum] += tmpR;
//         //         R = std::max(R, std::abs(tmpR));
//         //     }
//         //     if (R > 0.01 && h / 2 >= min_h) {
//         //         h /= 2;
//         //     } else if (R < 0.001 && h * 2 <= max_h) {
//         //         h *= 2;
//         //     }
//         // }

//         //изменение шага
//         if (butcher.size().m != butcher.size().n) {
//             float128_t delta_control = 0;
//             for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//                 delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
//             }
//             control.push_back(control.back() + delta_control);
//             float128_t R = std::abs(Yi[1].back() - control.back()); //norma(Yi[1], control);
//             float128_t eps = 0.01;
//             if (R > eps && h > min_h) {
//                 h /= 2;
//             } else if (R < eps / 64) {
//                 h *= 2;
//             }
//         }
//         // else if (Yi[0].size() > 2) {
//         //     std::vector<float128_t> tmpY(Yi.size());
//         //     float128_t R = 0, tmpR;
//         //     //float128_t h2Y;
//         //     K = ExplicitStep(K, Yi, stepNum - 2, func, butcher, h*2);
//         //     tmpY[0] = Yi[0][stepNum - 2] + 2 * h;
//         //     for (uint64_t j = 0; j < orderOfTask; ++j) {
//         //         float128_t delta = 0;
//         //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
//         //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
//         //         }
//         //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
//         //         tmpY[j + 1] = Yi[j + 1][stepNum - 2] + delta;
//         //     }
//         //     for (uint64_t j = 0; j < Yi.size(); ++j) {
//         //         tmpR = (Yi[j][stepNum - 1] - tmpY[j]) / (std::pow(2, orderOfApprox) - 1);
//         //         Yi[j][stepNum - 1] += tmpR;
//         //         R = std::max(R, std::abs(tmpR));
//         //     }
//         //     if (R > 0.01 && h / 2 > min_h) {
//         //         h /= 2;
//         //     } else if (R < 0.001) {
//         //         h *= 2;
//         //     }
//         // }
//         //корректировка шага
//         if (Yi[0].back() + h > Xn) {
//             if (Xn - Yi[0].back() > min_h) {
//                 h = Xn - Yi[0].back();
//             } else {
//                 break;
//             }
//         }
//         ++stepNum;
//     }
//     std::cout << "Tough coeff = " << tough << "\n";
//     //return std::make_pair(Yi[0], Yi[1]);
//     return Yi;
// }

std::vector<std::vector<float128_t>> KoshiSolver (SolveMethod method, const KoshiTask &task, const Matrix<float128_t> butcher, float128_t h, IterationAlgo iter_alg, float128_t approx) {
    // switch (method) {
    //     case SolveMethod::RUNGE_KUTTA:
    //         return RungeKutta(task, butcher, h);
    //     case SolveMethod::FALBERG:
    //     case SolveMethod::CHESKINO:
    //     case SolveMethod::MERSON:
    //         return Falberg(task, butcher, h);
    //     default:
    //         return NonExpl(task, butcher, h, iter_alg, approx);
    // }
    float128_t X0, Xn;
    std::tie(X0, Xn) = task.getBorders();
    const auto &func = task.getODE();
    const auto &Y0 = task.getY0();
    uint64_t orderOfTask = func.size(), orderOfApprox = butcher.size().m - 1;
    float128_t min_h = (Xn - X0) / MIN_H_SIZE;
    float128_t max_h = (Xn - X0) / MAX_H_SIZE;
    float128_t tough = 0;

    std::vector<std::vector<float128_t>> Yi(orderOfTask + 1); //X Y Y' Y'' ...
    std::vector<float128_t> control;

    //инициализация начальных условий
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(Y0[i]);
    }
    control.push_back(Y0[0]);

    //order of approx
    //----->
    //K1 K2 K3 K4 | order of task
    //L1 L2 L3 L4 |
    //M1 M2 M3 M4 V
    std::vector<std::vector<float128_t>> K(orderOfTask, std::vector<float128_t>(orderOfApprox, 0));

    uint64_t stepNum = 1;
    bool expl;
    switch (method) {
        case SolveMethod::RUNGE_KUTTA:
        case SolveMethod::FALBERG:
        case SolveMethod::CHESKINO:
        case SolveMethod::MERSON:
            expl = true;
            break;
        default:
            expl = false;
            break;
    }
    if (expl) {
        std::cout << "explicit\n";
    } else {
        std::cout << "not explicit\n";
    }


    //пока не достигли конца интегрирования
    while (true) {
        //std::cout << "kek\n";
        
        //итерация или явный шаг
        // switch (method) {
        //     case SolveMethod::RUNGE_KUTTA:
        //     case SolveMethod::FALBERG:
        //     case SolveMethod::CHESKINO:
        //     case SolveMethod::MERSON:
        //         K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
        //         break;
        //     default:
        //         if (stepNum == 1) {
        //             K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
        //         }
        //         K = IterationStep(K, Yi, stepNum - 1, func, butcher, h, approx, iter_alg);
        //         break;
        // }
        if (!expl) {
            if (stepNum == 1) {
                K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
            }
            K = IterationStep(K, Yi, stepNum - 1, func, butcher, h, approx, iter_alg);
        } else {
            K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
        }
        //tough = std::max(tough, getToughK(K));

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            float128_t delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            if (!expl) {
                delta *= h;
            }
            //std::cout << "delta: " << delta << "\n";
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }
        // std::vector<float128_t> Y(Yi.size());
        // for (uint64_t j = 0; j < Yi.size(); ++j) {
        //     Y[j] = Yi[j].back();
        // }
        //tough = std::max(tough, ToughCoeff(func, Y));
        //std::cout << "Tough " << stepNum << ": " << tough << "\n";

        //изменение шага
        if (butcher.size().m != butcher.size().n) {
            float128_t delta_control = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
            }
            control.push_back(control.back() + delta_control);
            float128_t R = std::abs(Yi[1].back() - control.back()); //norma(Yi[1], control);
            if (R > approx && h > min_h) {
                h /= 2;
            } else if (R < approx / 64) {
                h *= 2;
            }
        }
        //  else if (Yi[0].size() > 1) {
        //     std::vector<float128_t> tmpY2h(Yi.size());
        //     for (uint64_t i = 0; i < Yi.size(); ++i) {
        //         tmpY2h[i] = Yi[i][stepNum];
        //     }
        //     float128_t R = 0, tmpR;
        //     //float128_t h2Y;
        //     float128_t half_step = h / 2;
        //     for (uint64_t i = 0; i < 2; ++i) {}

        //     K = IterationStep(K, Yi, stepNum - 1, func, butcher, half_step, approx, iter_alg);
        //     Yi[0][stepNum] = Yi[0][stepNum - 1] + half_step;
        //     for (uint64_t j = 0; j < orderOfTask; ++j) {
        //         float128_t delta = 0;
        //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
        //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
        //         }
        //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        //         Yi[j + 1][stepNum] = Yi[j + 1][stepNum - 1] + delta;
        //     }

        //     //K = ExplicitStep(K, Yi, stepNum, func, butcher, half_step);
        //     K = IterationStep(K, Yi, stepNum, func, butcher, half_step, approx, iter_alg);
        //     Yi[0][stepNum] = Yi[0][stepNum] + half_step;
        //     for (uint64_t j = 0; j < orderOfTask; ++j) {
        //         float128_t delta = 0;
        //         for (uint64_t k = 1; k <= orderOfApprox; ++k) {
        //             delta += butcher(orderOfApprox, k) * K[j][k - 1];
        //         }
        //         //Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        //         Yi[j + 1][stepNum] = Yi[j + 1][stepNum] + delta;
        //     }
        //     // for (uint64_t i = 0; i < Yi.size(); ++i) {
        //     //     Yi[i][stepNum] = tmpY1h[i];
        //     // }
        //     for (uint64_t j = 0; j < Yi.size(); ++j) {
        //         tmpR = (Yi[j][stepNum] - tmpY2h[j]) / (std::pow(2, orderOfApprox) - 1);
        //         Yi[j][stepNum] += tmpR;
        //         R = std::max(R, std::abs(tmpR));
        //     }
        //     if (R > approx && h / 2 >= min_h) {
        //         h /= 2;
        //     } else if (R < approx / 10 && h * 2 <= max_h) {
        //         h *= 2;
        //     }
        // }

        //корректировка шага
        if (Yi[0].back() + h > Xn) {
            if (Xn - Yi[0].back() > min_h) {
                h = Xn - Yi[0].back();
            } else {
                break;
            }
        }
        ++stepNum;
    }
    std::cout << "Tough coeff = " << tough << "\n";
    return Yi;
}