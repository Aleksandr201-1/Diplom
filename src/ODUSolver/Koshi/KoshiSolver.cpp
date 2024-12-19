#include "KoshiSolver.hpp"

std::vector<std::vector<double>> KoshiSolver (SolveMethod method, const KoshiTask &task, const Matrix<double> butcher, double h, IterationAlgo iter_alg, double approx) {
    double X0, Xn;
    std::tie(X0, Xn) = task.getBorders();
    const auto &func = task.getODE();
    const auto &Y0 = task.getY0();
    uint64_t orderOfTask = func.size(), orderOfApprox = butcher.size().m - 1;
    double min_h = (Xn - X0) / MIN_H_SIZE;
    double max_h = (Xn - X0) / MAX_H_SIZE;
    double tough = 0;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //X Y Y' Y'' ...
    std::vector<double> control;

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
    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));

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
    // if (expl) {
    //     std::cout << "explicit\n";
    // } else {
    //     std::cout << "not explicit\n";
    // }


    //пока не достигли конца интегрирования
    while (true) {
        //std::cout << "kek\n";
        if (!expl) {
            if (stepNum == 1) {
                K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
            }
            K = IterationStep(K, Yi, stepNum - 1, func, butcher, h, approx, iter_alg);
        } else {
            K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
        }
        //std::cout << stepNum << "\n";
        //tough = std::max(tough, getToughK(K));

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            double delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            if (!expl) {
                delta *= h;
            }
            //std::cout << "delta: " << delta << "\n";
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

        //изменение шага
        if (butcher.size().m != butcher.size().n) {
            double delta_control = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
            }
            control.push_back(control.back() + delta_control);
            double R = std::abs(Yi[1].back() - control.back()); //norma(Yi[1], control);
            if (R > approx && h > min_h) {
                h /= 2;
            } else if (R < approx / 64) {
                h *= 2;
            }
        }
        //  else if (Yi[0].size() > 1) {
        //     std::vector<double> tmpY2h(Yi.size());
        //     for (uint64_t i = 0; i < Yi.size(); ++i) {
        //         tmpY2h[i] = Yi[i][stepNum];
        //     }
        //     double R = 0, tmpR;
        //     //double h2Y;
        //     double half_step = h / 2;
        //     for (uint64_t i = 0; i < 2; ++i) {}

        //     K = IterationStep(K, Yi, stepNum - 1, func, butcher, half_step, approx, iter_alg);
        //     Yi[0][stepNum] = Yi[0][stepNum - 1] + half_step;
        //     for (uint64_t j = 0; j < orderOfTask; ++j) {
        //         double delta = 0;
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
        //         double delta = 0;
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