#include "KoshiSolver.hpp"

//модуль определения жёсткости
//сохранение функционального дерева в файл
//добавить в метод рунге-кутты возможность вставлять критерий выбора шага
//итерационные методы (неявные, то есть с полной таблицей бутчера)
//печать отчёта
//графический интерфейс

float128_t getToughK (const std::vector<std::vector<float128_t>> &K) {
    float128_t max = 0;
    
    for (uint64_t i = 0; i < K.size(); ++i) {
        float128_t tmp = 0;
        for (uint64_t j = 0; j < K[i].size(); ++j) {
            tmp += std::abs(K[i][j]);
        }
        max = std::max(max, tmp);
    }

    for (uint64_t i = 0; i < K[0].size(); ++i) {
        float128_t tmp = 0;
        for (uint64_t j = 0; j < K.size(); ++j) {
            tmp += std::abs(K[j][i]);
        }
        max = std::max(max, tmp);
    }

    return max;
}

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
    uint64_t stepNum = 1;
    while (true) {
        K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
        tough = std::max(tough, getToughK(K));

        //добавление нового значения
        Yi[0].push_back(Yi[0].back() + h);
        for (uint64_t j = 0; j < orderOfTask; ++j) {
            float128_t delta = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

        //изменение шага
        if (Yi[0].back() + h > Xn) {
            if (Xn - Yi[0].back() >= min_h) {
                h = Xn - Yi[0].back();
            } else {
                break;
            }
        }
        ++stepNum;
    }
    return Yi;
}

std::vector<std::vector<float128_t>> KoshiSolver (SolveMethod method, const KoshiTask &task, const Matrix<float128_t> butcher, float128_t h, IterationAlgo iter_alg, float128_t approx) {
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
        if (0) {
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
                delta += butcher(orderOfApprox, k) * K[j][k - 1] * h;
            }
            Yi[j + 1].push_back(Yi[j + 1].back() + delta);
        }

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