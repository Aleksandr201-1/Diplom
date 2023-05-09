#include "ChemicalSolver.hpp"

float128_t NewtonFind (const std::function<float128_t (float128_t)> &f, float128_t x, float128_t approx) {
    uint64_t count = 1;
    float128_t x2 = x - f(x) / derivative(f, x, 0.001, DiffConfig::POINTS4_ORDER1_WAY2);
    while (std::abs(x2 - x) > approx) {
        x = x2;
        x2 = x - f(x) / derivative(f, x, 0.001, DiffConfig::POINTS4_ORDER1_WAY2);
        ++count;
        if (count > ITERATION_CAP) {
            throw std::runtime_error("NewtonFind: the maximum number of iterations has been reached");
        }
    }
    return x2;
}

float128_t updateRho (const std::vector<std::vector<float128_t>> &Yi, float128_t P, float128_t R, float128_t T, uint64_t i) {
    float128_t rho = 0;
    for (uint64_t j = 1; j < Yi.size() - 2; ++j) {
        rho += Yi[j][i];
    }
    rho = P / (R * T * rho);
    return rho;
}

std::vector<std::function<float128_t (float128_t)>> generateUi (const std::vector<std::function<float128_t (float128_t)>> &GFunc,  float128_t R) {
    std::vector<std::function<float128_t (float128_t)>> Ui;
    for (uint64_t j = 0; j < GFunc.size(); ++j) {
        auto u = [=] (float128_t T) -> float128_t {
            return GFunc[j](T) - T * derivative(GFunc[j], T, 0.001, DiffConfig::POINTS4_ORDER1_WAY1) - R * T;
        };
        Ui.push_back(u);
    }
    return Ui;
}

float128_t getU (const std::vector<std::function<float128_t (float128_t)>> &Ui, const std::vector<float128_t> &Y, float128_t T) {
    float128_t U = 0;
    for (uint64_t i = 0; i < Ui.size(); ++i) {
        U += Y[i] * Ui[i](T);
    }
    return U;
}

float128_t updateT (const std::vector<std::vector<float128_t>> &Yi, const std::vector<std::function<float128_t (float128_t)>> &Ui, float128_t T, float128_t U, uint64_t i) {
    std::vector<float128_t> Y(Yi.size() - 3);
    for (uint64_t j = 0; j < Y.size(); ++j) {
        Y[j] = Yi[j + 1][i];
    }
    auto f = [=] (float128_t t) -> float128_t {
        return std::abs(getU(Ui, Y, t) - U);
    };
    return NewtonFind(f, T, 0.001);
}

std::vector<std::vector<float128_t>> ChemicalSolver (SolveMethod method,
                                                     const ChemicalSystem &task,
                                                     const Matrix<float128_t> butcher,
                                                     float128_t h,
                                                     IterationAlgo iter_alg,
                                                     float128_t approx,
                                                     bool Tconst) {
    float128_t X0 = 0, rho = 0, T = task.getTemperature(), P = task.getPressure(), R = ChemicalSystem::R, U;
    auto func = task.getODE();
    auto GFunc = task.getGFunc();
    auto Y = task.getY0();
    auto Ui = generateUi(GFunc, R);
    uint64_t orderOfTask = func.size(), orderOfApprox = butcher.size().m - 1;
    //float128_t min_h = (Xn - X0) / MIN_H_SIZE;
    //float128_t max_h = (Xn - X0) / MAX_H_SIZE;
    float128_t tough = 0;

    std::vector<std::vector<float128_t>> Yi(orderOfTask + 3); //t W1 W2 W3 ... Wn rho T (n = orderOfTask)
    std::vector<float128_t> control;

    //std::cout << Yi.size() << " " << Yi[0].size() << " " << orderOfTask << " " << orderOfApprox << "\n";
    //exit(0);

    //инициализация начальных условий
    Yi[0].push_back(X0);
    for (uint64_t i = 0; i < orderOfTask; ++i) {
        Yi[i + 1].push_back(Y[i]);
    }
    Yi[orderOfTask + 1].push_back(updateRho(Yi, P, R, T, 0));
    Yi[orderOfTask + 2].push_back(T);
    control.push_back(Y[0]);

    U = getU(Ui, Y, T);
    std::cout << "U: " << U << "\n";

    // for (uint64_t j = 1; j < task.Y.size() - 2; ++j) {
    //     rho += task.Y[j];
    // }
    // rho = P / (R * T * rho);

    //order of approx
    //----->
    //K1 K2 K3 K4 | order of task
    //L1 L2 L3 L4 |
    //M1 M2 M3 M4 V
    std::vector<std::vector<float128_t>> K(orderOfTask, std::vector<float128_t>(orderOfApprox, 0));
    //std::vector<float128_t> args(orderOfTask + 2, 0);

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
    //T = 3285
    auto Ytough = Y;
    Ytough.insert(Ytough.begin(), X0);
    Ytough.insert(Ytough.end(), updateRho(Yi, P, R, T, 0));
    Ytough.insert(Ytough.end(), T);
    printVector(Ytough);
    tough = ToughCoeff(func, Ytough);
    //exit(0);

    //пока не достигли конца интегрирования
    while (true) {
        //инициализация аргументов
        //for (uint64_t i = 0; i < args.size() - 1; ++i) {
        //    args[i] = Yi[i].back();
        //}
        
        //итерация или явный шаг
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
        if (Tconst) {
            Yi[orderOfTask + 1].push_back(updateRho(Yi, P, R, T, stepNum));
            Yi[orderOfTask + 2].push_back(Yi[orderOfTask + 2].back());
        } else {
            float128_t tmp = Yi[orderOfTask + 2].back();
            Yi[orderOfTask + 1].push_back(Yi[orderOfTask + 1].back());
            Yi[orderOfTask + 2].push_back(updateT(Yi, Ui, tmp, U, stepNum));
        }
        //exit(0);
        // std::vector<float128_t> Y(Yi.size());
        // for (uint64_t j = 0; j < Yi.size(); ++j) {
        //     Y[j] = Yi[j].back();
        // }
        //tough = std::max(tough, ToughCoeff(func, Y));
        //std::cout << "Tough " << stepNum << ": " << tough << "\n";
        //exit(0);

        //изменение шага
        if (butcher.size().m != butcher.size().n) {
            float128_t delta_control = 0;
            for (uint64_t k = 1; k <= orderOfApprox; ++k) {
                delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
            }
            control.push_back(control.back() + delta_control);
            float128_t R = std::abs(Yi[1].back() - control.back()); //norma(Yi[1], control);
            if (R > approx) {
                h /= 2;
            } else if (R < approx / 64) {
                h *= 2;
            }
        }
        //окончание вычисления
        bool stop = true;
        for (uint64_t i = 1; i < Yi.size() - 2; ++i) {
            if (std::abs(Yi[i][stepNum] - Yi[i][stepNum - 1]) > 1e-6) {
                stop = false;
                break;
            }
        }
        if (stop || stepNum > 200) {
            break;
        }
        ++stepNum;
    }
    for (uint64_t i = 0; i < Y.size(); ++i) {
        Y[i] = Yi[i + 1].back();
    }
    T = Yi.back().back();
    U = getU(Ui, Y, T);
    std::cout << "U end: " << U << "\n";
    std::cout << "Tough coeff = " << tough << "\n";
    return Yi;
}