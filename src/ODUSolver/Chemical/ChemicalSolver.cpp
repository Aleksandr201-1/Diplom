#include "ChemicalSolver.hpp"

std::vector<std::vector<float128_t>> drop (const std::vector<std::vector<float128_t>> &vv) {
    uint64_t dropSize = 100;
    std::vector<std::vector<float128_t>> ans(vv.size());
    for (uint64_t i = 0; i < vv[0].size(); ++i) {
        if (i % dropSize == 0) {
            for (uint64_t j = 0; j < vv.size(); ++j) {
                ans[j].push_back(vv[j][i]);
            }
        }
    }
    for (uint64_t j = 0; j < vv.size(); ++j) {
        ans[j].push_back(vv[j].back());
    }
    return ans;
}

float128_t NewtonFind (const std::function<float128_t (float128_t)> &f, float128_t x, float128_t approx) {
    uint64_t count = 1;
    float128_t x2 = x - f(x) / derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1);
    while (std::abs(x2 - x) > approx) {
        x = x2;
        x2 = x - f(x) / derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1);
        ++count;
        //if (count > ITERATION_CAP) {
        //    throw std::runtime_error("NewtonFind: the maximum number of iterations has been reached");
        //}
    }
    return x2;
}

// float128_t SIFind (const std::function<float128_t (float128_t)> &f, float128_t x, float128_t approx) {
//     uint64_t count = 1;
//     float128_t x2 = x - f(x) / derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1);
//     while (std::abs(x2 - x) > approx) {
//         x = x2;
//         x2 = x - f(x) / derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1);
//         ++count;
//         //if (count > ITERATION_CAP) {
//         //    throw std::runtime_error("NewtonFind: the maximum number of iterations has been reached");
//         //}
//     }
//     return x2;
// }

bool checkSI (float128_t a, float128_t b, const std::function<float128_t (float128_t)> &f, float128_t q) {
    float128_t e = 0.000000001, x = a;
    while (x + e < b) {
        x += e;
        if (derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1) > q) {
            return false;
        }
        if (f(x) < a || f(x) > b) {
            return false;
        }
    }
    return true;
}

float128_t SIFind (const std::function<float128_t (float128_t)> &f, float128_t x, float128_t a, float128_t b, float128_t approx) {
    float128_t q = std::max(std::abs(derivative(f, a, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1)), std::abs(derivative(f, b, 0.001, DiffConfig::POINTS2_ORDER1_WAY1)));
    if (!checkSI(a, b, f, q)) {
        std::cout << "fi(x) = " << derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1) << "\n";
        throw std::runtime_error("SI: Incorrect function fi(x). fi'(x) > q, x in (a, b) or fi(x) not in (a, b)");
    }
    uint64_t count = 1;
    float128_t x2 = f(x);
    while ((q / (1 - q)) * (std::abs(x2 - x)) > approx) {
        x = x2;
        x2 = f(x2);
        ++count;
        if (count > ITERATION_CAP) {
            throw std::runtime_error("SI: the maximum number of iterations has been reached");
        }
    }
    return f(x2);
}

float128_t updateRho (const std::vector<std::vector<float128_t>> &Yi, float128_t P, float128_t R, float128_t T, uint64_t i) {
    float128_t rho = 0;
    for (uint64_t j = 1; j < Yi.size() - 2; ++j) {
        rho += Yi[j][i];
    }
    rho = P / (R * T * rho);
    return rho;
}

std::vector<std::function<float128_t (float128_t)>> generateUi (const std::vector<std::function<float128_t (float128_t)>> &GFunc, const std::vector<PhiFunction> &PhiFunc, float128_t R) {
    std::vector<std::function<float128_t (float128_t)>> Ui;
    for (uint64_t j = 0; j < GFunc.size(); ++j) {
        auto u = [=] (float128_t T) -> float128_t {
            //return GFunc[j](T) - T * derivative(GFunc[j], T, 0.1, DiffConfig::POINTS3_ORDER1_WAY1) - R * T;
            return GFunc[j](T) + T * (PhiFunc[j](T) + T*PhiFunc[j].der1(T)) - R * T;
        };
        Ui.push_back(u);
    }
    return Ui;
}

float128_t getU (const std::vector<std::function<float128_t (float128_t)>> &Ui, const std::vector<float128_t> &Y, float128_t T) {
    float128_t U = 0;
    for (uint64_t i = 0; i < Ui.size(); ++i) {
        U += Y[i + 1] * Ui[i](T);
    }
    return U;
}

float128_t updateT (const std::vector<std::vector<float128_t>> &Yi, const std::vector<std::function<float128_t (float128_t)>> &Ui, float128_t U, uint64_t i) {
    std::vector<float128_t> Y(Yi.size());
    for (uint64_t j = 0; j < Y.size(); ++j) {
        Y[j] = Yi[j][i];
    }
    Y[0] = Yi[0][i - 1];
    Y[Y.size() - 2] = Yi[Y.size() - 2][i - 1];
    Y[Y.size() - 1] = Yi[Y.size() - 1][i - 1];
    float128_t R = 8.3144;
    auto f = [=] (float128_t t) -> float128_t {
        return std::abs(getU(Ui, Y, t) - U);
        //return std::abs(U - getU(Ui, Y, t)) * (0.000005) + t;
    };
    //std::cout << Y.back() << " " << Yi.back().back() << "\n";
    return NewtonFind(f, Y.back(), 0.0001);
    //return SIFind(f, Y.back(), Y.back() - 1, Y.back() + 1, 0.0001);
}

std::vector<std::vector<float128_t>> ChemicalSolver (SolveMethod method,
                                                     const ChemicalSystem &task,
                                                     const Matrix<float128_t> butcher,
                                                     float128_t h,
                                                     IterationAlgo iter_alg,
                                                     float128_t approx,
                                                     ReactionType type) {
    float128_t X0 = 0, rho = 0, T = task.getTemperature(), P = task.getPressure(), R = ChemicalSystem::R, U;
    auto func = task.getODE();
    auto GFunc = task.getGFunc();
    auto PhiFunc = task.getPhiFunc();
    // for (uint64_t i = 0; i < PhiFunc.size(); ++i) {
    //     std::cout << i << "\n";
    //     std::cout << "Phi(" << T << ") = " << PhiFunc[i](T) << "\n";
    //     std::cout << "Phi'(" << T << ") = " << PhiFunc[i].der1(T) << "\n";
    //     std::cout << "Phi''(" << T << ") = " << PhiFunc[i].der2(T) << "\n\n";

    //     std::cout << "Phi(" << T << ") = " << PhiFunc[i](T) << "\n";
    //     std::cout << "Numeric Phi'(" << T << ") = " << derivative( PhiFunc[i], T, 0.00001, DiffConfig::POINTS5_ORDER1_WAY1) << "\n";
    //     std::cout << "Numeric Phi''(" << T << ") = " << derivative( PhiFunc[i], T, 0.00001, DiffConfig::POINTS5_ORDER2_WAY1) << "\n\n";
    // }
    //exit(0);
    auto Y = task.getY0();
    //Y = {0.688402337883, 10.28972285, 22.0855861377, 8.537531053, 2.9605809371, 9.34790964455};
    //T = 3400;
    printVector(Y);
    auto Ui = generateUi(GFunc, PhiFunc, R);
    float128_t tough = 0;
    uint64_t orderOfTask = func.size()+2, orderOfApprox = butcher.size().m - 1;

    std::vector<std::vector<float128_t>> Yi(orderOfTask + 1); //t W1 W2 W3 ... Wn-2 rho T (n = orderOfTask)
    std::vector<float128_t> control;

    //инициализация начальных условий
    Y.insert(Y.begin(), 0);
    for (uint64_t i = 0; i < Y.size(); ++i) {
        Yi[i].push_back(Y[i]);
    }
    Y.push_back(updateRho(Yi, P, R, T, 0));
    Y.push_back(T);
    Yi[Yi.size() - 2].push_back(updateRho(Yi, P, R, T, 0));
    Yi[Yi.size() - 1].push_back(T);
    control.push_back(Y[0]);

    //Y.push_back(Yi[Yi.size() - 2][0]);
    //Y.push_back(Yi[Yi.size() - 1][0]);
    U = getU(Ui, Y, T);

    std::cout << "U: " << U << "\n";
    uint64_t stepNum = 1;

    std::vector<std::function<float128_t (const std::vector<float128_t> &)>> nu;
    for (uint64_t i = 0; i < Ui.size(); ++i) {
        auto nui = [=] (const std::vector<float128_t> &X) -> float128_t {
            float128_t Temp = X[X.size() - 1], Rho = X[X.size() - 2];
            return R * Temp * std::log((Rho * R * Temp * X[i + 1]) / P) + GFunc[i](Temp);
        };
        nu.push_back(nui);
    }

    //выбор типа реакции
    switch(type) {
        case ReactionType::ISOTERM_CONST_RHO: //изотермическая реакция с постоянной плотностью
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    return 0;
                }
            );
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    return 0;
                }
            );
            break;
        case ReactionType::ADIABAT_CONST_RHO: //адиабатическая реакция с постоянной плотностью (не работает)
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    return 0;
                }
            );
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    static uint64_t idx = 0;
                    float128_t ans = 0, cV = 0, cp = 0, currTemp = X[X.size() - 1];
                    //std::cout << X[X.size() - 1] << "fuck\n";
                    //exit(0);
                    for (uint64_t i = 0; i < Ui.size(); ++i) {
                        ans += func[i](X) * Ui[i](currTemp);
                        //return GFunc[j](T) + T * (PhiFunc[j](T) + T*PhiFunc[j].der1(T)) - R * T;
                        //cV += X[i + 1] * derivative(Ui[i], X[X.size() - 1], 0.00001, DiffConfig::POINTS5_ORDER1_WAY1);
                        cV += X[i + 1] * (currTemp * (2*PhiFunc[i].der1(currTemp) + currTemp*PhiFunc[i].der2(currTemp)) - R);
                        cp += R * X[i + 1];
                    }
                    if (idx % 100 == 0) {
                        std::cout << "cp/cv = " << (cp + cV) / cV << '\n';
                    }
                    ans /= cV;
                    ++idx;
                    return -ans;
                    // auto f = [=] (const std::vector<float128_t> &X) -> float128_t {
                    //     return std::abs(getU(Ui, X, X.back()) - U);
                    // };
                    // return derivative(f, X, 0.0001, X.size() - 1, DiffConfig::POINTS4_ORDER1_WAY2);
                    //return NewtonFind(f, X.back(), 0.0001);
                }
            );
            break;
        case ReactionType::ISOTERM_CONST_P: //изотермическая реакция с постоянным давлением
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    float128_t ans = 0, div = 0;
                    for (uint64_t i = 0; i < orderOfTask - 2; ++i) {
                        ans += func[i](X);
                        div += X[i + 1];
                    }
                    //return 0;
                    return (-ans / div) * X[X.size() - 2];
                }
            );
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    return 0;
                }
            );
            break;
        case ReactionType::ADIABAT_CONST_P: //адиабатическая реакция с постоянным давлением
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    return 0;
                }
            );
            func.push_back(
                [=] (const std::vector<float128_t> &X) -> float128_t {
                    float128_t ans = 0, cV = 0;
                    for (uint64_t i = 0; i < orderOfTask - 2; ++i) {
                        ans += func[i](X) * Ui[i](X[X.size() - 1]);
                        cV += X[i + 1] * derivative(Ui[i], X[X.size() - 1], 0.0001, DiffConfig::POINTS4_ORDER1_WAY2);
                    }
                    ans /= X[X.size() - 2] * cV;
                    return ans;
                    // auto f = [=] (const std::vector<float128_t> &X) -> float128_t {
                    //     return std::abs(getU(Ui, X, X.back()) - U);
                    // };
                    // return derivative(f, X, 0.0001, X.size() - 1, DiffConfig::POINTS4_ORDER1_WAY2);
                    // return NewtonFind(f, X.back(), 0.0001);
                }
            );
            break;
    }

    //std::cout << "U: " << U << "\n";

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
    auto Ytough = Y;
    Ytough.insert(Ytough.begin(), X0);
    Ytough.insert(Ytough.end(), updateRho(Yi, P, R, T, 0));
    Ytough.insert(Ytough.end(), T);

    //пока не достигли конца интегрирования
    while (true) {
        //инициализация аргументов
        //for (uint64_t i = 0; i < args.size() - 1; ++i) {
        //    args[i] = Yi[i].back();
        //}

        for (uint64_t i = 0; i < Yi.size(); ++i) {
            Yi[i].push_back(0);
        }
        //Yi[Yi.size() - 2][stepNum] = Yi[Yi.size() - 2][stepNum - 1];
        //Yi[Yi.size() - 1][stepNum] = Yi[Yi.size() - 1][stepNum - 1];

        //float128_t oldT = 0, newT = Yi[Yi.size() - 1][stepNum - 1];

        //if (stepNum % 10 == 0) {
            //std::cout << stepNum << "\n";
        //}

        //while (std::abs(oldT - newT) > 0.0001) {

            //oldT = newT;
        
            //итерация или явный шаг
            if (!expl) {
                if (stepNum == 1) {
                    K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
                }
                K = IterationStep(K, Yi, stepNum - 1, func, butcher, h, approx, iter_alg);
            } else {
                K = ExplicitStep(K, Yi, stepNum - 1, func, butcher, h);
            }
            //K = IterationStep(K, Yi, stepNum - 1, func, butcher, h, approx, iter_alg);
            //tough = std::max(tough, getToughK(K));

            //добавление нового значения
            Yi[0][stepNum] = Yi[0][stepNum - 1] + h;
            for (uint64_t j = 0; j < K.size(); ++j) { //orderOfTask
                float128_t delta = 0;
                for (uint64_t k = 1; k <= K[0].size(); ++k) { //orderOfApprox
                    delta += butcher(orderOfApprox, k) * K[j][k - 1];
                }
                if (!expl) {
                    delta *= h;
                }
                //std::cout << "delta: " << delta << "\n";
                Yi[j + 1][stepNum] = Yi[j + 1][stepNum - 1] + delta;
            }
            // if (type != ReactionType::ADIABAT_CONST_RHO) {
            //     Yi[Yi.size() - 2][stepNum] = updateRho(Yi, P, R, T, stepNum);
            //     Yi[Yi.size() - 1][stepNum] = Yi[Yi.size() - 1][stepNum - 1];
            // } else {
            //     //float128_t tmp = Yi[orderOfTask + 2].back();
            //     Yi[Yi.size() - 2][stepNum] = Yi[Yi.size() - 2][stepNum - 1];
            //     Yi[Yi.size() - 1][stepNum] = updateT(Yi, Ui, U, stepNum);
            // }
            //newT = Yi[Yi.size() - 1][stepNum];
            // if (Tconst) {
            //     Yi[orderOfTask + 1].push_back(updateRho(Yi, P, R, T, stepNum));
            //     Yi[orderOfTask + 2].push_back(Yi[orderOfTask + 2].back());
            // } else {
            //     //float128_t tmp = Yi[orderOfTask + 2].back();
            //     Yi[orderOfTask + 1].push_back(Yi[orderOfTask + 1].back());
            //     Yi[orderOfTask + 2].push_back(updateT(Yi, Ui, U, stepNum));
            // }
        //}

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
        // in M: 4.85803
        // out M: 4.85803
        //std::cout << "\n";
        // if (stepNum == 3) {
        //     exit(0);
        // }
        //exit(0);
        //окончание вычисления
        bool stop = true;
        for (uint64_t i = 1; i < Yi.size() - 2; ++i) {
            if (std::abs(Yi[i][stepNum] - Yi[i][stepNum - 1]) > 1e-6) {
                stop = false;
                break;
            }
        }
        if (stepNum > 200) {
            break;
        }
        ++stepNum;
    }
    for (uint64_t i = 0; i < Y.size(); ++i) {
        Y[i] = Yi[i].back();
    }
    std::cout << "\n";
    for (uint64_t i = 0; i < task.getCount(); ++i) {
        float128_t qin = 0, qout = 0;
        auto gammaIn = task[i].getInput();
        auto gammaOut = task[i].getOutput();
        for (uint64_t j = 0; j < func.size() - 2; ++j) {
            qin += nu[j](Y) * gammaIn[j];
            qout += nu[j](Y) * gammaOut[j];
            //qin += func[j](Y) * gammaIn[j];
            //qout += func[j](Y) * gammaOut[j];
        }
        std::cout << "In: " << qin << "\nOut: " << qout << "\n";
    }
    T = Yi.back().back();
    U = getU(Ui, Y, T);
    auto subs = task.getSubstanceList();
    auto subs_mass = task.getSubstanceMasses();
    std::cout << "ending concentrations:\n";
    float128_t rho_gamma = 0;
    for (uint64_t j = 0; j < Yi.size() - 3; ++j) {
        rho_gamma += Yi[j + 1].back() * subs_mass.find(subs[j])->second;
    }
    for (uint64_t i = 0; i < Yi.size() - 3; ++i) {
        float128_t curr_conc = Yi[i + 1].back();
        curr_conc *= rho_gamma;
        std::cout << curr_conc << "\n";
    }
    std::cout << "U end: " << U << "\n";
    //std::cout << "Tough coeff = " << tough << "\n";
    return Yi;
}