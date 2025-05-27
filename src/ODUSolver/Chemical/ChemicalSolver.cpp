#include "ChemicalSolver.hpp"

std::vector<std::vector<double>> drop (const std::vector<std::vector<double>> &vv, uint64_t newSize) {
    uint64_t dropSize = vv[0].size() / newSize;
    std::vector<std::vector<double>> ans(vv.size());
    for (uint64_t i = 0; i < vv[0].size(); ++i) {
        if ((i + 1) % dropSize == 0) {
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

double NewtonFindT (const std::function<double (double)> &f, double x, double approx, DiffConfig conf) {
    uint64_t count = 1;
    double x2 = x - f(x) / derivative(f, x, approx, conf);
    while (std::abs(x2 - x) > approx) {
        x = x2;
        x2 = x - f(x) / derivative(f, x, approx, conf);
        ++count;
        if (count > 3000) {
            std::cout << "NewtonFindT: the maximum number of iterations has been reached\n";
            std::cout.flush();
            exit(1);
            return x2;
            //throw std::runtime_error("NewtonFindT: the maximum number of iterations has been reached");
        }
    }
    return x2;
}

// double SIFind (const std::function<double (double)> &f, double x, double approx) {
//     uint64_t count = 1;
//     double x2 = x - f(x) / derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1);
//     while (std::abs(x2 - x) > approx) {
//         x = x2;
//         x2 = x - f(x) / derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1);
//         ++count;
//         //if (count > ITERATION_CAP) {
//         //    throw std::runtime_error("NewtonFindT: the maximum number of iterations has been reached");
//         //}
//     }
//     return x2;
// }

bool checkSI (double a, double b, const std::function<double (double)> &f, double q) {
    double e = 0.000000001, x = a;
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

double SIFind (const std::function<double (double)> &f, double x, double a, double b, double approx) {
    double q = std::max(std::abs(derivative(f, a, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1)), std::abs(derivative(f, b, 0.001, DiffConfig::POINTS2_ORDER1_WAY1)));
    if (!checkSI(a, b, f, q)) {
        //std::cout << "fi(x) = " << derivative(f, x, 0.0001, DiffConfig::POINTS2_ORDER1_WAY1) << "\n";
        throw std::runtime_error("SI: Incorrect function fi(x). fi'(x) > q, x in (a, b) or fi(x) not in (a, b)");
    }
    uint64_t count = 1;
    double x2 = f(x);
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

double updateRho (const std::vector<std::vector<double>> &Yi, double P, double R, double T, uint64_t i) {
    double rho = 0;
    for (uint64_t j = 1; j < Yi.size() - 2; ++j) {
        rho += Yi[j][i];
    }
    rho = P / (R * T * rho);
    return rho;
}

std::vector<std::function<double (double)>> generateUi (const std::vector<std::function<double (double)>> &GFunc, const std::vector<PhiFunction> &PhiFunc, double R) {
    std::vector<std::function<double (double)>> Ui;
    for (uint64_t j = 0; j < GFunc.size(); ++j) {
        auto u = [=] (double T) -> double {
            //return GFunc[j](T) - T * derivative(GFunc[j], T, 0.01, DiffConfig::POINTS2_ORDER1_WAY3) - R * T;
            return GFunc[j](T) + T * (PhiFunc[j](T) + T * PhiFunc[j].der1(T)) - R * T;
            //return T * (2 * PhiFunc[j].der1(T) + T * PhiFunc[j].der2(T)) - R;
        };
        Ui.push_back(u);
    }
    return Ui;
}

double getU (const std::vector<std::function<double (double)>> &Ui, const std::vector<double> &Y, double T) {
    double U = 0;
    for (uint64_t i = 0; i < Ui.size(); ++i) {
        U += Y[i + 1] * Ui[i](T);
    }
    return U;
}

enum RecationPeriod {
    PRE_REACT = 0,
    ACTIVE_REACT,
    POST_REACT
};

double updateT (const std::vector<std::vector<double>> &Yi, const std::vector<std::function<double (double)>> &Ui, double U, uint64_t i) {
    std::vector<double> Y(Yi.size());
    for (uint64_t j = 0; j < Y.size(); ++j) {
        Y[j] = Yi[j][i];
    }
    Y[0] = Yi[0][i - 1];
    Y[Y.size() - 2] = Yi[Y.size() - 2][i - 1];
    Y[Y.size() - 1] = Yi[Y.size() - 1][i - 1];
    double R = 8.3144;
    auto f = [=] (double t) -> double {
        return std::abs(getU(Ui, Y, t) - U);
    };
    return NewtonFindT(f, Y.back(), 0.0001);
}

std::vector<std::vector<double>> ChemicalSolver (SolveMethod method,
                                                     const ChemicalSystem &task,
                                                     const Matrix<double> butcher,
                                                     double h_min,
                                                     double h_max,
                                                     double h_last,
                                                     IterationAlgo iter_alg,
                                                     double approx,
                                                     ReactionType type) {
    double X0 = 0, rho = 0, T = task.getTemperature(), P = task.getPressure(), R = ChemicalSystem::R, U;
    //std::cout << "getting ODE\n";
    double currLen = 0;
    auto func = task.getODE();
    //std::cout << "getting G-funcs\n";
    auto GFunc = task.getGFunc();
    //std::cout << "getting Phi funcs\n";
    auto PhiFunc = task.getPhiFunc();
    // for (uint64_t i = 0; i < PhiFunc.size(); ++i) {
    //     //std::cout << i << "\n";
    //     //std::cout << "Phi(" << T << ") = " << PhiFunc[i](T) << "\n";
    //     //std::cout << "Phi'(" << T << ") = " << PhiFunc[i].der1(T) << "\n";
    //     //std::cout << "Phi''(" << T << ") = " << PhiFunc[i].der2(T) << "\n\n";

    //     //std::cout << "Phi(" << T << ") = " << PhiFunc[i](T) << "\n";
    //     //std::cout << "Numeric Phi'(" << T << ") = " << derivative( PhiFunc[i], T, 0.00001, DiffConfig::POINTS5_ORDER1_WAY1) << "\n";
    //     //std::cout << "Numeric Phi''(" << T << ") = " << derivative( PhiFunc[i], T, 0.00001, DiffConfig::POINTS5_ORDER2_WAY1) << "\n\n";
    // }
    //exit(0);
    auto Y = task.getY0();
    auto ental = task.getEnthalpy();
    auto Ui = generateUi(GFunc, PhiFunc, R);
    double tough = 0;
    uint64_t orderOfTask = func.size()+2, orderOfApprox = butcher.size().m - 1;
    RecationPeriod period = PRE_REACT;

    std::vector<std::vector<double>> Yi(orderOfTask + 1); //t W1 W2 W3 ... Wn-2 rho T (n = orderOfTask)
    std::vector<double> control, args(orderOfTask + 1);

    std::ofstream output("chem_log.txt"), tough_output("tough_output.txt");
    output << "#\tt\t";
    auto names = task.getSubstanceList();
    for (auto el : names) {
        output << el << "\t";
    }
    output << "rho\tT\n";

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
    std::cout.precision(15);
    std::cout.setf(std::ios_base::fixed);
    U = getU(Ui, Y, T);
    double Ustart = U;

    //std::cout << "U: " << U << "\n";
    //std::cout << "order of task: " << orderOfTask << "\norder of approx: " << orderOfApprox << "\n";
    uint64_t stepNum = 1;

    std::vector<std::function<double (const std::vector<double> &)>> nu;
    for (uint64_t i = 0; i < Ui.size(); ++i) {
        auto nui = [=] (const std::vector<double> &X) -> double {
            double Temp = X[X.size() - 1], Rho = X[X.size() - 2];
            return R * Temp * std::log((Rho * R * Temp * X[i + 1]) / P) + GFunc[i](Temp);
        };
        nu.push_back(nui);
    }

    //выбор типа реакции
    switch(type) {
        case ReactionType::ISOTERM_CONST_RHO: //изотермическая реакция с постоянной плотностью
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    return 0;
                }
            );
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    return 0;
                }
            );
            break;
        case ReactionType::ADIABAT_CONST_RHO: //адиабатическая реакция с постоянной плотностью
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    return 0;
                }
            );
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    double ans = 0, cV = 0, currTemp = X[X.size() - 1], currRho = X[X.size() - 2];
                    for (uint64_t i = 0; i < Ui.size(); ++i) {
                        double ui_der = currTemp * (2 * PhiFunc[i].der1(currTemp) + currTemp * PhiFunc[i].der2(currTemp)) - R;
                        ans += func[i](X) * Ui[i](currTemp);
                        //cV += X[i + 1] * (currTemp * (2*PhiFunc[i].der1(currTemp) + currTemp*PhiFunc[i].der2(currTemp)) - R);
                        cV += X[i + 1] * ui_der;//derivative(Ui[i], currTemp, 0.0001, DiffConfig::POINTS4_ORDER1_WAY2);
                    }
                    //ans /= (cV * currRho); //* rho ???
                    ans /= cV;
                    return -ans;
                }
            );
            break;
        case ReactionType::ISOTERM_CONST_P: //изотермическая реакция с постоянным давлением
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    double ans = 0, div = 0;
                    for (uint64_t i = 0; i < orderOfTask - 2; ++i) {
                        ans += func[i](X);
                        div += X[i + 1]; //???
                    }
                    return (-ans / div) * X[X.size() - 2];
                }
            );
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    return 0;
                }
            );
            break;
        case ReactionType::ADIABAT_CONST_P: //адиабатическая реакция с постоянным давлением
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    double ans = 0, cp = 0, T = X.back(), rho = X[X.size() - 2];
                    for (uint64_t i = 0; i < orderOfTask - 2; ++i) {
                        ans += func[i](X) * ental[i](T);
                        cp += X[i + 1] * derivative(ental[i], T, 0.001, DiffConfig::POINTS4_ORDER1_WAY2);
                    }
                    return -1.0 / (rho * cp) * ans;
                }
            );
            func.push_back(
                [&] (const std::vector<double> &X) -> double {
                    double ans = 0, cP = 0, nuSum = 0, T = X.back(), rho = X[X.size() - 2];;
                    for (uint64_t i = 0; i < orderOfTask - 2; ++i) {
                        nuSum += X[i + 1];
                        cP += X[i + 1] * derivative(ental[i], T, 0.001, DiffConfig::POINTS4_ORDER1_WAY2);
                    }
                    for (uint64_t i = 0; i < orderOfTask - 2; ++i) {
                        ans += (ental[i](T) - cP * T / nuSum) * func[i](X);
                    }
                    ans /= (cP * T);
                    return ans;
                }
            );
            break;
        default:
            exit(2);
    }

    //order of approx
    //----->
    //K1 K2 K3 K4 | order of task
    //L1 L2 L3 L4 |
    //M1 M2 M3 M4 V
    std::vector<std::vector<double>> K(orderOfTask, std::vector<double>(orderOfApprox, 0));

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
    auto Yprev = Ytough;
    double tough_coeff = 0;//ToughCoeff(func, Ytough, Yprev);
    //std::cout << "Tough coeff: " << tough_coeff << "\n";
    //tough_output << tough_coeff << "\n";
    double h = h_min;

    for (uint64_t i = 0; i < Yi.size(); ++i) {
        Yi[i].push_back(Yi[i][0]);
    }
    output << 0 << "\t";
    for (uint64_t j = 0; j < Yi.size() - 1; ++j) {
        output << Yi[j][1] << "\t";
    }
    output << Yi.back()[1] << "\n";

    //пока не достигли конца интегрирования
    while (true) {

        //добавление нового значения
        for (uint64_t j = 0; j < Yi.size(); ++j) { //orderOfTask
            Yi[j][0] = Yi[j][1];
        }

        //итерация или явный шаг
        if (!expl) {
            if (stepNum == 1) {
                K = ExplicitStep(K, Yi, 0, func, butcher, h);
            }
            K = IterationStep(K, Yi, 0, func, butcher, h, approx, iter_alg);
        } else {
            K = ExplicitStep(K, Yi, 0, func, butcher, h);
        }

        //добавление нового значения
        Yi[0][1] = Yi[0][0] + h;
        for (uint64_t j = 0; j < K.size(); ++j) { //orderOfTask
            double delta = 0;
            for (uint64_t k = 1; k <= K[0].size(); ++k) { //orderOfApprox
                delta += butcher(orderOfApprox, k) * K[j][k - 1];
            }
            if (!expl) {
                delta *= h;
            }
            Yi[j + 1][1] = Yi[j + 1][0] + delta;
        }

        //изменение шага
        // if (butcher.size().m != butcher.size().n) {
        //     double delta_control = 0;
        //     for (uint64_t k = 1; k <= orderOfApprox; ++k) {
        //         delta_control += butcher(orderOfApprox + 1, k) * K[0][k - 1];
        //     }
        //     control.push_back(control.back() + delta_control);
        //     double R = std::abs(Yi[1].back() - control.back()); //norma(Yi[1], control);
        //     if (R > approx) {
        //         h /= 2;
        //     } else if (R < approx / 64) {
        //         h *= 2;
        //     }
        // }
        //double diff = std::abs(Yi[Yi.size() - 1].back() - Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 6]);
        //std::cout << "T1 = " << Yi[Yi.size() - 1].back() << "\n T2 = " << Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 2] << "\n";
        //std::cout << "h = " << h << "\n";
        // if (stepNum % 5 == 0 && diff > 1.0 && h > 1e-9) {
        //     if (oldH < h) {
        //         h = oldH * 0.9;
        //     } else {
        //         h *= 0.9;
        //     }
        //     std::cout << "T1 = " << Yi[Yi.size() - 1].back() << "\n T2 = " << Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 2] << "\n";
        //     std::cout << "step " << stepNum << ": h * 0.75 = " << h << "\n";
        // } else if (currLen < h_last / 2 && stepNum % 5 == 0 && diff < 1.e-6 && h < 1e-8) {
        //     h *= 1.1;
        //     std::cout << "T1 = " << Yi[Yi.size() - 1].back() << "\n T2 = " << Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 2] << "\n";
        //     std::cout << "step " << stepNum << ": h * 1.25 = " << h << "\n";
        // }
        //инициализация аргументов
        for (uint64_t k = 0; k < args.size(); ++k) {
           args[k] = Yi[k][1];
        }
        double maxSpeed = 0;
        for (uint64_t k = 0; k < func.size() - 2; ++k) {
            maxSpeed = std::max(std::abs(func[k](args)), maxSpeed);
        }
        currLen += h;
        if (maxSpeed > 1.e+5 && h > h_min  && period != POST_REACT) {
            //if (oldH < h) {
            h = h_min;
            //} else {
            //    h *= 0.9;
            //}
            if (period == PRE_REACT) {
                period = ACTIVE_REACT;
            }
            //std::cout << "T1 = " << Yi[Yi.size() - 1].back() << "\n T2 = " << Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 2] << "\n";
            //std::cout << "speed: " << maxSpeed << "\n";
            //std::cout << "step " << stepNum << ": h * 0.9 = " << h << "\n";
        } else if (maxSpeed < 1.e+2 && h < h_max && period != ACTIVE_REACT) {
            h *= 1.1;
            if (h > h_max) {
                h = h_max;
            }
            if (period == ACTIVE_REACT) {
                period = POST_REACT;
            }
            //std::cout << "T1 = " << Yi[Yi.size() - 1].back() << "\n T2 = " << Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 2] << "\n";
            //std::cout << "speed: " << maxSpeed << "\n";
            //std::cout << "step " << stepNum << ": h * 1.1 = " << h << "\n";
        }
        if (stepNum % 200 == 0) {
            std::cout << "\nsolution: " << currLen / h_last * 100 << "%\nsteps: " << stepNum << ": h = " << h << "\n";
            std::cout << "T1 = " << Yi[Yi.size() - 1].back() << "\n T2 = " << Yi[Yi.size() - 1][Yi[Yi.size() - 1].size() - 2] << "\n";
            std::cout << "max speed: " << maxSpeed << "\n";
            output << stepNum << "\t";
            for (uint64_t j = 0; j < Yi.size() - 1; ++j) {
                Ytough[j] = Yi[j][1];
                Yprev[j] = Yi[j][0];
                output << Yi[j][1] << "\t";
            }
            output << Yi.back()[1] << "\n";
            //exit(0);
            //tough_coeff = ToughCoeff(func, Ytough, Yprev);
            //std::cout << "Tough coeff: " << tough_coeff << "\n";
            //tough_output << tough_coeff << "\n";
        }
        //окончание вычисления
        // bool stop = true;
        // for (uint64_t i = 1; i < Yi.size() - 2; ++i) {
        //     if (std::abs(Yi[i][stepNum] - Yi[i][stepNum - 1]) > 1e-6) {
        //         stop = false;
        //         break;
        //     }
        // }
        // if (stop && ) {
        //     break;
        // }
        ++stepNum;
        // if (stepNum == 2) {
        //     stepNum = 1;
        // }
        if (currLen > h_last || maxSpeed < 1.e-3) {
            break;
        }
    }
    output.close();
    tough_output.close();
    //std::cout << "Total steps: " << stepNum << "\n";
    for (uint64_t i = 0; i < Y.size(); ++i) {
        Y[i] = Yi[i].back();
    }
    ////std::cout << "\n";
    for (uint64_t i = 0; i < task.getCount(); ++i) {
        double qin = 0, qout = 0;
        auto gammaIn = task[i].getInput();
        auto gammaOut = task[i].getOutput();
        for (uint64_t j = 0; j < func.size() - 2; ++j) {
            qin += nu[j](Y) * gammaIn[j];
            qout += nu[j](Y) * gammaOut[j];
        }
    }
    T = Yi.back().back();
    U = getU(Ui, Y, T);
    auto subs = task.getSubstanceList();
    auto subs_mass = task.getSubstanceMasses();
    double rho_gamma = 0;
    for (uint64_t j = 0; j < Yi.size() - 3; ++j) {
        rho_gamma += Yi[j + 1].back();// * subs_mass.find(subs[j])->second;
    }
    // for (uint64_t j = 0; j < gamma.size(); ++j) {
    //     rho += gamma[j];
    // }
    // return P / (R * T * rho);
    for (uint64_t i = 0; i < Yi.size() - 3; ++i) {
        double curr_conc = Yi[i + 1].back();
        curr_conc *= rho_gamma;
    }
    double lastP = Yi[Yi.size() - 2][1] * (R * Yi[Yi.size() - 1][1] * rho_gamma);
    //std::cout << "U begin: " << Ustart << "\n";
    //std::cout << "U end: "   << U      << "\n";
    //std::cout << "P end: " << lastP << "\n";
    //std::cout << "T end: " << T << "\n";
    return Yi;
}