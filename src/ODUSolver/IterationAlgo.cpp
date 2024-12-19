#include "IterationAlgo.hpp"

const std::map<IterationAlgo, std::string> iteration_algos = {
    {IterationAlgo::NEWTON, "Newton"},
    {IterationAlgo::ZEIDEL, "Zeidel"},
    {IterationAlgo::SIMPLE_ITERATION, "SI"},
    {IterationAlgo::EXPLICIT_STEP, "Explicit"}
};

std::string IterationAlgoToString (IterationAlgo algo) {
    return enumToString(algo, iteration_algos);
}

IterationAlgo stringToIterationAlgo (const std::string &str) {
    return stringToEnum(str, iteration_algos);
}

double norma (const std::vector<double> &a, const std::vector<double> &b) {
    uint64_t size = std::min(a.size(), b.size());
    double ans = 0;
    for (uint64_t i = 0; i < size; ++i) {
        ans = std::max(ans, std::abs(a[i] - b[i]));
    }
    return ans;
}

double derivative (const std::function<double (const std::vector<double> &, const std::vector<std::vector<double>> &, const Matrix<double> &)> &f, const std::vector<double> &X, const std::vector<std::vector<double>> &K, const Matrix<double> &butcher, double h, uint64_t idx, uint64_t idy) {
    //std::vector<double> args(X);
    //Matrix<double> tmp(butcher);
    std::vector<std::vector<double>> Ktmp(K);
    double a, b, c, d;
    //args[idx] += 2 * h;
    Ktmp[idx][idy] += 2 * h;
    a = f(X, Ktmp, butcher);
    Ktmp[idx][idy] -= h;
    b = f(X, Ktmp, butcher);
    Ktmp[idx][idy] -= 2 * h;
    c = f(X, Ktmp, butcher);
    Ktmp[idx][idy] -= h;
    d = f(X, Ktmp, butcher);
    return (-a + 8 * b - 8 * c + d) / (12 * h);
}

std::vector<std::vector<double>> ExplicitStep (const std::vector<std::vector<double>> &K, const std::vector<std::vector<double>> &Yi, uint64_t idx, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h) {
    std::vector<std::vector<double>> next = K;
    uint64_t orderOfTask = next.size(), orderOfApprox = next[0].size();
    std::vector<double> args(Yi.size(), 0);
    for (uint64_t j = 0; j < args.size(); ++j) {
        args[j] = Yi[j][idx];
    }
    //std::cout << args.size() << " " << orderOfTask << " " << orderOfApprox << "\n";
    for (uint64_t j = 0; j < orderOfApprox; ++j) {
        //аргументы для функций
        args[0] = Yi[0][idx] + butcher(j, 0) * h;
        for (uint64_t k = 1; k < orderOfTask + 1; ++k) {
            args[k] = Yi[k][idx];
            for (uint64_t l = 0; l < j; ++l) { //изменить orderOfApprox на j?
                args[k] += butcher(j, l + 1) * next[k - 1][l];
            }
        }
        //коэффициенты
        // for (auto el : args) {
        //     std::cout << el << "\n";
        // }
        // for (uint64_t i = 0; i < f.size(); ++i) {
        //     std::cout << "\tFunc " << i + 1 << ": " << f[i](args) << "\n";
        // }
        for (uint64_t k = 0; k < orderOfTask; ++k) {
            next[k][j] = h * f[k](args);
        }
    }
    return next;
}

std::vector<std::vector<double>> SimpleIteration (const std::vector<std::vector<double>> &K, const std::vector<std::vector<double>> &Yi, uint64_t idx, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx, IterationAlgo algo) {
    //std::cout << "eeee\n";
    std::vector<double> args(Yi.size(), 0);
    for (uint64_t j = 0; j < args.size(); ++j) {
        args[j] = Yi[j][idx];
    }
    std::vector<std::vector<double>> prev = K, next = K;
    uint64_t orderOfTask = f.size(), orderOfApprox = butcher.size().m - 1;
    std::vector<double> new_args;
    double prevDiff, diff = INFINITY;
    uint64_t iter = 0;
    //std::cout << "eeee\n";
    do {
        prevDiff = diff;
        diff = 0;
        prev = next;
        for (uint64_t j = 0; j < orderOfApprox; ++j) {
            new_args = args;
            //аргументы для функций
            new_args[0] += butcher(j, 0) * h;
            for (uint64_t k = 1; k < orderOfTask + 1; ++k) {
                for (uint64_t l = 0; l < orderOfApprox; ++l) {
                    if (algo == IterationAlgo::SIMPLE_ITERATION) {
                        new_args[k] += butcher(j, l + 1) * prev[k - 1][l] * h;
                    } else {
                        new_args[k] += butcher(j, l + 1) * next[k - 1][l] * h;
                    }
                }
            }
            //коэффициенты
            for (uint64_t k = 0; k < orderOfTask; ++k) {
                next[k][j] = f[k](new_args);
            }
        }
        //std::cout << "eeee\n";
        for (uint64_t i = 0; i < prev.size(); ++i) {
            diff = std::max(diff, norma(prev[i], next[i]));
        }
        //std::cout << "diff: " << diff << "\n";
        //std::cout << "iter num: " << iter << "\n";
        ++iter;
        if (prevDiff < diff) {
            return prev;
        }
    } while (diff > approx);
    //std::cout << "Count of iteration: " << iter << "\n";
    return next;
}

std::vector<std::vector<double>> NewtonIteration (const std::vector<std::vector<double>> &K, const std::vector<std::vector<double>> &Yi, uint64_t idx, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx) {
    std::vector<std::vector<double>> prev = K, next = K;
    uint64_t orderOfTask = K.size(), orderOfApprox = butcher.size().m - 1;
    std::vector<double> args(Yi.size(), 0);
    for (uint64_t j = 0; j < args.size(); ++j) {
        args[j] = Yi[j][idx];
    }
    std::vector<double> new_args(args.size(), 0);
    double prevDiff, diff = INFINITY;
    uint64_t iter = 0;
    Matrix<std::function<double (const std::vector<double> &, const std::vector<std::vector<double>> &, const Matrix<double> &)>> F(orderOfTask, orderOfApprox);
    do {
        prevDiff = diff;
        diff = 0;
        prev = next;
        for (uint64_t i = 0; i < orderOfApprox; ++i) {
            //аргументы для функций
            //new_args = args;
            new_args[0] = butcher(i, 0) * h;
            for (uint64_t k = 1; k < orderOfTask + 1; ++k) {
                new_args[k] = 0;
                for (uint64_t l = 0; l < orderOfApprox; ++l) {
                    new_args[k] += butcher(i, l + 1) * next[k - 1][l] * h;
                }
            }
            //std::cout << "size: " << X.size() << "\n";
            //std::cout << "approx size: " << orderOfApprox << "\n";
            //std::cout << "task order size: " << orderOfTask << "\n";
            //std::cout << "eeee\n";
            //F
            //std::vector<std::vector<std::function<double (const std::vector<double> &, const std::vector<std::vector<double>> &, const Matrix<double> &)>>> F(orderOfTask);
            for (uint64_t j = 0; j < orderOfTask; ++j) {
                //args.size = orderOfTask * orderOfApprox + orderOfTask + 2
                // std::vector<double> Kargs(X.size(), 0);
                // for (uint64_t k = 0; k < X.size(); ++k) {
                //     Kargs[k] = new_args[k];
                //     for (uint64_t l = 0; l < X[0].size(); ++l) {
                //         Kargs[k] += prev[k][l] * butcher(k, l);
                //     }
                // }
                auto lambda = [=] (const std::vector<double> &args, const std::vector<std::vector<double>> &K, const Matrix<double> &butcher) -> double {
                    //double ans = butcher(j, i);
                    double ans = K[j][i];
                    std::vector<double> arg = args;
                    arg[0] += butcher(i, 0) * h;
                    for (uint64_t k = 1; k < orderOfTask + 1; ++k) {
                        //arg[k] = 0;
                        for (uint64_t l = 0; l < orderOfApprox; ++l) {
                            arg[k] += butcher(i, l + 1) * next[k - 1][l] * h;
                        }
                    }
                    // for (uint64_t k = 0; k < X.size(); ++k) {
                    //     arg[k] += args[k];
                    // }
                    //for (ui)
                    return ans - f[j](arg);
                };
                F(j, i) = lambda;
                //F[j] = lambda;
            }
        }
        //std::cout << "eeee\n";
        //матрица Якоби
        Matrix<double> Yakobi(orderOfApprox * orderOfTask);
        for (uint64_t i = 0; i < orderOfApprox; ++i) {
            for (uint64_t j = 0; j < orderOfTask; ++j) {
                for (uint64_t k = 0; k < orderOfApprox; ++k) {
                    for (uint64_t l = 0; l < orderOfTask; ++l) {
                        Yakobi(orderOfTask*i + j, orderOfTask*k + l) = derivative(F(j, i), new_args, K, butcher, h, l, k);
                    }
                }
            }
        }
        //std::cout << "Yakobi matrix:\n" << Yakobi << "\n";
        //std::cout << "eeee\n";
        //реверс Матрицы Якоби
        Yakobi = LUReverseMatrix(Yakobi);
        //std::cout << "eeee\n";
        //новые коэффициенты
        for (uint64_t i = 0; i < orderOfTask; ++i) {
            for (uint64_t j = 0; j < orderOfApprox; ++j) {
                for (uint64_t k = 0; k < orderOfTask; ++k) {
                    for (uint64_t l = 0; l < orderOfApprox; ++l) {
                        next[i][j] -= Yakobi(orderOfApprox*i + j, orderOfApprox*k + l) * F(k, l)(args, prev, butcher);
                    }
                }
                //next[k][i] = f[k](new_args);
            }
        }
        //std::cout << "eeee\n";
        
        for (uint64_t i = 0; i < prev.size(); ++i) {
            diff = std::max(diff, norma(prev[i], next[i]));
        }
        if (prevDiff < diff) {
            return prev;
        }
        ++iter;
    } while (diff > approx);
    //std::cout << "Count of iteration: " << iter << "\n";
    return next;
}

std::vector<std::vector<double>> IterationStep (const std::vector<std::vector<double>> &X, const std::vector<std::vector<double>> &Yi, uint64_t idx, const std::vector<std::function<double (const std::vector<double> &)>> &f, const Matrix<double> &butcher, double h, double approx, IterationAlgo algo) {
    switch (algo) {
        case IterationAlgo::ZEIDEL:
        case IterationAlgo::SIMPLE_ITERATION:
            return SimpleIteration(X, Yi, idx, f, butcher, h, approx, algo);
        case IterationAlgo::NEWTON:
            return NewtonIteration(X, Yi, idx, f, butcher, h, approx);
        case IterationAlgo::EXPLICIT_STEP:
            return ExplicitStep(X, Yi, idx, f, butcher, h);
        default:
            return {};
    }
}