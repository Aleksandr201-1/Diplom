#include "IterationAlgo.hpp"

float128_t norma (const std::vector<float128_t> &a, const std::vector<float128_t> &b) {
    uint64_t size = std::min(a.size(), b.size());
    float128_t ans = 0;
    for (uint64_t i = 0; i < size; ++i) {
        ans = std::max(ans, std::abs(a[i] - b[i]));
    }
    return ans;
}

float128_t derivative (const std::function<float128_t (const std::vector<float128_t> &, const std::vector<std::vector<float128_t>> &, const Matrix<float128_t> &)> &f, const std::vector<float128_t> &X, const std::vector<std::vector<float128_t>> &K, const Matrix<float128_t> &butcher, float128_t h, uint64_t idx, uint64_t idy) {
    //std::vector<float128_t> args(X);
    //Matrix<float128_t> tmp(butcher);
    std::vector<std::vector<float128_t>> Ktmp(K);
    float128_t a, b, c, d;
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

std::vector<std::vector<float128_t>> ExplicitStep (const std::vector<std::vector<float128_t>> &K, const std::vector<std::vector<float128_t>> &Yi, uint64_t idx, const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f, const Matrix<float128_t> &butcher, float128_t h) {
    std::vector<std::vector<float128_t>> next = K;
    uint64_t orderOfTask = next.size(), orderOfApprox = next[0].size();
    std::vector<float128_t> args(Yi.size(), 0);
    for (uint64_t j = 0; j < args.size(); ++j) {
        args[j] = Yi[j][idx];
    }
    //std::cout << args.size() << " " << orderOfTask << " " << orderOfApprox << "\n";
    //exit(0);
    for (uint64_t j = 0; j < orderOfApprox; ++j) {
        //аргументы для функций
        args[0] = Yi[0][idx] + butcher(j, 0) * h;
        for (uint64_t k = 1; k < orderOfTask + 1; ++k) {
            args[k] = Yi[k][idx];
            for (uint64_t l = 0; l < orderOfApprox; ++l) { //изменить orderOfApprox на j?
                args[k] += butcher(j, l + 1) * next[k - 1][l];
            }
        }
        //exit(0);
        //коэффициенты
        // for (auto el : args) {
        //     std::cout << el << "\n";
        // }
        // for (uint64_t i = 0; i < f.size(); ++i) {
        //     std::cout << "\tFunc " << i + 1 << ": " << f[i](args) << "\n";
        // }
        ///exit(0);
        for (uint64_t k = 0; k < orderOfTask; ++k) {
            next[k][j] = h * f[k](args);
        }
    }
    //exit(0);
    return next;
}

std::vector<std::vector<float128_t>> SimpleIteration (const std::vector<std::vector<float128_t>> &K, const std::vector<std::vector<float128_t>> &Yi, uint64_t idx, const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f, const Matrix<float128_t> &butcher, float128_t h, float128_t approx, IterationAlgo algo) {
    //std::cout << "eeee\n";
    std::vector<float128_t> args(Yi.size(), 0);
    for (uint64_t j = 0; j < args.size(); ++j) {
        args[j] = Yi[j][idx];
    }
    std::vector<std::vector<float128_t>> prev = K, next = K;
    uint64_t orderOfTask = f.size(), orderOfApprox = butcher.size().m - 1;
    std::vector<float128_t> new_args;
    float128_t prevDiff, diff = INFINITY;
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

std::vector<std::vector<float128_t>> NewtonIteration (const std::vector<std::vector<float128_t>> &K, const std::vector<std::vector<float128_t>> &Yi, uint64_t idx, const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f, const Matrix<float128_t> &butcher, float128_t h, float128_t approx) {
    std::vector<std::vector<float128_t>> prev = K, next = K;
    uint64_t orderOfTask = K.size(), orderOfApprox = butcher.size().m - 1;
    std::vector<float128_t> args(Yi.size(), 0);
    for (uint64_t j = 0; j < args.size(); ++j) {
        args[j] = Yi[j][idx];
    }
    std::vector<float128_t> new_args(args.size(), 0);
    float128_t prevDiff, diff = INFINITY;
    uint64_t iter = 0;
    Matrix<std::function<float128_t (const std::vector<float128_t> &, const std::vector<std::vector<float128_t>> &, const Matrix<float128_t> &)>> F(orderOfTask, orderOfApprox);
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
            //std::vector<std::vector<std::function<float128_t (const std::vector<float128_t> &, const std::vector<std::vector<float128_t>> &, const Matrix<float128_t> &)>>> F(orderOfTask);
            for (uint64_t j = 0; j < orderOfTask; ++j) {
                //args.size = orderOfTask * orderOfApprox + orderOfTask + 2
                // std::vector<float128_t> Kargs(X.size(), 0);
                // for (uint64_t k = 0; k < X.size(); ++k) {
                //     Kargs[k] = new_args[k];
                //     for (uint64_t l = 0; l < X[0].size(); ++l) {
                //         Kargs[k] += prev[k][l] * butcher(k, l);
                //     }
                // }
                auto lambda = [=] (const std::vector<float128_t> &args, const std::vector<std::vector<float128_t>> &K, const Matrix<float128_t> &butcher) -> float128_t {
                    //float128_t ans = butcher(j, i);
                    float128_t ans = K[j][i];
                    std::vector<float128_t> arg = args;
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
        Matrix<float128_t> Yakobi(orderOfApprox * orderOfTask);
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

std::vector<std::vector<float128_t>> IterationStep (const std::vector<std::vector<float128_t>> &X, const std::vector<std::vector<float128_t>> &Yi, uint64_t idx, const std::vector<std::function<float128_t (const std::vector<float128_t> &)>> &f, const Matrix<float128_t> &butcher, float128_t h, float128_t approx, IterationAlgo algo) {
    switch (algo) {
        case IterationAlgo::ZEIDEL:
        case IterationAlgo::SIMPLE_ITERATION:
            return SimpleIteration(X, Yi, idx, f, butcher, h, approx, algo);
        case IterationAlgo::NEWTON:
            return NewtonIteration(X, Yi, idx, f, butcher, h, approx);
        default:
            return ExplicitStep(X, Yi, idx, f, butcher, h);
    }
}