#include "KrancNicolas.hpp"

enum ARGS {
    RHO = 0,
    P,
    X,
    V,
    U,
    T,
    Y,
    PSI,
    VI
};

double diff (const std::function<double (const std::vector<double> &)> &func, const std::vector<double> &args, int arg) {
    return 0;
}

//auto ux = [] (const std::vector<double> &args) -> double {};
void sysGen () {
    const std::vector<std::function<double (const std::vector<double> &)>> funcs;
    auto ux = [=] (const std::vector<double> &args) -> double {
        return -1.0 / args[RHO] * diff(funcs[P], args, X) + funcs[V](args) / funcs[U](args) * std::pow(args[Y] / args[PSI], args[VI]) * diff(funcs[P], args, PSI);
    };
}

std::vector<double> GaussSolveSLAE (const Matrix<double> &matrix, const std::vector<double> &ans) {
    Matrix<double> A = matrix;
    std::vector<double> b(ans), x(ans), used(ans.size(), 0);
    uint64_t n = ans.size();
    for (uint64_t i = 0; i < n; ++i) {
        double maxAbsEl = 0;
        uint64_t idx = -1;
        for (uint64_t j = 0; j < n; ++j) {
            if (std::abs(maxAbsEl) < std::abs(A(j, i)) && used[j] == 0) {
                maxAbsEl = A(j, i);
                idx = j;
            }
            //maxAbsEl = std::max(maxAbsEl, std::abs(A(j, i)));
        }
        if (idx == -1) {
            std::cout << "Undefined x\n";
            exit(-1);
        }
        used[idx] = 1;
        for  (uint64_t j = 0; j < n; ++j) {
            if (A(j, i) == 0) {
                continue;
            }
            if (j != idx) {
                double div = A(j, i) / maxAbsEl;
                b[j] -= div * b[idx];
                for (uint64_t k = i; k < n; ++k) {
                    A(j, k) -= div * A(idx, k);
                }
                A(j, i) = 0.0;
            }
        }
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            if (A(i, j) != 0.0) {
                x[j] = b[i] / A(i, j);
                break;
            }
        }
    }
    // for (uint64_t i = 0; i < b.size(); ++i) {
    //     b[i] /= A(i, i);
    // }
    return x;
}

double Point2Order1 (const std::vector<double> &coeff, double h, double u1, double f, uint64_t i) {
    double alpha = coeff[0], beta = coeff[1];
    double ans = 0;
    if (i == 0) {
        ans += f - alpha * u1 / h;
        ans /= beta - alpha / h;
    } else {
        ans += f + alpha * u1 / h;
        ans /= beta + alpha / h;
    }
    return ans;
}

double Point2Order2 (const std::vector<double> &ux, const std::vector<double> &coeff, double h, double t, double u0, double u1, double f, uint64_t i) {
    double alpha = ux[0], beta = ux[1];
    double a = coeff[0], b = coeff[1], c = coeff[2];
    double ans = 0;
    if (i == 0) {
        ans += h / t * u0 - f * (2 * a - b * h) / alpha + 2 * u1 * a / h;
        ans /= 2 * a / h + h / t - c * h - (beta / alpha) * (2 * a - b * h);
    } else {
        ans += h / t * u0 + f * (2 * a + b * h) / alpha + 2 * u1 * a / h;
        ans /= 2 * a / h + h / t - c * h + (beta / alpha) * (2 * a + b * h);
    }
    return ans;
}

double Point3Order2 (const std::vector<double> &coeff, double h, double u1, double u2, double f, uint64_t i) {
    double alpha = coeff[0], beta = coeff[1];
    double ans = 0;
    if (i == 0) {
        ans += f - alpha * (4 * u1 - u2) / (2 * h);
        ans /= beta - 3 * alpha / (2 * h);
    } else {
        ans += f + alpha * (4 * u1 - u2) / (2 * h);
        ans /= beta + 3 * alpha / (2 * h);
    }
    return ans;
}

//вычисление 4-х точек для аппроксимации граничных условий
void FiniteDifferenceExplicitApprox (std::vector<std::vector<double>> &u, double X0, double xh, double th, const std::vector<double> &coeff, const std::function<double(double, double)> &f, uint64_t id) {
    double a = coeff[0], b = coeff[1], c = coeff[2];
    double sigmaA = a * th / (xh * xh), sigmaB = b * th / (2 * xh), sigmaC = c * th;
    uint64_t idx[] = {1, u[id].size() - 3};
    for (uint64_t i = 0; i < 2; ++i) {
        for (uint64_t j = id; j < id + 1; ++j) {
            for (uint64_t k = idx[i]; k < idx[i] + 2; ++k) {
                double tmp = 0;
                tmp += sigmaA * (u[j - 1][k + 1] - 2 * u[j - 1][k] + u[j - 1][k - 1]);
                tmp += sigmaB * (u[j - 1][k + 1] - u[j - 1][k]);
                tmp += sigmaC * u[j - 1][k];
                tmp += f(X0 + k * xh, j * th) * th;
                u[j][k] = tmp + u[j - 1][k];
            }
        }
    }
}

/**
    * \brief Функция одной итерации для явно-неявной схемы
    *
    * \param theta Вес неявной части конечно-разностной схемы. При theta = 1 получаем неявную схему, 
    *   при theta = 0 --- явную, при theta = 0.5 --- схему Кранка-Николаса
    * \param u Таблица функции U(x, t), в которую будет записан ответ.
    * \param ux 2 пары значений коэффициентов для ux и u из левого и правого краевого условия.
    * \param X0 Левая граница (по умолчанию 0).
    * \param xh Размер шага для X.
    * \param th Размер шага для T.
    * \param coeff Коэффициенты при uxx, ux, u из уравнения.
    * \param f Функция-источник из уравнения ut = ux + f(x).
    * \param i номер временного слоя, для которого вычисляется решение
**/
void ExplNonExplIteration (double theta, std::vector<std::vector<double>> &u, const std::vector<std::vector<double>> &ux, double X0, double xh, double th, const std::vector<double> &coeff, const std::function<double(double, double)> &f, uint64_t i) {
    double a = coeff[0], b = coeff[1], c = coeff[2];
    //double sigmaA = a * th / (xh * xh), sigmaB = b * th / (2 * xh), sigmaC = c * th;
    double sigmaA = a * th, sigmaB = b * xh * th, sigmaC = c * xh * xh * th;
    double alpha, beta;


    uint64_t n = u[0].size();
    Matrix<double> M(n, n);
    std::vector<double> ans(n);


    //там, где идёт домножение на theta --- неявная часть, там, где на (1 - theta) --- явная часть
    // M(0, 0) = -1 - (2 * sigmaA - sigmaB + sigmaC) * theta;
    // M(0, 1) = (sigmaA + sigmaB) * theta;
    // ans[0] = -u[i - 1][1] - (sigmaA + sigmaB + sigmaC) * u[i][0] * theta;
    // ans[0] += -(sigmaA * (u[i - 1][2] - 2 * u[i - 1][1] + u[i - 1][0]) + sigmaB * (u[i - 1][2] - u[i - 1][1]) + sigmaC * u[i - 1][1]) * (1 - theta);
    // ans[0] += -f(X0 + xh, i * th) * th;

    // M(0, 0) = 1;
    // ans[0] = u[i][0];

    alpha = ux[0][1];
    beta = ux[0][0];
    M(0, 0) = 2 * xh * alpha - 3 * beta;
    M(0, 1) = 4 * beta;
    M(0, 2) = -beta;
    ans[0] = u[i][0] * 2 * xh;

    for (uint64_t j = 1; j < n - 1; ++j) {
        M(j, j - 1) = -sigmaA * theta;
        M(j, j) = xh * xh + (2 * sigmaA + sigmaB - sigmaC) * theta;
        M(j, j + 1) = -(sigmaA + sigmaB) * theta;
        ans[j] = u[i - 1][j - 1] * sigmaA * (1.0 - theta) + 
                 u[i - 1][j] * (xh * xh + (-2 * sigmaA - sigmaB + sigmaC) * (1.0 - theta)) +
                 u[i - 1][j + 1] * (sigmaA + sigmaB) * (1.0 - theta) +
                 f(X0 + j * xh, i * th) * xh * xh * th;
        // M(j, j - 1) = sigmaA * theta;
        // M(j, j) = -1 - (2 * sigmaA - sigmaB + sigmaC) * theta;
        // M(j, j + 1) = (sigmaA + sigmaB) * theta;
        // ans[j] = -u[i - 1][j + 1];
        // ans[j] += -(sigmaA * (u[i - 1][j + 2] - 2 * u[i - 1][j + 1] + u[i - 1][j]) + sigmaB * (u[i - 1][j + 2] - u[i - 1][j + 1]) + sigmaC * u[i - 1][j + 1]) * (1 - theta);
        // ans[j] += -f(X0 + (j + 1) * xh, i * th) * th;
    }

    // M(n - 1, n - 2) = sigmaA * theta;
    // M(n - 1, n - 1) = -1 - (2 * sigmaA - sigmaB + sigmaC) * theta;
    // ans[n - 1] = -u[i - 1][n] - (sigmaA + sigmaB + sigmaC) * u[i][n + 1] * theta;
    // ans[n - 1] += -(sigmaA * (u[i - 1][n + 1] - 2 * u[i - 1][n] + u[i - 1][n - 1]) + sigmaB * (u[i - 1][n + 1] - u[i - 1][n]) + sigmaC * u[i - 1][n]) * (1 - theta);
    // ans[n - 1] += -f(X0 + xh * n, i * th) * th;

    // M(n - 1, n - 1) = 1;
    // ans[n - 1] = u[i][n - 1];

    alpha = ux[1][1];
    beta = ux[1][0];
    M(n - 1, n - 3) = beta;
    M(n - 1, n - 2) = -4 * beta;
    M(n - 1, n - 1) = 2 * xh * alpha + 3 * beta;
    ans[n - 1] = u[i][n - 1] * 2 * xh;

    // double div;
    // div = M(0, 2) / M(1, 2);
    // M(0, 0) -= div * M(1, 0);
    // M(0, 1) -= div * M(1, 1);
    // M(0, 2) = 0;
    // ans[0] -= div * ans[1];

    // div = M(n - 1, n - 3) / M(n - 2, n - 3);
    // M(n - 1, n - 1) -= div * M(n - 2, n - 1);
    // M(n - 1, n - 2) -= div * M(n - 2, n - 2);
    // M(n - 1, n - 3) = 0;
    // ans[n - 1] -= div * ans[n - 2];

    //если theta = 0, то решать СЛАУ не нужно
    if (theta == 0.0) {
        // for (uint64_t j = 0; j < ans.size(); ++j) {
        //     ans[j] *= -1;
        // }
    } else {
        //ans = LUsolveSLAE(M, ans);
        ans = GaussSolveSLAE(M, ans);
        //ans = SISolveSLAE(M, ans, 0.1);
        //ans = RUNsolveSLAE(M, ans);
    }

    //записываем в следующий временной слой решение
    for (uint64_t j = 0; j < ans.size(); ++j) {
        u[i][j] = ans[j];
    }
}

/**
    * \brief Неявно-явная конечно-разностная схема решения краевой задачи
    *
    * \param theta Вес неявной части конечно-разностной схемы. При theta = 1 получаем неявную схему, 
    *   при theta = 0 --- явную, при theta = 0.5 --- схему Кранка-Николаса
    * \param u Таблица функции U(x, t), в которую будет записан ответ.
    * \param ux 2 пары значений коэффициентов для ux и u из левого и правого краевого условия.
    * \param X0 Левая граница (по умолчанию 0).
    * \param xh Размер шага для X.
    * \param th Размер шага для T.
    * \param coeff Коэффициенты при uxx, ux, u из уравнения.
    * \param f Функция-источник из уравнения ut = ux + f(x).
    * \param left Метод аппроксимации производной на левой границе.
    * \param right Метод аппроксимации производной на правой границе.
    *
**/
void ExplNonExpl (double theta, std::vector<std::vector<double>> &u, const std::vector<std::vector<double>> &ux, double X0, double xh, double th, const std::vector<double> &coeff, const std::function<double(double, double)> &f, ApproxLevel left, ApproxLevel right) {
    if (theta < 0 || theta > 1) {
        throw std::logic_error("ExplNonExpl: theta must be in range [0, 1]");
    }
    for (uint64_t i = 1; i < u.size(); ++i) {
        //применяя явную конечно-разностную схему находим 2 точки справа от левой границы и слева от правой
        //FiniteDifferenceExplicitApprox(u, X0, xh, th, coeff, f, i);
        uint64_t start = 0, end = u[i].size() - 1;
        //применение аппроксимации, если нужно, для левой и правой границы, используя точки, найденые ранее
        // switch (left) {
        //     case ApproxLevel::POINT2_ORDER1:
        //         u[i][start] = Point2Order1(ux[0], xh, u[i][start + 1], u[i][start], 0);
        //         break;
        //     case ApproxLevel::POINT2_ORDER2:
        //         u[i][start] = Point2Order2(ux[0], coeff, xh, th, u[i - 1][start], u[i][start + 1], u[i][start], 0);
        //         break;
        //     case ApproxLevel::POINT3_ORDER2:
        //         u[i][start] = Point3Order2(ux[0], xh, u[i][start + 1], u[i][start + 2], u[i][start], 0);
        //         break;
        //     default:
        //         break;
        // }
        // switch (right) {
        //     case ApproxLevel::POINT2_ORDER1:
        //         u[i][end] = Point2Order1(ux[1], xh, u[i][end - 1], u[i][end], 0);
        //         break;
        //     case ApproxLevel::POINT2_ORDER2:
        //         u[i][end] = Point2Order2(ux[1], coeff, xh, th, u[i - 1][end], u[i][end - 1], u[i][end], 0);
        //         break;
        //     case ApproxLevel::POINT3_ORDER2:
        //         u[i][end] = Point3Order2(ux[1], xh, u[i][end - 1], u[i][end - 2], u[i][end], 0);
        //         break;
        //     default:
        //         break;
        // }
        //используем итерацию для вычисления следующего временного слоя
        ExplNonExplIteration(theta, u, ux, X0, xh, th, coeff, f, i);
    }
}

std::vector<std::vector<double>> SolveIBVP (const Task &task, double timeLimit, double xh, double th, Method method, ApproxLevel approx) {
    double X1 = task.X[0], X2 = task.X[1];
    //задаём размер таблицы функции u
    std::vector<std::vector<double>> u(uint64_t(timeLimit / th), std::vector<double>(uint64_t((X2 - X1) / xh), 0));
    ApproxLevel left = ApproxLevel::NONE, right = ApproxLevel::NONE;
    //запоминаем аппроксимацию граничных условий, если они второго или третьего рода
    if (task.ux[0][0] != 0) { //если коэффициент при ux в краевом условии не 0, то либо 2й, либо 3й род
        left = approx;
    }
    if (task.ux[1][0] != 0) {
        right = approx;
    }
    //уравнение ut = a * uxx + b * ux + c * u + f(x)
    //возвращает значение f(x) из правой части
    auto f = [&] (double x, double t) -> double {
        return task.trees[0]({0, 0, 0, x, t});
    };
    //функция распределения на левой границе
    auto fx0 = [&] (double t) -> double {
        return task.trees[1](t);
    };
    //функция распределения на правой границе
    auto fxl = [&] (double t) -> double {
        return task.trees[2](t);
    };
    //функция распределения в начальный момент времени
    auto ft0 = [&] (double x) -> double {
        return task.trees[3](x);
    };
    //инициализация начальных условий
    for (uint64_t i = 0; i < u[0].size(); ++i) {
        u[0][i] = ft0(i * xh);
    }
    for (uint64_t i = 1; i < u.size(); ++i) {
        u[i].front() = fx0(i * th);
        u[i].back()  = fxl(i * th);
    }
    //решения задачи в зависимости от метода и метода аппроксимации
    switch (method) {
        case Method::EXPLICIT:
            ExplNonExpl(0, u, task.ux, X1, xh, th, task.coeff, f, left, right); //theta = 0, получаем явную схему
            break;
        case Method::NOT_EXPLICIT:
            ExplNonExpl(1, u, task.ux, X1, xh, th, task.coeff, f, left, right); //theta = 1, получаем неявную схему
            break;
        case Method::KRANK_NICOLAS:
            ExplNonExpl(0.5, u, task.ux, X1, xh, th, task.coeff, f, left, right); //theta = 0.5, получаем схему Кранка-Николоса
            break;
        default:
            break;
    }
    return u;
}