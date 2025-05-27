#include "CrankNicolson.hpp"

double HalfDivFind (const std::function<double (double)> &f, double a, double b, double approx) {
    if (f(a) * f(b) > 0) {
        //throw std::runtime_error("HalfDivFind: Incorrect begining point");
    }
    uint64_t count = 0;
    double x = (a + b) / 2;
    while (f(x) > approx) {
        if (f(x) * f(a) > 0) {
            a = x;
        }
        if (f(x) * f(b) > 0) {
            b = x;
        }
        x = (a + b) / 2;
        ++count;
        if (count > 100) {
            return x;
        }
    }
    return x;
}

const std::map<ARGS, std::string> strToArgsMap = {
    {U, "U"},
    {W, "W"},
    {J, "J"},
    {RHO, "RHO"},
    {Q, "Q"},
    {P, "P"},
    {V, "V"},
    {T, "T"},
    {Y, "Y"},
    {MU, "MU"},
    {C1, "C1"}
};

int strToArgs (const std::string &str) {
    return stringToEnum(str, strToArgsMap);
}

std::string argToStr (ARGS arg) {
    return enumToString(arg, strToArgsMap);
}

float LinearFunc (const std::vector<std::pair<double, double>> &points, float x, int method = 0) {
    float ans = 0;
    uint64_t i = 0;
    if (method == 0) {
        while (x > points[i + 1].first) {
            ++i;
            if (i + 1 == points.size()) {
                --i;
                break;
            }
        }
        float coeff = (points[i + 1].second - points[i].second) / (points[i + 1].first - points[i].first);
        return (points[i].second - points[i].first * coeff) + coeff * x;
    } else {
        while (x > points[i + 1].second) {
            ++i;
            if (i + 1 == points.size()) {
                --i;
                break;
            }
        }
        float coeff = (points[i + 1].first - points[i].first) / (points[i + 1].second - points[i].second);
        return (points[i].first - points[i].second * coeff) + coeff * x;
    }
}

double ValueProfile::toPhysical (double val) const {
    return val * in;
}

double ValueProfile::toNormal (double val) const {
    return val / in;
}

double ValueProfile::getPhysical (double x) const {
    return toPhysical(LinearFunc(profile, x));
}

ValueProfile::ValueProfile () {
    in = 0;
}

ValueProfile::ValueProfile (uint64_t size) {
    in = 0;
    profile.resize(size, {0, 0});
}

ValueProfile::ValueProfile (double in, uint64_t size) : in(in) {
    profile.resize(size, {0, 0});
}

double ChemicalProfile::toPhysical (double val) const {
    return in == out ? in : std::abs(val * (in - out) + out);
    //return val * in + (1.0 - val) * out;
}

double ChemicalProfile::toNormal (double val) const {
    return in == out ? val / in : std::abs((val - out) / (in - out));
}

double ChemicalProfile::getPhysical (double x) const {
    return toPhysical(LinearFunc(profile, x));
}

double ChemicalProfile::getNormal (double val) const {
    return LinearFunc(profile, val);
}

ChemicalProfile::ChemicalProfile () {
    in = out = 0;
}

ChemicalProfile::ChemicalProfile (uint64_t size) {
    in = out = 0;
    profile.resize(size, {0, 0});
}

ChemicalProfile::ChemicalProfile (double in, double out, uint64_t size) : in(in), out(out) {
    profile.resize(size, {0, 0});
}

uint64_t findDeltaU (const std::vector<Matrix<double>> &f, uint64_t i, double value) {
    uint64_t pos = 0;
    double prevDiff = 1, currDiff = 1;
    auto size = f[U].size();
    //double U0 = f[U](0, 0), Un = f[U](0, size.m - 1);
    while (currDiff > 0) {
        prevDiff = currDiff;
        currDiff = f[U](i, pos);
        //currDiff = (f[U](i, pos) - Un) / (U0 - Un);
        currDiff = std::abs(currDiff - value);
        if (currDiff > prevDiff) {
            --pos;
            break;
        }
        if (pos == size.m - 1) {
            break;
        }
        ++pos;
    }
    //exit(0);
    return pos;
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
void ExplNonExplIteration (ReportInfo &info, double alpha, double R0, std::vector<Matrix<double>> &f, double Xi_0, double x_h, double xi_h, uint64_t i, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs) {
    const double nu = 1;
    const double Sc = 0.7, Pr = 0.7, Le = Pr / Sc;
    uint64_t chemCount = f.size() - C1;
    uint64_t n = f[0].size().m;
    Matrix<double> M(n, n);
    std::vector<double> ans(n);
    Matrix<double> Mw(n - 1, n - 1);
    std::vector<double> answ(n - 1);
    const double &L = x_h, &H = xi_h;
    bool fireAvail = false;

    auto comp = [&] (const std::vector<uint64_t> &numerator, const std::vector<uint64_t> &denominator, uint64_t n, uint64_t m) -> double {
        double ans = 1.0;
        for (uint64_t i : numerator) {
            ans *= f[i](n, m);
        }
        for (uint64_t i : denominator) {
            ans /= f[i](n, m);
        }
        return ans;
    };
    auto writeAns = [&] (uint64_t funcIdx) -> void {
        for (uint64_t j = 0; j < ans.size(); ++j) {
            f[funcIdx](i, j) = ans[j];
        }
    };
    auto yp = [&] (uint64_t idx, uint64_t jdx) -> double {
        double xi = (jdx + 1) * xi_h + Xi_0;
        return std::pow(f[Y](idx, jdx) / xi, nu);
    };
    auto yp2 = [&] (uint64_t idx, uint64_t jdx) -> double {
        double xi = (jdx + 1) * xi_h + Xi_0;
        return std::pow(f[Y](idx, jdx) * f[Y](idx, jdx) / xi, nu);
    };
    auto F1 = [&] (uint64_t n, uint64_t m) -> double {
        double xi = (m + 1) * xi_h + Xi_0;
        return f[V](n, m) / f[U](n, m) * std::pow(f[Y](n, m) / xi, nu);
    };
    auto F2 = [&] (uint64_t idx, uint64_t jdx) -> double {
        double xi = (jdx + 1) * xi_h + Xi_0;
        if (jdx == f[MU].size().m - 1) {
            jdx -= 1;
        }
        return (comp({MU, RHO, U}, {}, idx, jdx + 1) * yp2(idx, jdx + 1) - comp({MU, RHO, U}, {}, idx, jdx) * yp2(idx, jdx)) / (2 * H);
    };
    auto F3 = [&] (uint64_t n, uint64_t m) -> double {
        double xi = (m + 1) * xi_h + Xi_0;
        if (m == f[MU].size().m - 1) {
            m -= 1;
        }
        return (comp({MU, RHO, U}, {}, n, m + 1) * yp(n, m + 1) - comp({MU, RHO, U}, {}, n, m) * yp(n, m)) / (2 * H);
    };
    auto F4 = [&] (uint64_t n, uint64_t m) -> double {
        double xi = (m + 1) * xi_h + Xi_0;
        if (m == f[MU].size().m - 1) {
            m -= 1;
        }
        return (comp({MU, RHO, U, U}, {}, n, m + 1) * yp2(n, m + 1) - comp({MU, RHO, U, U}, {}, n, m) * yp2(n, m)) / (2 * H);
    };
    auto printIterData = [&] (uint64_t j) {
        if (j > C1) {
            std::cout << "\n" << i << " before C" << j - C1 + 1 << ":\n";
        } else {
            std::cout << "\n" << i << " before " << argToStr(ARGS(j)) << ":\n";
        }
        for (uint64_t k = 0; k < f[j].size().m; ++k) {
            std::cout << std::setw(13) << f[j](i - 1, k) << " ";
        }
        if (j > C1) {
            std::cout << "\n" << i << " after C" << j - C1 + 1 << ":\n";
        } else {
            std::cout << "\n" << i << " after " << argToStr(ARGS(j)) << ":\n";
        }
        for (uint64_t k = 0; k < f[j].size().m; ++k) {
            std::cout << std::setw(13) << f[j](i, k) << " ";
        }
        std::cout << "\n";
    };

    for (uint64_t j = 0; j < chemProfs.size(); ++j) {
        if (chemProfs[j].in != chemProfs[j].out) {
            fireAvail = true;
        }
    }

    //Y1
    f[Y](i, 0) = 0;
    for (uint64_t j = 1; j < n; ++j) {
        //5
        auto toFind = [&] (double x) -> double {
            return std::pow(x, nu) / (f[U](i, j) * f[RHO](i, j));
        };
        double tmp = IntegralTrapeze(toFind, Xi_0, j * xi_h + Xi_0, xi_h);
        tmp *= (nu + 1);
        f[Y](i, j) = std::pow(tmp, 1.0 / (nu + 1));
    }
    //printIterData(Y);

    //U
    // for (int j = 0; j < 5; ++j) {
    //     f[U](i, 0) = f[U](i - 1, 0) + 8.0 * L / (H * H * H) * (alpha * F2(i, 0) * f[U](i - 1, 1) + (1.0 - alpha) / H * F2(i - 1, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0)) -
    //                  (f[P](i, 0) - f[P](i - 1, 0)) * (alpha / (f[RHO](i, 0) * f[U](i, 0)) + (1.0 - alpha) / (f[RHO](i - 1, 0) * f[U](i - 1, 0))));
    //     f[U](i, 0) /= (1.0 + 8.0 * alpha * L / (H * H * H) * F2(i, 0));
    // }
    //main
    for (int k = 0; k < 5; ++k) {
        M(0, 0) = 1.0 / L + alpha * std::pow(xi_h + Xi_0, nu) / (2 * H) * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H + comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0) / H)
                + alpha * 4.0 / 3.0 * comp({V}, {U}, i, 0) * yp(i, 0) / (2 * H) * (comp({MU, RHO, V}, {}, i, 1) * yp(i, 1) / H + comp({MU, RHO, V}, {}, i, 0) * yp(i, 0) / H);
        M(0, 1) = -alpha * std::pow(xi_h + Xi_0, nu) / (2 * H) * comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H
                - alpha * 4.0 / 3.0 * comp({V}, {U}, i, 0) * yp(i, 0) / (2 * H) * comp({MU, RHO, V}, {}, i, 1) * yp(i, 1) / H
                - alpha * 2.0 / 3.0 * comp({V}, {U}, i, 0) * yp(i, 0) / (2 * H) * comp({MU, RHO}, {}, i, 1) * yp(i, 1) * (f[V](i, 1) - f[V](i, 0)) / H;
        ans[0]  = f[U](i - 1, 0) / L
                - (alpha * comp({}, {RHO, U}, i, 0) + (1.0 - alpha) * comp({}, {RHO, U}, i - 1, 0)) * (f[P](i, 0) - f[P](i - 1, 0)) / L
                + alpha * comp({V}, {U}, i, 0) * yp(i, 0) * (f[P](i, 1) - f[P](i, 0)) / (2 * H) + (1.0 - alpha) * comp({V}, {U}, i - 1, 0) * yp(i - 1, 0) * (f[P](i - 1, 1) - f[P](i - 1, 0)) / (2 * H)
                + (1.0 - alpha) / std::pow(xi_h + Xi_0, nu) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 1) * yp2(i - 1, 1) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H - comp({MU, RHO, U}, {}, i - 1, 0) * yp2(i - 1, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H)
                + (1.0 - alpha) * 4.0 / 3.0 * comp({V}, {U}, i - 1, 0) * yp(i - 1, 0) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, 1) * yp(i - 1, 1) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H - comp({MU, RHO, V}, {}, i - 1, 0) * yp(i - 1, 0) * (f[U](i - 1, 0) - f[U](i - 1, 0)) / H)
                + alpha * yp(i, 0) / (2 * H) * (f[MU](i, 1) * (f[V](i, 1) - f[V](i - 1, 1)) / L - f[MU](i, 0) * (f[V](i, 0) - f[V](i - 1, 0)) / L) + (1.0 - alpha) * yp(i - 1, 0) / (2 * H) * (f[MU](i - 1, 1) * (f[V](i, 1) - f[V](i - 1, 1)) / L - f[MU](i - 1, 0) * (f[V](i, 0) - f[V](i - 1, 0)) / L)
                + alpha * yp(i, 0) / (2 * H) * (comp({MU, RHO, V}, {}, i, 1) * yp(i, 1) * (f[V](i, 1) - f[V](i, 0)) / H - comp({MU, RHO, V}, {}, i, 0) * yp(i, 0) * (f[V](i, 0) - f[V](i, 0)) / H) + (1.0 - alpha) * yp(i - 1, 0) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, 1) * yp(i - 1, 1) * (f[V](i - 1, 1) - f[V](i - 1, 0)) / H - comp({MU, RHO, V}, {}, i - 1, 0) * yp(i - 1, 0) * (f[V](i - 1, 0) - f[V](i - 1, 0)) / H)
                + (1.0 - alpha) * 2.0 / 3.0 * comp({V}, {U}, i - 1, 0) * yp(i - 1, 0) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 1) * yp(i - 1, 1) * (f[V](i - 1, 1) - f[V](i - 1, 0)) / H - comp({MU, RHO, U}, {}, i - 1, 0) * yp(i - 1, 0) * (f[V](i - 1, 0) - f[V](i - 1, 0)) / H);
        for (uint64_t j = 1; j < n - 1; ++j) {
            double xi = (j + 1) * xi_h + Xi_0;
            M(j, j - 1) = -alpha * std::pow(xi, nu) / (2 * H) * comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H
                        - alpha * 4.0 / 3.0 * comp({V}, {U}, i, j) * yp(i, j) / (2 * H) * comp({MU, RHO, V}, {}, i, j - 1) * yp(i, j - 1) / H
                        + alpha * 2.0 / 3.0 * comp({V}, {U}, i, j) * yp(i, j) / (2 * H) * comp({MU, RHO}, {}, i, j - 1) * yp(i, j - 1) * (f[V](i, j) - f[V](i, j - 1)) / H;
            M(j, j + 1) = -alpha * std::pow(xi, nu) / (2 * H) * comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H
                        - alpha * 4.0 / 3.0 * comp({V}, {U}, i, j) * yp(i, j) / (2 * H) * comp({MU, RHO, V}, {}, i, j + 1) * yp(i, j + 1) / H
                        - alpha * 2.0 / 3.0 * comp({V}, {U}, i, j) * yp(i, j) / (2 * H) * comp({MU, RHO}, {}, i, j + 1) * yp(i, j + 1) * (f[V](i, j + 1) - f[V](i, j)) / H;
            M(j, j) = 1.0 / L + alpha * std::pow(xi, nu) / (2 * H) * (comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H + comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H)
                    + alpha * 4.0 / 3.0 * comp({V}, {U}, i, j) * yp(i, j) / (2 * H) * (comp({MU, RHO, V}, {}, i, j + 1) * yp(i, j + 1) / H + comp({MU, RHO, V}, {}, i, j - 1) * yp(i, j - 1) / H);
            ans[j]  = f[U](i - 1, j) / L
                    - (alpha * comp({}, {RHO, U}, i, j) + (1.0 - alpha) * comp({}, {RHO, U}, i - 1, j)) * (f[P](i, j) - f[P](i - 1, j)) / L
                    + alpha * comp({V}, {U}, i, j) * yp(i, j) * (f[P](i, j + 1) - f[P](i, j - 1)) / (2 * H) + (1.0 - alpha) * comp({V}, {U}, i - 1, j) * yp(i - 1, j) * (f[P](i - 1, j + 1) - f[P](i - 1, j - 1)) / (2 * H)
                    + (1.0 - alpha) / std::pow(xi, nu) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp2(i - 1, j + 1) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) / H - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp2(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H)
                    + (1.0 - alpha) * 4.0 / 3.0 * comp({V}, {U}, i - 1, j) * yp(i - 1, j) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) / H - comp({MU, RHO, V}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H)
                    + alpha * yp(i, j) / (2 * H) * (f[MU](i, j + 1) * (f[V](i, j + 1) - f[V](i - 1, j + 1)) / L - f[MU](i, j - 1) * (f[V](i, j - 1) - f[V](i - 1, j - 1)) / L) + (1.0 - alpha) * yp(i - 1, j) / (2 * H) * (f[MU](i - 1, j + 1) * (f[V](i, j + 1) - f[V](i - 1, j + 1)) / L - f[MU](i - 1, j - 1) * (f[V](i, j - 1) - f[V](i - 1, j - 1)) / L)
                    + alpha * yp(i, j) / (2 * H) * (comp({MU, RHO, V}, {}, i, j + 1) * yp(i, j + 1) * (f[V](i, j + 1) - f[V](i, j)) / H - comp({MU, RHO, V}, {}, i, j - 1) * yp(i, j - 1) * (f[V](i, j) - f[V](i, j - 1)) / H) + (1.0 - alpha) * yp(i - 1, j) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * (f[V](i - 1, j + 1) - f[V](i - 1, j)) / H - comp({MU, RHO, V}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * (f[V](i - 1, j) - f[V](i - 1, j - 1)) / H)
                    + (1.0 - alpha) * 2.0 / 3.0 * comp({V}, {U}, i - 1, j) * yp(i - 1, j) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * (f[V](i - 1, j + 1) - f[V](i - 1, j)) / H - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * (f[V](i - 1, j) - f[V](i - 1, j - 1)) / H);
        }
        M(n - 1, n - 2) = -alpha * std::pow(xi_h * (n - 1) + Xi_0, nu) / (2 * H) * comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H
                        - alpha * 4.0 / 3.0 * comp({V}, {U}, i, n - 1) * yp(i, n - 1) / (2 * H) * comp({MU, RHO, V}, {}, i, n - 2) * yp(i, n - 2) / H
                        + alpha * 2.0 / 3.0 * comp({V}, {U}, i, n - 1) * yp(i, n - 1) / (2 * H) * comp({MU, RHO}, {}, i, n - 2) * yp(i, n - 2) * (f[V](i, n - 1) - f[V](i, n - 2)) / H;
        M(n - 1, n - 1) = 1.0 / L + alpha * std::pow(xi_h * (n - 1) + Xi_0, nu) / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp2(i, n - 1) / H + comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H)
                        + alpha * 4.0 / 3.0 * comp({V}, {U}, i, n - 1) * yp(i, n - 1) / (2 * H) * (comp({MU, RHO, V}, {}, i, n - 1) * yp(i, n - 1) / H + comp({MU, RHO, V}, {}, i, n - 2) * yp(i, n - 2) / H);
        ans[n - 1]  = f[U](i - 1, n - 1) / L
                    - (alpha * comp({}, {RHO, U}, i, n - 1) + (1.0 - alpha) * comp({}, {RHO, U}, i - 1, n - 1)) * (f[P](i, n - 1) - f[P](i - 1, n - 1)) / L
                    + alpha * comp({V}, {U}, i, n - 1) * yp(i, n - 1) * (f[P](i, n - 1) - f[P](i, n - 2)) / (2 * H) + (1.0 - alpha) * comp({V}, {U}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[P](i - 1, n - 1) - f[P](i - 1, n - 2)) / (2 * H)
                    + (1.0 - alpha) / std::pow(xi_h * (n - 1) + Xi_0, nu) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp2(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp2(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H)
                    + (1.0 - alpha) * 4.0 / 3.0 * comp({V}, {U}, i - 1, n - 1) * yp(i - 1, n - 1) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1)) / H - comp({MU, RHO, V}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H)
                    + alpha * yp(i, n - 1) / (2 * H) * (f[MU](i, n - 1) * (f[V](i, n - 1) - f[V](i - 1, n - 1)) / L - f[MU](i, n - 2) * (f[V](i, n - 2) - f[V](i - 1, n - 2)) / L) + (1.0 - alpha) * yp(i - 1, n - 1) / (2 * H) * (f[MU](i - 1, n - 1) * (f[V](i, n - 1) - f[V](i - 1, n - 1)) / L - f[MU](i - 1, n - 2) * (f[V](i, n - 2) - f[V](i - 1, n - 2)) / L)
                    + alpha * yp(i, n - 1) / (2 * H) * (comp({MU, RHO, V}, {}, i, n - 1) * yp(i, n - 1) * (f[V](i, n - 1) - f[V](i, n - 1)) / H - comp({MU, RHO, V}, {}, i, n - 2) * yp(i, n - 2) * (f[V](i, n - 1) - f[V](i, n - 2)) / H) + (1.0 - alpha) * yp(i - 1, n - 1) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 1)) / H - comp({MU, RHO, V}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 2)) / H)
                    + (1.0 - alpha) * 2.0 / 3.0 * comp({V}, {U}, i - 1, n - 1) * yp(i - 1, n - 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 1)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 2)) / H);
        if (alpha != 0.0) {
            ans = GaussSolveSLAE(M, ans, "U");
        }
        for (uint64_t j = 0; j < ans.size(); ++j) {
            f[U](i, j) = ans[j];
        }
    }
    f[U](i, n - 1) = f[U](i - 1, n - 1);
    //printIterData(U);

    //W
    // Mw(0, 0) = 1.0 / L + alpha * comp({V, W}, {U, Y}, i, 1) + yp(i, 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, 2) * yp(i, 2) + yp(i, 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, 0) * yp(i, 0);
    // Mw(0, 1) = -yp(i, 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, 2) * yp(i, 2) + yp(i, 1) * alpha / (2 * H) * comp({MU}, {Y}, i, 2) - f[MU](i, 1) * yp(i, 1) * alpha / H / f[Y](i, 2);
    // answ[0] = f[W](i - 1, 1) / L + (1.0 - alpha) * 
    //         (-comp({V, W}, {U, Y}, i - 1, 1)
    //             + yp(0, 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 2) * yp(i - 1, 2) * (f[W](i - 1, 2) - f[W](i - 1, 1)) / H    -   comp({MU, RHO, U}, {}, i - 1, 0) * yp(i - 1, 0) * (f[W](i - 1, 1) - f[W](i - 1, 0)) / H)
    //             - yp(0, 1) / (2 * H) * (comp({MU, W}, {Y}, i - 1, 2) - comp({MU, W}, {Y}, i - 1, 1))
    //             + 2 * f[MU](i - 1, 1) * yp(i - 1, 1) / (2 * H) * (comp({W}, {Y}, i - 1, 2) - comp({W}, {Y}, i - 1, 1))
    //         );
    // for (uint64_t j = 1; j < n - 2; ++j) {
    //     double xi = (j + 1) * xi_h + Xi_0;
    //     // Mw(j, j - 1) = -yp(i, j) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1) - yp(i, j) * alpha / (2 * H) * comp({MU}, {Y}, i, j - 1) + f[MU](i, j) * yp(i, j) * alpha / H / f[Y](i, j - 1);
    //     // Mw(j, j) = 1.0 / L + alpha * comp({V}, {U, Y}, i, j) + yp(i, j) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) + yp(i, j) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1);
    //     // Mw(j, j + 1) = -yp(i, j) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) + yp(i, j) * alpha / (2 * H) * comp({MU}, {Y}, i, j + 1) - f[MU](i, j) * yp(i, j) * alpha / H / f[Y](i, j + 1);
    //     Mw(j, j - 1) = -yp(i, j + 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j) * yp(i, j) - yp(i, j + 1) * alpha / (2 * H) * comp({MU}, {Y}, i, j) + f[MU](i, j + 1) * yp(i, j + 1) * alpha / H / f[Y](i, j);
    //     Mw(j, j) = 1.0 / L + alpha * comp({V}, {U, Y}, i, j + 1) + yp(i, j + 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j + 2) * yp(i, j + 2) + yp(i, j + 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j) * yp(i, j);
    //     Mw(j, j + 1) = -yp(i, j + 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, j + 2) * yp(i, j + 2) + yp(i, j + 1) * alpha / (2 * H) * comp({MU}, {Y}, i, j + 2) - f[MU](i, j + 1) * yp(i, j + 1) * alpha / H / f[Y](i, j + 2);
    //     answ[j] = f[W](i - 1, j + 1) / L + (1.0 - alpha) * 
    //             (-comp({V, W}, {U, Y}, i - 1, j + 1)
    //              + yp(i - 1, j + 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 2) * yp(i - 1, j + 2) * (f[W](i - 1, j + 2) - f[W](i - 1, j + 1)) / H    -   comp({MU, RHO, U}, {}, i - 1, j) * yp(i - 1, j) * (f[W](i - 1, j + 1) - f[W](i - 1, j)) / H)
    //              - yp(i - 1, j + 1) / (2 * H) * (comp({MU, W}, {Y}, i - 1, j + 2) - comp({MU, W}, {Y}, i - 1, j))
    //              + 2 * f[MU](i - 1, j + 1) * yp(i - 1, j + 1) / (2 * H) * (comp({W}, {Y}, i - 1, j + 2) - comp({W}, {Y}, i - 1, j))
    //             );
    // }
    // Mw(n - 2, n - 3) = -yp(i, n - 2) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, n - 3) * yp(i, n - 3) - yp(i, n - 2) * alpha / (2 * H) * comp({MU}, {Y}, i, n - 3) + f[MU](i, n - 2) * yp(i, n - 2) * alpha / H / f[Y](i, n - 3);
    // Mw(n - 2, n - 2) = 1.0 / L + alpha * comp({V}, {U, Y}, i, n - 1) + yp(i, n - 1) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 1) + yp(i, n - 2) * alpha / (2 * H * H) * comp({MU, RHO, U}, {}, i, n - 3) * yp(i, n - 3);
    // answ[n - 2] = f[W](i - 1, n - 2) / L + (1.0 - alpha) * 
    //             (-comp({V, W}, {U, Y}, i - 1, n - 2)
    //              + yp(i - 1, n - 2) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[W](i - 1, n - 1) - f[W](i - 1, n - 2)) / H    -   comp({MU, RHO, U}, {}, i - 1, n - 3) * yp(i - 1, n - 3) * (f[W](i - 1, n - 2) - f[W](i - 1, n - 3)) / H)
    //              - yp(i - 1, n - 2) / (2 * H) * (comp({MU, W}, {Y}, i - 1, n - 1) - comp({MU, W}, {Y}, i - 1, n - 3))
    //              + 2 * f[MU](i - 1, n - 2) * yp(i - 1, n - 2) / (2 * H) * (comp({W}, {Y}, i - 1, n - 1) - comp({W}, {Y}, i - 1, n - 3))
    //             );
    // if (alpha != 0.0) {
    //     answ = GaussSolveSLAE(Mw, answ, "W");
    // }
    // for (uint64_t j = 0; j < answ.size(); ++j) {
    //     f[W](i, j + 1) = answ[j];
    // }

    //V
    //1
    // Mw(0, 1) = 4.0 / 3.0 * alpha * L / (H * H) * yp(i, 1) * F3(i, 1) + nu * alpha * L / (H * f[Y](i, 1)) * yp(i, 1) * (f[MU](i, 1) - f[MU](i, 2));
    // Mw(0, 0) = -1.0 - Mw(0, 1) - 4.0 / 3.0 * nu * alpha * comp({MU}, {Y, Y, RHO, U}, i, 1);
    // answ[0]  =  -f[V](i - 1, 1)
    //             + alpha * yp(i, 1) * (f[P](i, 2) - f[P](i, 0)) / (2 * H) + (1.0 - alpha) * yp(i - 1, 1) * (f[P](i - 1, 2) - f[P](i - 1, 0)) / (2 * H) //0 & 1
    //             - (1.0 - alpha) / (H * H) * L * yp(i - 1, 1) * (F3(i - 1, 1) * (f[V](i - 1, 2) - f[V](i - 1, 1)) - F3(i - 1, 0) * (f[V](i - 1, 1) - f[V](i - 1, 0)))
    //             - (1.0 - alpha) * nu * L * comp({MU}, {Y}, i - 1, 1) / (H * H) * yp(i - 1, 1) * (f[V](i - 1, 2) - f[V](i - 1, 0))
    //             + (1.0 - alpha) / (3 * H) * nu * L / f[Y](i - 1, 1) * yp(i - 1, 1) * (comp({MU, V}, {}, i - 1, 2) - comp({MU, V}, {}, i - 1, 0))
    //             + 4.0 / 3.0 * (1.0 - alpha) * nu * comp({V, MU}, {Y, Y, RHO, U}, i - 1, 1)
    //             ;
    // for (uint64_t j = 2; j < n - 1; ++j) {
    //     double xi = (j + 1) * xi_h + Xi_0;
    //     Mw(j - 1, j - 2) = 4.0 / 3.0 * alpha * L / (H * H) * yp(i, j) * F3(i, j - 1) - nu * alpha * L / (H * f[Y](i, j)) * yp(i, j) * (f[MU](i, j) - f[MU](i, j - 1));
    //     Mw(j - 1, j) = 4.0 / 3.0 * alpha * L / (H * H) * yp(i, j) * F3(i, j) + nu * alpha * L / (H * f[Y](i, j)) * yp(i, j) * (f[MU](i, j) - f[MU](i, j + 1));
    //     Mw(j - 1, j - 1) = -1.0 - Mw(j - 1, j - 2) - Mw(j - 1, j) - 4.0 / 3.0 * nu * alpha * comp({MU}, {Y, Y, RHO, U}, i, j);
    //     answ[j - 1] = -f[V](i - 1, j)
    //                 + alpha * yp(i, j) * (f[P](i, j + 1) - f[P](i, j - 1)) / (2 * H) + (1.0 - alpha) * yp(i - 1, j) * (f[P](i - 1, j + 1) - f[P](i - 1, j - 1)) / (2 * H) //0 & 1
    //                 - (1.0 - alpha) / (H * H) * L * yp(i - 1, j) * (F3(i - 1, j) * (f[V](i - 1, j + 1) - f[V](i - 1, j)) - F3(i - 1, j - 1) * (f[V](i - 1, j) - f[V](i - 1, j - 1)))
    //                 - (1.0 - alpha) * nu * L * comp({MU}, {Y}, i - 1, j) / (H * H) * yp(i - 1, j) * (f[V](i - 1, j + 1) - f[V](i - 1, j - 1))
    //                 + (1.0 - alpha) / (3 * H) * nu * L / f[Y](i - 1, j) * yp(i - 1, j) * (comp({MU, V}, {}, i - 1, j + 1) - comp({MU, V}, {}, i - 1, j - 1))
    //                 + 4.0 / 3.0 * (1.0 - alpha) * nu * comp({V, MU}, {Y, Y, RHO, U}, i - 1, j)
    //                 ;
    // }
    // Mw(n - 2, n - 3) = 4.0 / 3.0 * alpha * L / (H * H) * yp(i, n - 1) * F3(i, n - 2) - nu * alpha * L / (H * f[Y](i, n - 1)) * yp(i, n - 1) * (f[MU](i, n - 1) - f[MU](i, n - 2));
    // Mw(n - 2, n - 2) = -1.0 - Mw(n - 2, n - 3) - 4.0 / 3.0 * nu * alpha * comp({MU}, {Y, Y, RHO, U}, i, n - 1);
    // answ[n - 2] = -f[V](i - 1, n - 1)
    //             + alpha * yp(i, n - 1) * (f[P](i, n - 1) - f[P](i, n - 2)) / (2 * H) + (1.0 - alpha) * yp(i - 1, n - 1) * (f[P](i - 1, n - 1) - f[P](i - 1, n - 2)) / (2 * H) //0 & 1
    //             - (1.0 - alpha) / (H * H) * L * yp(i - 1, n - 1) * (F3(i - 1, n - 1) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 1)) - F3(i - 1, n - 2) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 2)))
    //             - (1.0 - alpha) * nu * L * comp({MU}, {Y}, i - 1, n - 1) / (H * H) * yp(i - 1, n - 1) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 2))
    //             + (1.0 - alpha) / (3 * H) * nu * L / f[Y](i - 1, n - 1) * yp(i - 1, n - 1) * (comp({MU, V}, {}, i - 1, n - 1) - comp({MU, V}, {}, i - 1, n - 2))
    //             + 4.0 / 3.0 * (1.0 - alpha) * nu * comp({V, MU}, {Y, Y, RHO, U}, i - 1, n - 1)
    //             ;
    //2
    // for (uint64_t j = 2; j < n - 1; ++j) {
    //     double xi = (j + 1) * xi_h + Xi_0;
    //     Mw(j - 1, j - 2) = 1;
    //     Mw(j - 1, j) = 1;
    //     Mw(j - 1, j - 1) = -1.0 / L - Mw(j - 1, j - 2) - Mw(j - 1, j);
    //     answ[j - 1] = 1;
    // }
    //main
    double midRhoU = alpha * comp({RHO, U}, {}, i, 1) + (1.0 - alpha) * comp({RHO, U}, {}, i - 1, 1);
    Mw(0, 0) = 1.0 / L + alpha * (4.0 / 3.0) * yp(i, 1) / (2 * H) * (comp({MU, RHO, U}, {}, i, 2) * yp(i, 2) / H + comp({MU, RHO, U}, {}, i, 0) * yp(i, 0) / H)
               + alpha * f[U](i, 1) * yp(i, 1) / (2 * H) * (comp({MU, RHO, U}, {}, i, 2) * yp(i, 2) * (f[U](i, 2) - f[U](i, 1)) / H - comp({MU, RHO, U}, {}, i, 0) * yp(i, 0) * (f[U](i, 1) - f[U](i, 0)) / H)
               ;
               //Mw(0, 0) /= midRhoU;
    Mw(0, 1) = - alpha * (4.0 / 3.0) * yp(i, 1) / (2 * H) * comp({MU, RHO, U}, {}, i, 2) * yp(i, 2) / H
               + alpha * 2.0 / 3.0 * nu * yp(i, 1) / (2 * H) * comp({MU}, {Y}, i, 2)
               - alpha * 2 * nu * f[MU](i, 1) * yp(i, 1) / (2 * H) / f[Y](i, 2)
               - alpha * (2.0 / 3.0) * yp(i, 1) / (2 * H) * comp({MU, RHO}, {}, i, 2) * yp(i, 2) * (f[U](i, 2) - f[U](i, 1)) / H
               ;
               //Mw(0, 1) /= midRhoU;
    answ[0] = 1.0 * f[V](i - 1, 1) / L - alpha * yp(i, 1) * (f[P](i, 2) - f[P](i, 0)) / (2 * H) - (1.0 - alpha) * yp(i - 1, 1) * (f[P](i - 1, 2) - f[P](i - 1, 0)) / (2 * H)
              + (1.0 - alpha) * (4.0 / 3.0) * yp(i - 1, 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 2) * yp(i - 1, 2) * (f[V](i - 1, 2) - f[V](i - 1, 1)) / H - comp({MU, RHO, U}, {}, i - 1, 0) * yp(i - 1, 0) * (f[V](i - 1, 1) - f[V](i - 1, 0)) / H)
              - (1.0 - alpha) * (2.0 / 3.0) * nu * yp(i - 1, 1) / H * (comp({MU, V}, {Y}, i - 1, 2) - comp({MU, V}, {Y}, i - 1, 1))
              + (1.0 - alpha) * 2 * nu * f[MU](i - 1, 1) * yp(i - 1, 1) / H * (comp({V}, {Y}, i - 1, 2) - comp({V}, {Y}, i - 1, 1))
              + (1.0 - alpha) * (2.0 / 3.0) * yp(i - 1, 1) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, 2) * yp(i - 1, 2) * (f[U](i - 1, 2) - f[U](i - 1, 1)) / H - comp({MU, RHO, V}, {}, i - 1, 0) * yp(i - 1, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H)
              - (1.0 - alpha) * comp({V}, {U}, i - 1, 1) * yp(i - 1, 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 2) * yp(i - 1, 2) * (f[U](i - 1, 2) - f[U](i - 1, 1)) / H - comp({MU, RHO, U}, {}, i - 1, 0) * yp(i - 1, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H)
              - 2.0 / 3.0 * yp(i - 1, 1) / (2 * H) * (f[MU](i - 1, 2) * (f[U](i, 2) - f[U](i - 1, 2)) / L - f[MU](i - 1, 0) * (f[U](i, 0) - f[U](i - 1, 0)) / L) //3
              //+ comp({}, {RHO, U}, i - 1, 1) / L * (comp({MU, RHO, U}, {}, i, 1) * yp(i, 1) * (f[U](i, 2) - f[U](i, 1)) / H - comp({MU, RHO, U}, {}, i - 1, 1) * yp(i - 1, 1) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H) //5
              ;
              //answ[0] /= midRhoU;
    for (uint64_t j = 2; j < n - 1; ++j) {
        double xi = (j + 1) * xi_h + Xi_0;
        midRhoU = alpha * comp({RHO, U}, {}, i, j) + (1.0 - alpha) * comp({RHO, U}, {}, i - 1, j);
        Mw(j - 1, j - 2) = - alpha * (4.0 / 3.0) * yp(i, j) / (2 * H) * comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1) / H
                           - alpha * 2.0 / 3.0 * nu * yp(i, j) / (2 * H) * comp({MU}, {Y}, i, j - 1) //7
                           + alpha * 2 * nu * f[MU](i, j) * yp(i, j) / (2 * H) / f[Y](i, j - 1)
                           + alpha * (2.0 / 3.0) * yp(i, j) / (2 * H) * comp({MU, RHO}, {}, i, j - 1) * yp(i, j - 1) * (f[U](i, j) - f[U](i, j - 1)) / H //4
                           ;
                            //Mw(j - 1, j - 2) /= midRhoU;
        Mw(j - 1, j - 1) = 1.0 / L + alpha * (4.0 / 3.0) * yp(i, j) / (2 * H) * (comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) / H + comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1) / H)
                           + alpha * f[U](i, j) * yp(i, j) / (2 * H) * (comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) * (f[U](i, j + 1) - f[U](i, j)) / H - comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1) * (f[U](i, j) - f[U](i, j - 1)) / H) //6
                           ;
                            //Mw(j - 1, j - 1) /= midRhoU;
        Mw(j - 1, j) = - alpha * (4.0 / 3.0) * yp(i, j) / (2 * H) * comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) / H
                       + alpha * 2.0 / 3.0 * nu * yp(i, j) / (2 * H) * comp({MU}, {Y}, i, j + 1) //7
                       - alpha * 2 * nu * f[MU](i, j) * yp(i, j) / (2 * H) / f[Y](i, j + 1)
                       - alpha * (2.0 / 3.0) * yp(i, j) / (2 * H) * comp({MU, RHO}, {}, i, j + 1) * yp(i, j + 1) * (f[U](i, j + 1) - f[U](i, j)) / H //4
                       ;
                        //Mw(j - 1, j) /= midRhoU;
        answ[j - 1] = 1.0 * f[V](i - 1, j) / L
                      - alpha * yp(i, j) * (f[P](i, j + 1) - f[P](i, j - 1)) / (2 * H) - (1.0 - alpha) * yp(i - 1, j) * (f[P](i - 1, j + 1) - f[P](i - 1, j - 1)) / (2 * H) //0 & 1
                      + (1.0 - alpha) * (4.0 / 3.0) * yp(i - 1, j) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * (f[V](i - 1, j + 1) - f[V](i - 1, j)) / H - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * (f[V](i - 1, j) - f[V](i - 1, j - 1)) / H) //2
                      - (1.0 - alpha) * (2.0 / 3.0) * nu * yp(i - 1, j) / (2 * H) * (comp({MU, V}, {Y}, i - 1, j + 1) - comp({MU, V}, {Y}, i - 1, j - 1)) //7
                      + (1.0 - alpha) * 2 * nu * f[MU](i - 1, j) * yp(i - 1, j) / (2 * H) * (comp({V}, {Y}, i - 1, j + 1) - comp({V}, {Y}, i - 1, j - 1)) //8
                      + (1.0 - alpha) * (2.0 / 3.0) * yp(i - 1, j) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) / H - comp({MU, RHO, V}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H) //4
                      - (1.0 - alpha) * comp({V}, {U}, i - 1, j) * yp(i - 1, j) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) / H - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H) //6
                      - 2.0 / 3.0 * yp(i - 1, j) / (2 * H) * (f[MU](i - 1, j + 1) * (f[U](i, j + 1) - f[U](i - 1, j + 1)) / L - f[MU](i - 1, j - 1) * (f[U](i, j - 1) - f[U](i - 1, j - 1)) / L) //3
                      //+ comp({}, {RHO, U}, i - 1, j) / L * (comp({MU, RHO, U}, {}, i, j) * yp(i, j) * (f[U](i, j + 1) - f[U](i, j)) / H - comp({MU, RHO, U}, {}, i - 1, j) * yp(i - 1, j) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H) //5
                      ;
                        //answ[j - 1] /= midRhoU;
    }
    midRhoU = alpha * comp({RHO, U}, {}, i, n - 1) + (1.0 - alpha) * comp({RHO, U}, {}, i - 1, n - 1);
    Mw(n - 2, n - 3) =  - alpha * (4.0 / 3.0) * yp(i, n - 1) / (2 * H) * comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 2) / H
                        - alpha * 2.0 / 3.0 * nu * yp(i, n - 1) / (2 * H) * comp({MU}, {Y}, i, n - 2) //7
                        + alpha * 2 * nu * f[MU](i, n - 1) * yp(i, n - 1) / (2 * H) / f[Y](i, n - 2)
                        + alpha * (2.0 / 3.0) * yp(i, n - 1) / (2 * H) * comp({MU, RHO}, {}, i, n - 2) * yp(i, n - 2) * (f[U](i, n - 1) - f[U](i, n - 2)) / H //4
                        ;
                        //Mw(n - 2, n - 3) /= midRhoU;
    Mw(n - 2, n - 2) = 1.0 / L + alpha * (4.0 / 3.0) * yp(i, n - 1) / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) / H + comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 2) / H)
                        + alpha * f[U](i, n - 1) * yp(i, n - 1) / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) * (f[U](i, n - 1) - f[U](i, n - 1)) / H - comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 2) * (f[U](i, n - 1) - f[U](i, n - 2)) / H) //6
                        ;
                        //Mw(n - 2, n - 2) /= midRhoU;
    answ[n - 2] = 1.0 * f[V](i - 1, n - 1) / L - alpha * yp(i, n - 1) * (f[P](i, n - 1) - f[P](i, n - 2)) / (2 * H) - (1.0 - alpha) * yp(i - 1, n - 1) * (f[P](i - 1, n - 1) - f[P](i - 1, n - 2)) / (2 * H) //0 & 1
                    + (1.0 - alpha) * (4.0 / 3.0) * yp(i - 1, n - 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 1)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 2)) / H) //2
                    - (1.0 - alpha) * (2.0 / 3.0) * nu * yp(i - 1, n - 1) / (2 * H) * (comp({MU, V}, {Y}, i - 1, n - 1) - comp({MU, V}, {Y}, i - 1, n - 2)) //7
                    + (1.0 - alpha) * 2 * nu * f[MU](i - 1, n - 1) * yp(i - 1, n - 1) / (2 * H) * (comp({V}, {Y}, i - 1, n - 1) - comp({V}, {Y}, i - 1, n - 2)) //8
                    + (1.0 - alpha) * (2.0 / 3.0) * yp(i - 1, n - 1) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1)) / H - comp({MU, RHO, V}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H) //4
                    - (1.0 - alpha) * comp({V}, {U}, i - 1, n - 1) * yp(i - 1, n - 1) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H) //6
                    - 2.0 / 3.0 * yp(i - 1, n - 1) / (2 * H) * (f[MU](i - 1, n - 1) * (f[U](i, n - 1) - f[U](i - 1, n - 1)) / L - f[MU](i - 1, n - 2) * (f[U](i, n - 2) - f[U](i - 1, n - 2)) / L) //3
                    //+ comp({}, {RHO, U}, i - 1, n - 1) / L * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) * (f[U](i, n - 1) - f[U](i, n - 1)) / H - comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H) //5
                    ;
                    //answ[n - 2] /= midRhoU;
    // Mw(n - 2, n - 3) = - alpha * (4.0 / 3.0) * yp(i, n - 2) / (2 * H) * comp({MU, RHO, U}, {}, i, n - 3) * yp(i, n - 3) / H
    //                    - alpha * 2.0 / 3.0 * nu * yp(i, n - 2) / (2 * H) * comp({MU}, {Y}, i, n - 3)
    //                    + alpha * 2 * nu * f[MU](i, n - 2) * yp(i, n - 2) / (2 * H) / f[Y](i, n - 3)
    //                    + alpha * (2.0 / 3.0) * yp(i, n - 2) / (2 * H) * comp({MU, RHO}, {}, i, n - 3) * yp(i, n - 3) * (f[U](i, n - 2) - f[U](i, n - 3)) / H
    //                    ;
    // Mw(n - 2, n - 2) = 1.0 / L + alpha * (4.0 / 3.0) * yp(i, n - 2) / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) / H + comp({MU, RHO, U}, {}, i, n - 3) * yp(i, n - 3) / H)
    //                    + alpha * f[U](i, n - 2) * yp(i, n - 2) / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) * (f[U](i, n - 1) - f[U](i, n - 2)) / H - comp({MU, RHO, U}, {}, i, n - 3) * yp(i, n - 3) * (f[U](i, n - 2) - f[U](i, n - 3)) / H)
    //                    ;
    // answ[n - 2] = f[V](i - 1, n - 2) / L - alpha * yp(i, n - 2) * (f[P](i, n - 1) - f[P](i, n - 3)) / (2 * H) - (1.0 - alpha) * yp(i - 1, n - 2) * (f[P](i - 1, n - 1) - f[P](i - 1, n - 3)) / (2 * H)
    //               + (1.0 - alpha) * (4.0 / 3.0) * yp(i - 1, n - 2) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[V](i - 1, n - 1) - f[V](i - 1, n - 2)) / H - comp({MU, RHO, U}, {}, i - 1, n - 3) * yp(i - 1, n - 3) * (f[V](i - 1, n - 2) - f[V](i - 1, n - 3)) / H)
    //               - (1.0 - alpha) * (2.0 / 3.0) * nu * yp(i - 1, n - 2) / (2 * H) * (comp({MU, V}, {Y}, i - 1, n - 1) - comp({MU, V}, {Y}, i - 1, n - 3))
    //               + (1.0 - alpha) * 2 * nu * f[MU](i - 1, n - 2) * yp(i - 1, n - 2) / (2 * H) * (comp({V}, {Y}, i - 1, n - 1) - comp({V}, {Y}, i - 1, n - 3))
    //               + (1.0 - alpha) * (2.0 / 3.0) * yp(i - 1, n - 2) / (2 * H) * (comp({MU, RHO, V}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H - comp({MU, RHO, V}, {}, i - 1, n - 3) * yp(i - 1, n - 3) * (f[U](i - 1, n - 2) - f[U](i - 1, n - 3)) / H)
    //               - (1.0 - alpha) * comp({V}, {U}, i - 1, n - 2) * yp(i - 1, n - 2) / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H - comp({MU, RHO, U}, {}, i - 1, n - 3) * yp(i - 1, n - 3) * (f[U](i - 1, n - 2) - f[U](i - 1, n - 3)) / H)
    //               - 2.0 / 3.0 * yp(i - 1, n - 2) / (2 * H) * (f[MU](i - 1, n - 1) * (f[U](i, n - 1) - f[U](i - 1, n - 1)) / L - f[MU](i - 1, n - 3) * (f[U](i, n - 3) - f[U](i - 1, n - 3)) / L) //3
    //               + comp({}, {RHO, U}, i - 1, n - 2) / L * (comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 2) * (f[U](i, n - 1) - f[U](i, n - 2)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * (f[U](i - 1, n - 2) - f[U](i - 1, n - 3)) / H) //5
    //               ;
    if (alpha != 0.0) {
        answ = GaussSolveSLAE(Mw, answ, "V");
    }
    for (uint64_t j = 0; j < answ.size(); ++j) {
        //f[V](i, j + 1) = std::abs(answ[j]);
        f[V](i, j + 1) = answ[j];
    }
    //printIterData(V);

    //J
    //1
    // double sum = 0;
    // for (int j = 0; j < 5; ++j) {
    //     f[J](i, 0) = f[J](i - 1, 0) + 8.0 * L / (H * H * H) * ((alpha / Pr) * F4(i, 0) * f[J](i - 1, 1) + (1.0 - alpha) / Pr * F4(i - 1, 0) * (f[J](i - 1, 1) - f[J](i - 1, 0)) +
    //                  alpha * (1.0 - 1.0 / Pr) * F4(i, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0)) +
    //                  (1.0 - alpha) * (1.0 - 1.0 / Pr) * F4(i - 1, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0))) -
    //                  (alpha * f[RHO](i, 0) * f[Q](i, 0) + (1.0 - alpha) * f[RHO](i - 1, 0) * f[Q](i - 1, 0));
    //     f[J](i, 0) /= (1.0 + 8.0 * alpha * L / (H * H * H * Pr) * F2(i, 0));
    // }
    // M(0, 1) = alpha * L / (H * H) * 1.0 / std::pow(Xi_0 + xi_h, nu) * F2(i, 0) / Pr;
    // M(0, 0) = -1.0 - M(0, 1);
    // ans[0] = -1.0 * f[J](i - 1, 0) - alpha / std::pow(Xi_0 + xi_h, nu) * L / (H * H) * (1.0 - 1.0 / Pr) * (F4(i, 0) * (f[U](i, 1) - f[U](i, 0)) - F4(i, 0) * (f[U](i, 0) - f[U](i, 0))) +
    //         1.0 / std::pow(Xi_0 + xi_h, nu) * (1.0 / Le - 1.0) * (alpha * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i, 0) * std::pow(f[Y](i, 0), 2 * nu - 1) / Pr);
    // double tmp = (alpha * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i, 0) * std::pow(f[Y](i, 0), 2 * nu - 1) / Pr);
    // sum = 0;
    // for (uint64_t k = 0; k < ental.size(); ++k) {
    //     sum += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) * (chemProfs[k].in - chemProfs[k].out);
    // }
    // ans[0] += tmp * sum * 1.0 / std::pow(Xi_0 + xi_h, nu) * (1.0 / Le - 1.0);
    // tmp = (1.0 - alpha) * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i - 1, 0) * std::pow(f[Y](i - 1, 0), 2 * nu - 1) / Pr;
    // sum = 0;
    // for (uint64_t k = 0; k < ental.size(); ++k) {
    //     sum += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) * (chemProfs[k].in - chemProfs[k].out);
    // }
    // ans[0] += tmp * sum * 1.0 / std::pow(Xi_0 + xi_h, nu) * (1.0 / Le - 1.0);
    // ans[0] += -1.0 / std::pow(Xi_0 + xi_h, nu) * (1.0 - alpha) * L / (H * H) * (1.0 - 1.0 / Pr) * (F4(i - 1, 0) * (f[U](i - 1, 1) - f[U](i - 1, 0)) - F4(i - 1, 0) * (f[U](i - 1, 0) - f[U](i - 1, 0))) -
    //           (1.0 - alpha) * L / (H * H) * (1.0 / std::pow(Xi_0 + xi_h, nu)) / Pr * (F2(i - 1, 0) * (f[J](i - 1, 1) - f[J](i - 1, 0)) - F2(i - 1, 0) * (f[J](i - 1, 0) - f[J](i - 1, 0)));
    // for (uint64_t j = 1; j < n - 1; ++j) {
    //     double xi = (j + 1) * xi_h + Xi_0;

    //     M(j, j - 1) = alpha * L / (H * H) * 1.0 / std::pow(xi, nu) * F2(i, j - 1) / Pr;
    //     M(j, j + 1) = alpha * L / (H * H) * 1.0 / std::pow(xi, nu) * F2(i, j) / Pr;
    //     M(j, j) = -1.0 - M(j, j - 1) - M(j, j + 1);
    //     ans[j] = -f[J](i - 1, j) - alpha / std::pow(xi, nu) * L / (H * H) * (1.0 - 1.0 / Pr) * (F4(i, j) * (f[U](i, j + 1) - f[U](i, j)) - F4(i, j - 1) * (f[U](i, j) - f[U](i, j - 1)));
    //     tmp = (alpha * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i, j) * std::pow(f[Y](i, j), 2 * nu - 1) / Pr);
    //     sum = 0;
    //     for (uint64_t k = 0; k < ental.size(); ++k) {
    //         sum += ental[k](f[T](i, j) * valProfs[T].in) * (f[C1 + k](i, j + 1) - f[C1 + k](i, j - 1)) * (chemProfs[k].in - chemProfs[k].out);
    //     }
    //     ans[j] += tmp * sum * 1.0 / std::pow(xi, nu) * (1.0 / Le - 1.0);
    //     tmp = (1.0 - alpha) * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i - 1, j) * std::pow(f[Y](i - 1, j), 2 * nu - 1) / Pr;
    //     sum = 0;
    //     for (uint64_t k = 0; k < ental.size(); ++k) {
    //         sum += ental[k](f[T](i - 1, j) * valProfs[T].in) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j - 1)) * (chemProfs[k].in - chemProfs[k].out);
    //     }
    //     ans[j] += tmp * sum * 1.0 / std::pow(xi, nu) * (1.0 / Le - 1.0);
    //     ans[j] += -1.0 / std::pow(xi, nu) * (1.0 - alpha) * L / (H * H) * (1.0 - 1.0 / Pr) * (F4(i - 1, j) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) - F4(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1))) -
    //               (1.0 - alpha) * L / (H * H) * (1.0 / std::pow(xi, nu)) / Pr * (F2(i - 1, j) * (f[J](i - 1, j + 1) - f[J](i - 1, j)) - F2(i - 1, j - 1) * (f[J](i - 1, j) - f[J](i - 1, j - 1)));
    // }
    // M(n - 1, n - 2) = alpha * L / (H * H) * 1.0 / std::pow((Xi_0 + n * xi_h), nu) * F2(i, n - 2) / Pr;
    // M(n - 1, n - 1) = -1.0 - M(n - 1, n - 2);
    // ans[n - 1] = -1.0 * f[J](i - 1, n - 1) - alpha / std::pow((Xi_0 + n * xi_h), nu) * L / (H * H) * (1.0 - 1.0 / Pr) * (F4(i, n - 1) * (f[U](i, n - 1) - f[U](i, n - 1)) - F4(i, n - 2) * (f[U](i, n - 1) - f[U](i, n - 2))) +
    //         1.0 / std::pow((Xi_0 + n * xi_h), nu) * (1.0 / Le - 1.0) * (alpha * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i, n - 1) * std::pow(f[Y](i, n - 1), 2 * nu - 1) / Pr);
    // tmp = (alpha * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i, n - 1) * std::pow(f[Y](i, n - 1), 2 * nu - 1) / Pr);
    // sum = 0;
    // for (uint64_t k = 0; k < ental.size(); ++k) {
    //     sum += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) * (chemProfs[k].in - chemProfs[k].out);
    // }
    // ans[n - 1] += tmp * sum * 1.0 / std::pow(Xi_0 + n * xi_h, nu) * (1.0 / Le - 1.0);
    // tmp = (1.0 - alpha) * L / (2 * H) * comp({MU, RHO, U, Y}, {}, i - 1, n - 1) * std::pow(f[Y](i - 1, n - 1), 2 * nu - 1) / Pr;
    // sum = 0;
    // for (uint64_t k = 0; k < ental.size(); ++k) {
    //     sum += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) * (chemProfs[k].in - chemProfs[k].out);
    // }
    // ans[n - 1] += tmp * sum * 1.0 / std::pow(Xi_0 + n * xi_h, nu) * (1.0 / Le - 1.0);
    // ans[n - 1] += -1.0 / std::pow((Xi_0 + n * xi_h), nu) * (1.0 - alpha) * L / (H * H) * (1.0 - 1.0 / Pr) * (F4(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1)) - F4(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1))) -
    //               (1.0 - alpha) * L / (H * H) * (1.0 / std::pow((Xi_0 + n * xi_h), nu)) / Pr * (F2(i - 1, n - 1) * (f[J](i - 1, n - 1) - f[J](i - 1, n - 1)) - F2(i - 1, n - 1) * (f[J](i - 1, n - 1) - f[J](i - 1, n - 1)));
    // if (alpha != 0.0) {
    //     ans = GaussSolveSLAE(M, ans, "J");
    // }
    // for (uint64_t j = 0; j < ans.size(); ++j) {
    //     f[J](i, j) = ans[j];
    // }
    //2
    double sumLeftPrev = 0.0, sumRightPrev = 0.0, sumLeftNext = 0.0, sumRightNext = 0.0;
    M(0, 0) = 1.0 / L
            + alpha / std::pow(Xi_0 + xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H + comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0) / H);
    M(0, 1) = -alpha / std::pow(Xi_0 + xi_h, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H;
    sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    for (uint64_t k = 0; k < ental.size(); ++k) {
        // sumLeftPrev  += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumRightPrev += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumLeftNext  += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumRightNext += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 0) - f[C1 + k](i, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumLeftPrev  += ental[k](f[T](i - 1, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) / H;
        // sumRightPrev += ental[k](f[T](i - 1, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) / H;
        // sumLeftNext  += ental[k](f[T](i, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) / H;
        // sumRightNext += ental[k](f[T](i, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, 0) - f[C1 + k](i, 0)) / H;
        sumLeftPrev  += ental[k](f[T](i - 1, 1) * valProfs[T].in) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) * chemProfs[k].in / H;
        sumRightPrev += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) * chemProfs[k].in / H;
        sumLeftNext  += ental[k](f[T](i, 1) * valProfs[T].in) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) * chemProfs[k].in / H;
        sumRightNext += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 0) - f[C1 + k](i, 0)) * chemProfs[k].in / H;
    }
    ans[0]  = f[J](i - 1, 0) / L
            + (1.0 - alpha) / std::pow(Xi_0 + xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 1) * yp2(i - 1, 1) * (f[J](i - 1, 1) - f[J](i - 1, 0)) / H - comp({MU, RHO, U}, {}, i - 1, 0) * yp2(i - 1, 0) * (f[J](i - 1, 0) - f[J](i - 1, 0)) / H)
            + (1.0 - alpha) / std::pow(Xi_0 + xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i - 1, 1) * yp2(i - 1, 1) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H - comp({MU, RHO, U, U}, {}, i - 1, 0) * yp2(i - 1, 0) * (f[U](i - 1, 0) - f[U](i - 1, 0)) / H)
            + alpha / std::pow(Xi_0 + xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i, 1) * yp2(i, 1) * (f[U](i, 1) - f[U](i, 0)) / H - comp({MU, RHO, U, U}, {}, i, 0) * yp2(i, 0) * (f[U](i, 0) - f[U](i, 0)) / H)
            - (1.0 - alpha) * (1.0 - 1.0 / Le) / std::pow(Xi_0 + xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i - 1, 1) * yp(i - 1, 1) * sumLeftPrev - comp({MU, RHO, U}, {}, i - 1, 0) * yp(i - 1, 0) * sumRightPrev)
            - alpha * (1.0 - 1.0 / Le) / std::pow(Xi_0 + xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i, 1) * yp(i, 1) * sumLeftNext - comp({MU, RHO, U}, {}, i, 0) * yp(i, 0) * sumRightNext)
            ;
    for (uint64_t j = 1; j < n - 1; ++j) {
        double xi = (j + 1) * xi_h + Xi_0;
        M(j, j - 1) = -alpha / std::pow(xi, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H;
        M(j, j + 1) = -alpha / std::pow(xi, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H;
        M(j, j) = 1.0 / L
                + alpha / std::pow(xi, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H + comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H);
        sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
        for (uint64_t k = 0; k < ental.size(); ++k) {
            //1
            //sumLeftPrev  += ental[k](f[T](i - 1, j) * valProfs[T].in) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) * (chemProfs[k].in - chemProfs[k].out) / H;
            //sumRightPrev += ental[k](f[T](i - 1, j) * valProfs[T].in) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
            //sumLeftNext  += ental[k](f[T](i, j) * valProfs[T].in) * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) * (chemProfs[k].in - chemProfs[k].out) / H;
            //sumRightNext += ental[k](f[T](i, j) * valProfs[T].in) * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
            //2 (main)
            // sumLeftPrev  += ental[k](f[T](i - 1, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) / H;
            // sumRightPrev += ental[k](f[T](i - 1, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) / H;
            // sumLeftNext  += ental[k](f[T](i, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) / H;
            // sumRightNext += ental[k](f[T](i, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) / H;
            //3
            sumLeftPrev  += ental[k](f[T](i - 1, j + 1) * valProfs[T].in) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) * chemProfs[k].in / H;
            sumRightPrev += ental[k](f[T](i - 1, j - 1) * valProfs[T].in) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) * chemProfs[k].in / H;
            sumLeftNext  += ental[k](f[T](i, j + 1) * valProfs[T].in) * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) * chemProfs[k].in / H;
            sumRightNext += ental[k](f[T](i, j - 1) * valProfs[T].in) * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) * chemProfs[k].in / H;
        }
        ans[j]  = f[J](i - 1, j) / L
                + (1.0 - alpha) / std::pow(xi, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp2(i - 1, j + 1) * (f[J](i - 1, j + 1) - f[J](i - 1, j)) / H - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp2(i - 1, j - 1) * (f[J](i - 1, j) - f[J](i - 1, j - 1)) / H)
                + (1.0 - alpha) / std::pow(xi, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i - 1, j + 1) * yp2(i - 1, j + 1) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) / H - comp({MU, RHO, U, U}, {}, i - 1, j - 1) * yp2(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H)
                + alpha / std::pow(xi, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i, j + 1) * yp2(i, j + 1) * (f[U](i, j + 1) - f[U](i, j)) / H - comp({MU, RHO, U, U}, {}, i, j - 1) * yp2(i, j - 1) * (f[U](i, j) - f[U](i, j - 1)) / H)
                - (1.0 - alpha) * (1.0 - 1.0 / Le) / std::pow(xi, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * sumLeftPrev - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * sumRightPrev)
                - alpha * (1.0 - 1.0 / Le) / std::pow(xi, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) * sumLeftNext - comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1) * sumRightNext)
                ;
    }
    M(n - 1, n - 2) = -alpha / std::pow(Xi_0 + n * xi_h, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H;
    M(n - 1, n - 1) = 1.0 / L
                    + alpha / std::pow(Xi_0 + n * xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp2(i, n - 1) / H + comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H);
    sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    for (uint64_t k = 0; k < ental.size(); ++k) {
        // sumLeftPrev  += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumRightPrev += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumLeftNext  += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumRightNext += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) * (chemProfs[k].in - chemProfs[k].out) / H;
        // sumLeftPrev  += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) / H;
        // sumRightPrev += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) / H;
        // sumLeftNext  += ental[k](f[T](i, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 1)) / H;
        // sumRightNext += ental[k](f[T](i, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) / H;
        sumLeftPrev  += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) * chemProfs[k].in / H;
        sumRightPrev += ental[k](f[T](i - 1, n - 2) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) * chemProfs[k].in / H;
        sumLeftNext  += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 1)) * chemProfs[k].in / H;
        sumRightNext += ental[k](f[T](i, n - 2) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) * chemProfs[k].in / H;
    }
    ans[n - 1]  = f[J](i - 1, n - 1) / L
                + (1.0 - alpha) / std::pow(Xi_0 + n * xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp2(i - 1, n - 1) * (f[J](i - 1, n - 1) - f[J](i - 1, n - 1)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp2(i - 1, n - 2) * (f[J](i - 1, n - 1) - f[J](i - 1, n - 2)) / H)
                + (1.0 - alpha) / std::pow(Xi_0 + n * xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i - 1, n - 1) * yp2(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1)) / H - comp({MU, RHO, U, U}, {}, i - 1, n - 2) * yp2(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H)
                + alpha / std::pow(Xi_0 + n * xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i, n - 1) * yp2(i, n - 1) * (f[U](i, n - 1) - f[U](i, n - 1)) / H - comp({MU, RHO, U, U}, {}, i, n - 2) * yp2(i, n - 2) * (f[U](i, n - 1) - f[U](i, n - 2)) / H)
                - (1.0 - alpha) * (1.0 - 1.0 / Le) / std::pow(Xi_0 + n * xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * sumLeftPrev - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * sumRightPrev)
                - alpha * (1.0 - 1.0 / Le) / std::pow(Xi_0 + n * xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) * sumLeftNext - comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 2) * sumRightNext)
                ;
    if (alpha != 0.0) {
        ans = GaussSolveSLAE(M, ans, "J");
    }
    for (uint64_t j = 0; j < ans.size(); ++j) {
        f[J](i, j) = ans[j];
    }
    //printIterData(J);

    //P
    for (uint64_t j = 0; j < n; ++j) {
        double muSum = 0;
        f[P](i, j) = f[RHO](i, j) * f[T](i, j) * ChemicalSystem::R;
        for (uint64_t k = C1; k < f.size(); ++k) {
            muSum += chemProfs[k - C1].toPhysical(f[k](i, j));
        }
        f[P](i, j) *= valProfs[RHO].in * valProfs[T].in * muSum;
        f[P](i, j) = valProfs[P].toNormal(f[P](i, j));
    }
    //printIterData(P);

    //T
    for (uint64_t j = 0; j < f[T].size().m - 1; ++j) {
        auto func = [&] (double T) -> double {
            double ans = 0;
            for (uint64_t k = 0; k < ental.size(); ++k) {
                ans += ental[k](T) * chemProfs[k].toPhysical(f[C1 + k](i, j));
            }
            double up = f[U](i, j) * valProfs[U].in;
            double vp = f[V](i, j) * valProfs[V].in;
            double wp = f[W](i, j) * valProfs[W].in;
            ans += 1.0 / 2.0 * (up*up + vp*vp + wp*wp);
            return ans - f[J](i, j) * valProfs[J].in;
        };
        f[T](i, j) = valProfs[T].toNormal(NewtonFindT(func, f[T](i, j) * valProfs[T].in, 0.01, DiffConfig::POINTS2_ORDER1_WAY2));
        //f[T](i, j) = valProfs[T].toNormal(HalfDivFind(func, (f[T](i, j) - 0.001) * valProfs[T].in, (f[T](i, j) + 0.001) * valProfs[T].in, 0.01));
        if (std::isnan(f[T](i, j))) {
            f[T](i, j) = f[T](i - 1, j);
        }
    }
    //printIterData(T);

    //rho
    double tmpRho = 0;
    for (uint64_t j = 1; j < n - 1; ++j) {
        double xi = (j + 1) * xi_h + Xi_0;
        //3
        tmpRho  = alpha / std::pow(xi, nu) * (comp({V}, {U}, i, j + 1) - comp({V}, {U}, i, j - 1)) / (2 * H)
                + (1.0 - alpha) / std::pow(xi, nu) * (comp({V}, {U}, i - 1, j + 1) - comp({V}, {U}, i - 1, j - 1)) / (2 * H);
        tmpRho *= L;
        tmpRho += comp({}, {RHO, U}, i - 1, j) / yp(i - 1, j);
        tmpRho *= comp({U}, {}, i, j) * yp(i, j);
        f[RHO](i, j) = 1.0 / tmpRho;
    }
    f[RHO](i, 0) = (4 * f[RHO](i, 1) - f[RHO](i, 2)) / 3;
    f[RHO](i, n - 1) = f[RHO](i - 1, n - 1);
    //printIterData(RHO);

    //Ci
    for (uint64_t k = 0; k < chemCount; ++k) {
        ChemicalSystem *sys = (ChemicalSystem *)info.task;
        auto ode = sys->getODE();
        std::vector<double> oldConc(chemCount + 3, 0), newConc(chemCount + 3, 0);
        for (int j = 0; j < 5; ++j) {
            f[C1 + k](i, 0) = f[C1 + k](i - 1, 0) + 8.0 * L / (H * H * H) * ((alpha / Sc) * F2(i, 0) * f[C1 + k](i - 1, 1) + (1.0 - alpha) / Sc * F2(i - 1, 0) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)));
            f[C1 + k](i, 0) /= (1.0 + 8.0 * alpha * L / (H * H * H * Sc) * F2(i, 0));
        }
        //1
        // double coeff = L / std::pow(xi_h + Xi_0, nu) * alpha / (H*H * Sc);
        // M(0, 1) = coeff * F2(i, 0);
        // M(0, 0) = -1.0 - M(0, 1);
        // ans[0] = -f[C1 + k](i - 1, 0) - coeff / alpha * (1.0 - alpha) * (F2(i - 1, 0) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) - F2(i - 1, 0) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)));
        // for (uint64_t j = 1; j < n - 1; ++j) {
        //     double xi = (j + 1) * xi_h + Xi_0;
        //     double coeff = L / std::pow(xi, nu) * alpha / (H*H * Sc);
        //     M(j, j - 1) = coeff * F2(i, j - 1);
        //     M(j, j + 1) = coeff * F2(i, j);
        //     M(j, j) = -1.0 - M(j, j - 1) - M(j, j + 1);
        //     ans[j] = -f[C1 + k](i - 1, j) - coeff / alpha * (1.0 - alpha) * (F2(i - 1, j) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) - F2(i - 1, j - 1) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)));
        // }
        // M(n - 1, n - 2) = coeff * F2(i, n - 2);
        // M(n - 1, n - 1) = -1.0 - M(n - 1, n - 2);
        // ans[n - 1] = -f[C1 + k](i - 1, n - 1) - coeff / alpha * (1.0 - alpha) * (F2(i - 1, n - 1) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) - F2(i - 1, n - 2) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)));
        // if (alpha != 0.0) {
        //     ans = GaussSolveSLAE(M, ans, "C" + std::to_string(k + 1));
        // }
        // //записываем в следующий временной слой решение
        // for (uint64_t j = 0; j < ans.size(); ++j) {
        //     f[C1 + k](i, j) = ans[j];
        // }
        //printIterData(C1 + k);
        //////////////////////////////////////////////////////////////////
        //2
        // M(0, 0) = 1 + 8 * alpha * L / (Sc * H*H*H) * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) - comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0)) / H;
        // ans[0] = f[C1 + k](i - 1, 0) + 8 * L / (H*H*H) * (alpha / Sc * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) - comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0)) / H * f[C1 + k](i - 1, 1) * f[C1 + k](i - 1, 1)
        //          + (1.0 - alpha) / Sc * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) - comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0)) / H * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)))
        //          ;//+ chemProfs[k].toNormal(ode[k](conc)) * L;
        // for (uint64_t j = 1; j < n - 1; ++j) {
        //     double xi = (j + 1) * xi_h + Xi_0;
        //     double coeff = L / std::pow(xi, nu) * alpha / (H*H * Sc);
        //     M(j, j - 1) = (comp({MU, RHO, U}, {}, i, j) * yp2(i, j) - comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1)) / H * alpha * coeff;
        //     M(j, j + 1) = -(comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) - comp({MU, RHO, U}, {}, i, j) * yp2(i, j)) / H * alpha * coeff;
        //     M(j, j) = 1.0 + M(j, j - 1) - M(j, j + 1);
        //     ans[j] = (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp2(i - 1, j + 1) - comp({MU, RHO, U}, {}, i - 1, j) * yp2(i - 1, j)) / H * (1.0 - alpha) * coeff * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j))
        //              - (comp({MU, RHO, U}, {}, i - 1, j) * yp2(i - 1, j) - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp2(i - 1, j - 1)) / H * (1.0 - alpha) * coeff * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1))
        //              //+ chemProfs[k].toNormal(ode[k](conc)) * L
        //              + f[C1 + k](i - 1, j);
        // }
        // M(n - 1, n - 2) = M(n - 2, n - 2);
        // M(n - 1, n - 1) = M(n - 2, n - 1);
        // ans[n - 1] = ans[n - 2];
        // if (alpha != 0.0) {
        //     ans = GaussSolveSLAE(M, ans, "C" + std::to_string(k + 1) + " #2");
        // }
        // //записываем в следующий временной слой решение
        // for (uint64_t j = 0; j < ans.size(); ++j) {
        //     f[C1 + k](i, j) = ans[j];
        // }
        //3
        double coeff = 1.0 / std::pow(xi_h + Xi_0, nu) / (2*H * Sc);
        M(0, 1) = -alpha * coeff * comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H;
        M(0, 0) = 1.0 / L + alpha * coeff * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) + comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0)) / H;
        for (uint64_t l = 0; l < chemCount; ++l) {
            oldConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i - 1, 0));
            newConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i, 0));
        }
        oldConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i - 1, 0));
        oldConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i - 1, 0));
        newConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i, 0));
        newConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i, 0));
        ans[0]  = f[C1 + k](i - 1, 0) / L
                + (1.0 - alpha) * coeff * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) / H - comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) / H)
                //+ alpha * ode[k](newConc) / comp({RHO, U}, {}, i, 0) / (valProfs[RHO].in * valProfs[U].in) + (1.0 - alpha) * ode[k](oldConc) / comp({RHO, U}, {}, i - 1, 0) / (valProfs[RHO].in * valProfs[U].in)
                //+ alpha * chemProfs[k].toNormal(ode[k](newConc) / comp({RHO, U}, {}, i, 0) / (valProfs[RHO].in * valProfs[U].in)) + (1.0 - alpha) * chemProfs[k].toNormal(ode[k](oldConc) / comp({RHO, U}, {}, i - 1, 0) / (valProfs[RHO].in * valProfs[U].in))
                ;
        for (uint64_t j = 1; j < n - 1; ++j) {
            double xi = (j + 1) * xi_h + Xi_0;
            double coeff = 1.0 / std::pow(xi, nu) / (2*H * Sc);
            M(j, j - 1) = -alpha * coeff * comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H;
            M(j, j + 1) = -alpha * coeff * comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H;
            M(j, j) = 1.0 / L + alpha * coeff * (comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) + comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1)) / H;
            for (uint64_t l = 0; l < chemCount; ++l) {
                oldConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i - 1, j));
                newConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i, j));
            }
            oldConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i - 1, j));
            oldConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i - 1, j));
            newConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i, j));
            newConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i, j));
            ans[j]  = f[C1 + k](i - 1, j) / L
                    + (1.0 - alpha) * coeff * (comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) / H - comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) / H)
                    //+ alpha * ode[k](newConc) / comp({RHO, U}, {}, i, j) / (valProfs[RHO].in * valProfs[U].in) + (1.0 - alpha) * ode[k](oldConc) / comp({RHO, U}, {}, i - 1, j) / (valProfs[RHO].in * valProfs[U].in)
                    //+ alpha * chemProfs[k].toNormal(ode[k](newConc) / comp({RHO, U}, {}, i, j) / (valProfs[RHO].in * valProfs[U].in)) + (1.0 - alpha) * chemProfs[k].toNormal(ode[k](oldConc) / comp({RHO, U}, {}, i - 1, j) / (valProfs[RHO].in * valProfs[U].in))
                    ;
            
        }
        coeff = 1.0 / std::pow(xi_h * n + Xi_0, nu) / (2*H * Sc);
        M(n - 1, n - 2) = -alpha * coeff * comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H;
        M(n - 1, n - 1) = 1.0 / L + alpha * coeff * (comp({MU, RHO, U}, {}, i, n - 1) * yp2(i, n - 1) + comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2)) / H;
        for (uint64_t l = 0; l < chemCount; ++l) {
            oldConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i - 1, n - 1));
            newConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i, n - 1));
        }
        oldConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i - 1, n - 1));
        oldConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i - 1, n - 1));
        newConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i, n - 1));
        newConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i, n - 1));
        ans[n - 1]  = f[C1 + k](i - 1, n - 1) / L
                    + (1.0 - alpha) * coeff * (comp({MU, RHO, U}, {}, i, n - 1) * yp2(i, n - 1) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) / H - comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) / H)
                    //+ alpha * ode[k](newConc) / comp({RHO, U}, {}, i, n - 1) / (valProfs[RHO].in * valProfs[U].in) + (1.0 - alpha) * ode[k](oldConc) / comp({RHO, U}, {}, i - 1, n - 1) / (valProfs[RHO].in * valProfs[U].in)
                    //+ alpha * chemProfs[k].toNormal(ode[k](newConc) / comp({RHO, U}, {}, i, n - 1) / (valProfs[RHO].in * valProfs[U].in)) + (1.0 - alpha) * chemProfs[k].toNormal(ode[k](oldConc) / comp({RHO, U}, {}, i - 1, n - 1) / (valProfs[RHO].in * valProfs[U].in))
                    ;
        // std::cout << "Ci M:\n" << M << "\n";
        // for (uint64_t j = 0; j < ans.size(); ++j) {
        //     std::cout << ans[j] << " ";
        // }
        // std::cout << "\n";
        if (alpha != 0.0) {
            ans = GaussSolveSLAE(M, ans, "C" + std::to_string(k + 1));
        }
        //записываем в следующий временной слой решение
        for (uint64_t j = 0; j < ans.size(); ++j) {
            f[C1 + k](i, j) = ans[j];
        }
    }
    //chem 2
    if (fireAvail) {
        for (uint64_t j = 0; j < n; ++j) {
            ChemicalSystem *sys = (ChemicalSystem *)info.task;
            const auto &ode = sys->getODE();
            std::vector<double> conc(chemCount + 3, 0);
            for (uint64_t k = 0; k < chemCount; ++k) {
                conc[k + 1] = chemProfs[k].toPhysical(f[C1 + k](i - 1, j));
            }
            conc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i - 1, j));
            conc[chemCount + 2] = valProfs[T].toPhysical(f[T](i - 1, j));
            // std::cout << "conc\n";
            // for (auto el : conc) {
            //     std::cout << el << " ";
            // }
            // std::cout << "\n";
            // for (uint64_t l = 0; l < chemCount; ++l) {
            //     oldConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i - 1, j));
            //     newConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i, j));
            // }
            // oldConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i - 1, j));
            // oldConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i - 1, j));
            // newConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i, j));
            // newConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i, j));
            //std::cout << "ode conc\n";
            double an = 0;
            for (uint64_t k = 0; k < chemCount; ++k) {
                //f[C1 + k](i, j) += chemProfs[k].toNormal(ode[k](conc) / 1'000'000);
                //f[C1 + k](i, j) += chemProfs[k].toNormal(ode[k](conc) / comp({RHO, U}, {}, i, j) / (valProfs[RHO].in * valProfs[U].in));
                //an += ode[k](conc) / 1000000;
                //std::cout << ode[k](conc) / 1000000 << " ";
            }
            //std::cout << "\n";
            //std::cout << an << "\n";
            //f[T](i, j) = valProfs[T].toNormal(ans.back().back());
        }
    }
    for (uint64_t k = 0; k < chemCount; ++k) {
        //printIterData(C1 + k);
    }

    //J+
    // sumLeftPrev = 0.0, sumRightPrev = 0.0, sumLeftNext = 0.0, sumRightNext = 0.0;
    // M(0, 0) = 1.0 / L
    //         + alpha / std::pow(Xi_0 + xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H + comp({MU, RHO, U}, {}, i, 0) * yp2(i, 0) / H);
    // M(0, 1) = -alpha / std::pow(Xi_0 + xi_h, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, 1) * yp2(i, 1) / H;
    // sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    // for (uint64_t k = 0; k < ental.size(); ++k) {
    //     // sumLeftPrev  += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumRightPrev += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumLeftNext  += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumRightNext += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 0) - f[C1 + k](i, 0)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumLeftPrev  += ental[k](f[T](i - 1, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) / H;
    //     // sumRightPrev += ental[k](f[T](i - 1, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) / H;
    //     // sumLeftNext  += ental[k](f[T](i, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) / H;
    //     // sumRightNext += ental[k](f[T](i, 0) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, 0) - f[C1 + k](i, 0)) / H;
    //     sumLeftPrev  += ental[k](f[T](i - 1, 1) * valProfs[T].in) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) * chemProfs[k].in / H;
    //     sumRightPrev += ental[k](f[T](i - 1, 0) * valProfs[T].in) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 0)) * chemProfs[k].in / H;
    //     sumLeftNext  += ental[k](f[T](i, 1) * valProfs[T].in) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) * chemProfs[k].in / H;
    //     sumRightNext += ental[k](f[T](i, 0) * valProfs[T].in) * (f[C1 + k](i, 0) - f[C1 + k](i, 0)) * chemProfs[k].in / H;
    // }
    // ans[0]  = f[J](i - 1, 0) / L
    //         + (1.0 - alpha) / std::pow(Xi_0 + xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, 1) * yp2(i - 1, 1) * (f[J](i - 1, 1) - f[J](i - 1, 0)) / H - comp({MU, RHO, U}, {}, i - 1, 0) * yp2(i - 1, 0) * (f[J](i - 1, 0) - f[J](i - 1, 0)) / H)
    //         + (1.0 - alpha) / std::pow(Xi_0 + xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i - 1, 1) * yp2(i - 1, 1) * (f[U](i - 1, 1) - f[U](i - 1, 0)) / H - comp({MU, RHO, U, U}, {}, i - 1, 0) * yp2(i - 1, 0) * (f[U](i - 1, 0) - f[U](i - 1, 0)) / H)
    //         + alpha / std::pow(Xi_0 + xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i, 1) * yp2(i, 1) * (f[U](i, 1) - f[U](i, 0)) / H - comp({MU, RHO, U, U}, {}, i, 0) * yp2(i, 0) * (f[U](i, 0) - f[U](i, 0)) / H)
    //         - (1.0 - alpha) * (1.0 - 1.0 / Le) / std::pow(Xi_0 + xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i - 1, 1) * yp(i - 1, 1) * sumLeftPrev - comp({MU, RHO, U}, {}, i - 1, 0) * yp(i - 1, 0) * sumRightPrev)
    //         - alpha * (1.0 - 1.0 / Le) / std::pow(Xi_0 + xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i, 1) * yp(i, 1) * sumLeftNext - comp({MU, RHO, U}, {}, i, 0) * yp(i, 0) * sumRightNext)
    //         ;
    // for (uint64_t j = 1; j < n - 1; ++j) {
    //     double xi = (j + 1) * xi_h + Xi_0;
    //     M(j, j - 1) = -alpha / std::pow(xi, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H;
    //     M(j, j + 1) = -alpha / std::pow(xi, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H;
    //     M(j, j) = 1.0 / L
    //             + alpha / std::pow(xi, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i, j + 1) * yp2(i, j + 1) / H + comp({MU, RHO, U}, {}, i, j - 1) * yp2(i, j - 1) / H);
    //     sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    //     for (uint64_t k = 0; k < ental.size(); ++k) {
    //         //1
    //         //sumLeftPrev  += ental[k](f[T](i - 1, j) * valProfs[T].in) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //         //sumRightPrev += ental[k](f[T](i - 1, j) * valProfs[T].in) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //         //sumLeftNext  += ental[k](f[T](i, j) * valProfs[T].in) * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //         //sumRightNext += ental[k](f[T](i, j) * valProfs[T].in) * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //         //2 (main)
    //         // sumLeftPrev  += ental[k](f[T](i - 1, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) / H;
    //         // sumRightPrev += ental[k](f[T](i - 1, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) / H;
    //         // sumLeftNext  += ental[k](f[T](i, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) / H;
    //         // sumRightNext += ental[k](f[T](i, j) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) / H;
    //         //3
    //         sumLeftPrev  += ental[k](f[T](i - 1, j + 1) * valProfs[T].in) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) * chemProfs[k].in / H;
    //         sumRightPrev += ental[k](f[T](i - 1, j - 1) * valProfs[T].in) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) * chemProfs[k].in / H;
    //         sumLeftNext  += ental[k](f[T](i, j + 1) * valProfs[T].in) * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) * chemProfs[k].in / H;
    //         sumRightNext += ental[k](f[T](i, j - 1) * valProfs[T].in) * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) * chemProfs[k].in / H;
    //     }
    //     ans[j]  = f[J](i - 1, j) / L
    //             + (1.0 - alpha) / std::pow(xi, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp2(i - 1, j + 1) * (f[J](i - 1, j + 1) - f[J](i - 1, j)) / H - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp2(i - 1, j - 1) * (f[J](i - 1, j) - f[J](i - 1, j - 1)) / H)
    //             + (1.0 - alpha) / std::pow(xi, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i - 1, j + 1) * yp2(i - 1, j + 1) * (f[U](i - 1, j + 1) - f[U](i - 1, j)) / H - comp({MU, RHO, U, U}, {}, i - 1, j - 1) * yp2(i - 1, j - 1) * (f[U](i - 1, j) - f[U](i - 1, j - 1)) / H)
    //             + alpha / std::pow(xi, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i, j + 1) * yp2(i, j + 1) * (f[U](i, j + 1) - f[U](i, j)) / H - comp({MU, RHO, U, U}, {}, i, j - 1) * yp2(i, j - 1) * (f[U](i, j) - f[U](i, j - 1)) / H)
    //             - (1.0 - alpha) * (1.0 - 1.0 / Le) / std::pow(xi, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i - 1, j + 1) * yp(i - 1, j + 1) * sumLeftPrev - comp({MU, RHO, U}, {}, i - 1, j - 1) * yp(i - 1, j - 1) * sumRightPrev)
    //             - alpha * (1.0 - 1.0 / Le) / std::pow(xi, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i, j + 1) * yp(i, j + 1) * sumLeftNext - comp({MU, RHO, U}, {}, i, j - 1) * yp(i, j - 1) * sumRightNext)
    //             ;
    // }
    // M(n - 1, n - 2) = -alpha / std::pow(Xi_0 + n * xi_h, nu) / Pr / (2 * H) * comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H;
    // M(n - 1, n - 1) = 1.0 / L
    //                 + alpha / std::pow(Xi_0 + n * xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i, n - 1) * yp2(i, n - 1) / H + comp({MU, RHO, U}, {}, i, n - 2) * yp2(i, n - 2) / H);
    // sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    // for (uint64_t k = 0; k < ental.size(); ++k) {
    //     // sumLeftPrev  += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumRightPrev += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumLeftNext  += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 1)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumRightNext += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) * (chemProfs[k].in - chemProfs[k].out) / H;
    //     // sumLeftPrev  += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) / H;
    //     // sumRightPrev += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) / H;
    //     // sumLeftNext  += ental[k](f[T](i, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 1)) / H;
    //     // sumRightNext += ental[k](f[T](i, n - 1) * valProfs[T].in) / valProfs[T].in * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) / H;
    //     sumLeftPrev  += ental[k](f[T](i - 1, n - 1) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 1)) * chemProfs[k].in / H;
    //     sumRightPrev += ental[k](f[T](i - 1, n - 2) * valProfs[T].in) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) * chemProfs[k].in / H;
    //     sumLeftNext  += ental[k](f[T](i, n - 1) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 1)) * chemProfs[k].in / H;
    //     sumRightNext += ental[k](f[T](i, n - 2) * valProfs[T].in) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) * chemProfs[k].in / H;
    // }
    // ans[n - 1]  = f[J](i - 1, n - 1) / L
    //             + (1.0 - alpha) / std::pow(Xi_0 + n * xi_h, nu) / Pr / (2 * H) * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp2(i - 1, n - 1) * (f[J](i - 1, n - 1) - f[J](i - 1, n - 1)) / H - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp2(i - 1, n - 2) * (f[J](i - 1, n - 1) - f[J](i - 1, n - 2)) / H)
    //             + (1.0 - alpha) / std::pow(Xi_0 + n * xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i - 1, n - 1) * yp2(i - 1, n - 1) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 1)) / H - comp({MU, RHO, U, U}, {}, i - 1, n - 2) * yp2(i - 1, n - 2) * (f[U](i - 1, n - 1) - f[U](i - 1, n - 2)) / H)
    //             + alpha / std::pow(Xi_0 + n * xi_h, nu) * (1.0 - 1.0 / Pr) / (2 * H) * (comp({MU, RHO, U, U}, {}, i, n - 1) * yp2(i, n - 1) * (f[U](i, n - 1) - f[U](i, n - 1)) / H - comp({MU, RHO, U, U}, {}, i, n - 2) * yp2(i, n - 2) * (f[U](i, n - 1) - f[U](i, n - 2)) / H)
    //             - (1.0 - alpha) * (1.0 - 1.0 / Le) / std::pow(Xi_0 + n * xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i - 1, n - 1) * yp(i - 1, n - 1) * sumLeftPrev - comp({MU, RHO, U}, {}, i - 1, n - 2) * yp(i - 1, n - 2) * sumRightPrev)
    //             - alpha * (1.0 - 1.0 / Le) / std::pow(Xi_0 + n * xi_h, nu) / (2 * H) / Pr * (comp({MU, RHO, U}, {}, i, n - 1) * yp(i, n - 1) * sumLeftNext - comp({MU, RHO, U}, {}, i, n - 2) * yp(i, n - 2) * sumRightNext)
    //             ;
    // if (alpha != 0.0) {
    //     ans = GaussSolveSLAE(M, ans, "J");
    // }
    // for (uint64_t j = 0; j < ans.size(); ++j) {
    //     f[J](i, j) = ans[j];
    // }
    // //T+
    // for (uint64_t j = 0; j < f[T].size().m - 1; ++j) {
    //     auto func = [&] (double T) -> double {
    //         double ans = 0;
    //         for (uint64_t k = 0; k < ental.size(); ++k) {
    //             ans += ental[k](T) * chemProfs[k].toPhysical(f[C1 + k](i, j));
    //         }
    //         double up = f[U](i, j) * valProfs[U].in;
    //         double vp = f[V](i, j) * valProfs[V].in;
    //         double wp = f[W](i, j) * valProfs[W].in;
    //         ans += 1.0 / 2.0 * (up*up + vp*vp + wp*wp);
    //         return ans - f[J](i, j) * valProfs[J].in;
    //     };
    //     f[T](i, j) = valProfs[T].toNormal(NewtonFindT(func, f[T](i, j) * valProfs[T].in, 0.01, DiffConfig::POINTS2_ORDER1_WAY2));
    //     //f[T](i, j) = valProfs[T].toNormal(HalfDivFind(func, (f[T](i, j) - 0.001) * valProfs[T].in, (f[T](i, j) + 0.001) * valProfs[T].in, 0.01));
    //     // if (std::isnan(f[T](i, j))) {
    //     //     f[T](i, j) = f[T](i - 1, j);
    //     // }
    // }

    //Y2
    f[Y](i, 0) = 0;
    for (uint64_t j = 1; j < n; ++j) {
        //5
        auto toFind = [&] (double x) -> double {
            return std::pow(x, nu) / (f[U](i, j) * f[RHO](i, j));
        };
        double tmp = IntegralTrapeze(toFind, Xi_0, j * xi_h + Xi_0, xi_h);
        tmp *= (nu + 1);
        f[Y](i, j) = std::pow(tmp, 1.0 / (nu + 1));
    }
    //printIterData(Y);

    //MU
    for (uint64_t j = 0; j < n; ++j) {
        f[MU](i, j) = turbulence[usingTurb + 1](f, valProfs, chemProfs, R0, i, j);
    }
    //printIterData(MU);
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
void ExplNonExpl (ReportInfo &info, double theta, double R0, std::vector<Matrix<double>> &f, double Xi_0, double x_h, double xi_h, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs) {
    if (theta < 0 || theta > 1) {
        throw std::logic_error("ExplNonExpl: \"theta\" must be in range [0, 1]");
    }
    for (uint64_t i = 1; i < f[0].size().n; ++i) {
        for (int j = 0; j < f.size(); ++j) {
            for (int k = 0; k < f[0].size().m; ++k) {
                f[j](i, k) = f[j](i - 1, k);
            }
        }
        //используем итерацию для вычисления следующего временного слоя
        for (int j = 0; j < 10; ++j) {
            ExplNonExplIteration(info, theta, R0, f, Xi_0, x_h, xi_h, i, ental, valProfs, chemProfs);
        }
        std::cout << "\n////////////////////////////////////////////////////////\n";
        for (int j = 0; j < f.size(); ++j) {
            if (j > C1) {
                std::cout << "\n" << i << " before C" << j - C1 + 1 << ":\n";
            } else {
                std::cout << "\n" << i << " before " << argToStr(ARGS(j)) << ":\n";
            }
            for (uint64_t k = 0; k < f[j].size().m; ++k) {
                std::cout << std::setw(13) << f[j](i - 1, k) << " ";
            }
            if (j > C1) {
                std::cout << "\n" << i << " after C" << j - C1 + 1 << ":\n";
            } else {
                std::cout << "\n" << i << " after " << argToStr(ARGS(j)) << ":\n";
            }
            for (uint64_t k = 0; k < f[j].size().m; ++k) {
                std::cout << std::setw(13) << f[j](i, k) << " ";
            }
            std::cout << "\n";
        }
        //return;
    }
}

std::vector<Matrix<double>> SolveIBVP (ReportInfo &info, std::vector<Matrix<double>> &u, double alpha, double R0, double xh, double th, Method method, ApproxLevel approx, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs) {
    //решения задачи в зависимости от метода и метода аппроксимации
    switch (method) {
        case Method::EXPLICIT:
            ExplNonExpl(info, 0, R0, u, 0, xh, th, ental, valProfs, chemProfs); //theta = 0, получаем явную схему
            break;
        case Method::NOT_EXPLICIT:
            ExplNonExpl(info, 1.0, R0, u, 0, xh, th, ental, valProfs, chemProfs); //theta = 1, получаем неявную схему
            break;
        case Method::KRANK_NICOLAS:
            ExplNonExpl(info, alpha, R0, u, 0, xh, th, ental, valProfs, chemProfs); //theta = 0.5, получаем схему Кранка-Николоса
            break;
        default:
            break;
    }
    return u;
}