#include <ODUSolver/CrankNicolson/LameSolver.hpp>

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
void ExplNonExplIteration (ReportInfo &info, double alpha, double R0, std::vector<Matrix<double>> &f, double x_h, double y_h, double t_h, uint64_t i, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs) {
    const double nu = 1;
    const double Sc = 0.7, Pr = 0.7, Le = Pr / Sc;
    uint64_t chemCount = f.size() - C1;
    uint64_t n = f[0].size().m;
    Matrix<double> M(n, n);
    std::vector<double> ans(n);
    Matrix<double> Mw(n - 2, n - 2);
    std::vector<double> answ(n - 2);
    const double &L = x_h, &H = y_h;
    bool fireAvail = false;

    auto ymu = [&] (uint64_t jdx) -> double {
        return std::pow(f[Y](0, jdx), nu);
    };
    auto comp = [&] (const std::vector<uint64_t> &numerator, const std::vector<uint64_t> &denominator, uint64_t n, uint64_t m) -> double {
        double ans = 1.0;
        for (uint64_t i : numerator) {
            ans *= i == Y ? std::pow(f[i](n, m), nu) : f[i](n, m);
        }
        for (uint64_t i : denominator) {
            ans /= i == Y ? std::pow(f[i](n, m), nu) : f[i](n, m);
        }
        return ans;
    };
    auto compAlpha = [&] (const std::vector<uint64_t> &numerator, const std::vector<uint64_t> &denominator, uint64_t n, uint64_t m) -> double {
        return (1.0 - alpha) * comp(numerator, denominator, n - 1, m) + alpha * comp(numerator, denominator, n, m);
    };
    auto compd = [&] (char firstDiff, const std::vector<uint64_t> &num1, const std::vector<uint64_t> &denum1, uint64_t n, uint64_t m) -> double {
        double ans = 0.0;
        if (firstDiff == 'x') {
            ans = 1.0 / (2 * L) * (comp(num1, denum1, n + 1, m) - comp(num1, denum1, n - 1, m));
        }
        if (firstDiff == 'y') {
            if (m == 0) {
                ans = 1.0 / H * (comp(num1, denum1, n, m + 1) - comp(num1, denum1, n, m));
            } else if (m == f[0].size().m - 1) {
                ans = 1.0 / H * (comp(num1, denum1, n, m) - comp(num1, denum1, n, m - 1));
            } else {
                ans = 1.0 / (2 * H) * (comp(num1, denum1, n, m + 1) - comp(num1, denum1, n, m - 1));
            }
        }
        return ans;
    };
    auto compdd = [&] (char firstDiff, const std::vector<uint64_t> &num1, const std::vector<uint64_t> &denum1, char secondDiff, const std::vector<uint64_t> &num2, const std::vector<uint64_t> &denum2, uint64_t n, uint64_t m) -> double {
        double ans = 0.0;
        if (firstDiff == 'x' && secondDiff == 'x') {
            ans = 1.0 / L * (comp(num1, denum1, n + 1, m) * (comp(num2, denum2, n + 1, m) - comp(num2, denum2, n, m)) / L - comp(num1, denum1, n - 1, m) * (comp(num2, denum2, n, m) - comp(num2, denum2, n - 1, m)) / L);
        }
        if (firstDiff == 'x' && secondDiff == 'y') {
            ans = 1.0 / (2 * L) * (comp(num1, denum1, n + 1, m) * (comp(num2, denum2, n + 1, m + 1) - comp(num2, denum2, n + 1, m - 1)) / (2 * H) - comp(num1, denum1, n - 1, m) * (comp(num2, denum2, n - 1, m + 1) - comp(num2, denum2, n - 1, m - 1)) / (2 * H));
        }
        if (firstDiff == 'y' && secondDiff == 'x') {
            ans = 1.0 / (2 * H) * (comp(num1, denum1, n, m + 1) * (comp(num2, denum2, n + 1, m + 1) - comp(num2, denum2, n - 1, m + 1)) / (2 * L) - comp(num1, denum1, n, m - 1) * (comp(num2, denum2, n + 1, m - 1) - comp(num2, denum2, n - 1, m - 1)) / (2 * L));
        }
        if (firstDiff == 'y' && secondDiff == 'y') {
            ans = 1.0 / H * (comp(num1, denum1, n, m + 1) * (comp(num2, denum2, n, m + 1) - comp(num2, denum2, n, m)) / H - comp(num1, denum1, n, m - 1) * (comp(num2, denum2, n, m) - comp(num2, denum2, n, m - 1)) / H);
        }
        return ans;
    };
    auto writeAns = [&] (uint64_t funcIdx) -> void {
        for (uint64_t j = 0; j < ans.size(); ++j) {
            f[funcIdx](i, j) = ans[j];
        }
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

    //U (V)
    for (int k = 0; k < 1; ++k) {
        // M(0, 0) = compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
        //         - 4.0 / 3.0 / (L * L) * f[MU](i, 0) //2
        //         - alpha / ymu(0) / (H * H) * (comp({MU, Y}, {}, i, 1) + comp({MU, Y}, {}, i, 1)) //3
        //         ;
        // M(0, 1) = alpha * comp({RHO, V}, {}, i, 0) / (2 * H) //0.2
        //         + alpha / ymu(0) / (H * H) * comp({MU, Y}, {}, i, 1) //3
        //         ;
        // ans[0]  = f[U](i - 2, 0) * compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
        //         - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, 0) * compd('y', {U}, {}, i - 1, 0) //0.2
        //         - compd('x', {P}, {}, i - 1, 0) //1 
        //         + 4.0 / 3.0 / (L * L) * (-f[U](i - 1, 0) * (f[MU](i, 0) + f[MU](i - 2, 0)) + f[U](i - 2, 0) * f[MU](i - 2, 0)) //2
        //         - (1.0 - alpha) / ymu(0) * compdd('y', {Y, MU}, {}, 'y', {U}, {}, i - 1, 0) //3
        //         + compdd('y', {MU}, {}, 'x', {V}, {}, i - 1, 0) //4
        //         - 2.0 / 3.0 * compdd('x', {MU}, {}, 'y', {V}, {}, i - 1, 0) //5
        //         - 2.0 / 3.0 * nu * compd('x', {MU, V}, {Y}, i - 1, 0) //6
        //         + nu * compAlpha({MU}, {Y}, i - 1, 0) * compd('x', {V}, {}, i - 1, 0) //7
        //         ;
        M(0, 0) = 1;
        M(0, 1) = -1;
        ans[0]  = 0;
        for (uint64_t j = 1; j < n - 1; ++j) {
            M(j, j - 1) = -alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                        + alpha / ymu(j) / (H * H) * comp({MU, Y}, {}, i, j - 1) //3
                        ;
            M(j, j) = compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                    - 4.0 / 3.0 / (L * L) * f[MU](i, j) //2
                    - alpha / ymu(j) / (H * H) * (comp({MU, Y}, {}, i, j + 1) + comp({MU, Y}, {}, i, j - 1)) //3
                    ;
            M(j, j + 1) = alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                        + alpha / ymu(j) / (H * H) * comp({MU, Y}, {}, i, j + 1) //3
                        ;
            ans[j] = f[U](i - 2, j) * compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                     - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, j) * compd('y', {U}, {}, i - 1, j) //0.2
                     - compd('x', {P}, {}, i - 1, j) //1 
                     + 4.0 / 3.0 / (L * L) * (-f[U](i - 1, j) * (f[MU](i, j) + f[MU](i - 2, j)) + f[U](i - 2, j) * f[MU](i - 2, j)) //2
                     - (1.0 - alpha) / ymu(j) * compdd('y', {Y, MU}, {}, 'y', {U}, {}, i - 1, j) //3
                     + compdd('y', {MU}, {}, 'x', {V}, {}, i - 1, j) //4
                     - 2.0 / 3.0 * compdd('x', {MU}, {}, 'y', {V}, {}, i - 1, j) //5
                     - 2.0 / 3.0 * nu * compd('x', {MU, V}, {Y}, i - 1, j) //6
                     + nu * compAlpha({MU}, {Y}, i - 1, j) * compd('x', {V}, {}, i - 1, j) //7
                     ;
        }
        M(n - 1, n - 2) = -alpha * comp({RHO, V}, {}, i, n - 1) / (2 * H) //0.2
                        + alpha / ymu(n - 1) / (H * H) * comp({MU, Y}, {}, i, n - 2) //3
                        ;
        M(n - 1, n - 1) = compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
                        - 4.0 / 3.0 / (L * L) * f[MU](i, n - 1) //2
                        - alpha / ymu(n - 1) / (H * H) * (comp({MU, Y}, {}, 0, n - 1) + comp({MU, Y}, {}, i, n - 2)) //3
                        ;
        ans[n - 1]  = f[U](i - 2, n - 1) * compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
                    - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, n - 1) * compd('y', {U}, {}, i - 1, n - 1) //0.2
                    - compd('x', {P}, {}, i - 1, n - 1) //1 
                    + 4.0 / 3.0 / (L * L) * (-f[U](i - 1, n - 1) * (f[MU](i, n - 1) + f[MU](i - 2, n - 1)) + f[U](i - 2, n - 1) * f[MU](i - 2, n - 1)) //2
                    - (1.0 - alpha) / ymu(n - 1) * compdd('y', {Y, MU}, {}, 'y', {U}, {}, i - 1, n - 1) //3
                    + compdd('y', {MU}, {}, 'x', {V}, {}, i - 1, n - 1) //4
                    - 2.0 / 3.0 * compdd('x', {MU}, {}, 'y', {V}, {}, i - 1, n - 1) //5
                    - 2.0 / 3.0 * nu * compd('x', {MU, V}, {Y}, i - 1, n - 1) //6
                    + nu * compAlpha({MU}, {Y}, i - 1, n - 1) * compd('x', {V}, {}, i - 1, n - 1) //7
                    ;
        if (alpha != 0.0) {
            ans = GaussSolveSLAE(M, ans, "U");
        }
        for (uint64_t j = 0; j < ans.size(); ++j) {
            f[U](i, j) = ans[j];
        }
    }

    //W (X)

    //V (V)
    M(0, 0) = 1;
    M(0, 1) = 0;
    ans[0]  = 0;
    // M(0, 0) = compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
    //         - 1.0 / (L * L) * f[MU](i, 0) //2
    //         + alpha * 4.0 / 3.0 / (H * H) * (f[MU](i, 1) + f[MU](i, 1)) //3
    //         - alpha * 2.0 / 3.0 * nu / (H * H) * (f[MU](i, 1) + f[MU](i, 1)) / ymu(0) //6
    //         ;
    // M(0, 1) = alpha * comp({RHO, V}, {}, i, 0) / (2 * H) //0.2
    //         - 4.0 / 3.0 / (H * H) * f[MU](i, 1) //3
    //         + alpha * 2.0 / 3.0 * nu / (H * H) * f[MU](i, 1) / ymu(0) //6
    //         - 2 * nu * f[MU](i, 0) / (2 * H) / ymu(1) //7
    //         ;
    // ans[0]  = f[V](i - 2, 0) * compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
    //         - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, 0) * compd('y', {V}, {}, i - 1, 0) //0.2
    //         - (1.0 - alpha) * compd('y', {P}, {}, i - 1, 0) - alpha * compd('y', {P}, {}, i, 0) //1
    //         + 1.0 / (L * L) * (-f[V](i - 1, 0) * (f[MU](i, 0) + f[MU](i - 2, 0)) + f[MU](i - 2, 0) * f[V](i - 2, 0)) //2
    //         + (1.0 - alpha) * 4.0 / 3.0 * compdd('y', {MU}, {}, 'y', {V}, {}, i - 1, 0) //3
    //         - 2.0 / 3.0 * compdd('y', {MU}, {}, 'x', {U}, {}, i - 1, 0) //4
    //         + compdd('x', {MU}, {}, 'y', {U}, {}, i - 1, 0) //5
    //         - (1.0 - alpha) * 2.0 / 3.0 * nu * compd('y', {MU, V}, {Y}, i - 1, 0) //6
    //         + (1.0 - alpha) * 2 * nu * f[MU](i - 1, 0) * compd('y', {V}, {Y}, i - 1, 0) //7
    //         + (1.0 - alpha) * nu * comp({RHO, W, W}, {Y}, i - 1, 0) + alpha * nu * comp({RHO, W, W}, {Y}, i, 0) //8
    //         ;
    for (uint64_t j = 1; j < n - 1; ++j) {
        M(j, j - 1) = -alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                    - 4.0 / 3.0 / (H * H) * f[MU](i, j - 1) //3
                    + alpha * 2.0 / 3.0 * nu / (H * H) * f[MU](i, j - 1) / ymu(j) //6
                    - 2 * nu * f[MU](i, j) / (2 * H) / ymu(j - 1) //7
                    ;
        M(j, j) = compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                - 1.0 / (L * L) * f[MU](i, j) //2
                + alpha * 4.0 / 3.0 / (H * H) * (f[MU](i, j + 1) + f[MU](i, j - 1)) //3
                - alpha * 2.0 / 3.0 * nu / (H * H) * (f[MU](i, j + 1) + f[MU](i, j - 1)) / ymu(j) //6
                ;
        M(j, j + 1) = alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                    - 4.0 / 3.0 / (H * H) * f[MU](i, j + 1) //3
                    + alpha * 2.0 / 3.0 * nu / (H * H) * f[MU](i, j + 1) / ymu(j) //6
                    - 2 * nu * f[MU](i, j) / (2 * H) / ymu(j + 1) //7
                    ;
        ans[j]  = f[V](i - 2, j) * compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, j) * compd('y', {V}, {}, i - 1, j) //0.2
                - (1.0 - alpha) * compd('y', {P}, {}, i - 1, j) - alpha * compd('y', {P}, {}, i, j) //1
                + 1.0 / (L * L) * (-f[V](i - 1, j) * (f[MU](i, j) + f[MU](i - 2, j)) + f[MU](i - 2, j) * f[V](i - 2, j)) //2
                + (1.0 - alpha) * 4.0 / 3.0 * compdd('y', {MU}, {}, 'y', {V}, {}, i - 1, j) //3
                - 2.0 / 3.0 * compdd('y', {MU}, {}, 'x', {U}, {}, i - 1, j) //4
                + compdd('x', {MU}, {}, 'y', {U}, {}, i - 1, j) //5
                - (1.0 - alpha) * 2.0 / 3.0 * nu * compd('y', {MU, V}, {Y}, i - 1, j) //6
                + (1.0 - alpha) * 2 * nu * f[MU](i - 1, j) * compd('y', {V}, {Y}, i - 1, j) //7
                + (1.0 - alpha) * nu * comp({RHO, W, W}, {Y}, i - 1, j) + alpha * nu * comp({RHO, W, W}, {Y}, i, j) //8
                ;
    }
    M(n - 1, n - 2) = 0;
    M(n - 1, n - 1) = 1;
    ans[n - 1]  = 0;
    // M(n - 1, n - 2) = -alpha * comp({RHO, V}, {}, i, n - 1) / (2 * H) //0.2
    //                 - 4.0 / 3.0 / (H * H) * f[MU](i, n - 2) //3
    //                 + alpha * 2.0 / 3.0 * nu / (H * H) * f[MU](i, n - 2) / ymu(n - 1) //6
    //                 - 2 * nu * f[MU](i, n - 1) / (2 * H) / ymu(n - 2) //7
    //                 ;
    // M(n - 1, n - 1) = compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
    //                 - 1.0 / (L * L) * f[MU](i, n - 1) //2
    //                 + alpha * 4.0 / 3.0 / (H * H) * (f[MU](i, n - 1) + f[MU](i, n - 2)) //3
    //                 - alpha * 2.0 / 3.0 * nu / (H * H) * (f[MU](i, n - 1) + f[MU](i, n - 2)) / ymu(n - 1) //6
    //                 ;
    // ans[n - 1]  = f[V](i - 2, n - 1) * compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
    //             - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, n - 1) * compd('y', {V}, {}, i - 1, n - 1) //0.2
    //             - (1.0 - alpha) * compd('y', {P}, {}, i - 1, n - 1) - alpha * compd('y', {P}, {}, i, n - 1) //1
    //             + 1.0 / (L * L) * (-f[V](i - 1, n - 1) * (f[MU](i, n - 1) + f[MU](i - 2, n - 1)) + f[MU](i - 2, n - 1) * f[V](i - 2, n - 1)) //2
    //             + (1.0 - alpha) * 4.0 / 3.0 * compdd('y', {MU}, {}, 'y', {V}, {}, i - 1, n - 1) //3
    //             - 2.0 / 3.0 * compdd('y', {MU}, {}, 'x', {U}, {}, i - 1, n - 1) //4
    //             + compdd('x', {MU}, {}, 'y', {U}, {}, i - 1, n - 1) //5
    //             - (1.0 - alpha) * 2.0 / 3.0 * nu * compd('y', {MU, V}, {Y}, i - 1, n - 1) //6
    //             + (1.0 - alpha) * 2 * nu * f[MU](i - 1, n - 1) * compd('y', {V}, {Y}, i - 1, n - 1) //7
    //             + (1.0 - alpha) * nu * comp({RHO, W, W}, {Y}, i - 1, n - 1) + alpha * nu * comp({RHO, W, W}, {Y}, i, n - 1) //8
    //             ;
    if (alpha != 0.0) {
        ans = GaussSolveSLAE(M, ans, "V");
    }
    for (uint64_t j = 0; j < ans.size(); ++j) {
        f[V](i, j) = ans[j];
    }

    //J (V)
    double sumLeftPrev = 0.0, sumRightPrev = 0.0, sumLeftNext = 0.0, sumRightNext = 0.0, sumMidNext = 0.0, sumMidPrev = 0.0;
    M(0, 0) = compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
            - 1.0 / Pr / (L * L) * f[MU](i, 0) //1
            + alpha / Pr / ymu(0) / (H * H) * (comp({Y, MU}, {}, i, 1) + comp({Y, MU}, {}, i, 1)) //2
            ;
    M(0, 1) = alpha * comp({RHO, V}, {}, i, 0) / (2 * H) //0.2
            - alpha / Pr / ymu(0) / (H * H) * comp({Y, MU}, {}, i, 1) //2
            ;
    sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    for (uint64_t k = 0; k < ental.size(); ++k) {
        sumLeftPrev  += ental[k](f[T](i - 1, 1)) * (f[C1 + k](i - 1, 1) - f[C1 + k](i - 1, 0)) / H;
        sumRightPrev += ental[k](f[T](i - 1, 1)) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 1, 1)) / H;
        sumLeftNext  += ental[k](f[T](i, 1)) * (f[C1 + k](i, 1) - f[C1 + k](i, 0)) / H;
        sumRightNext += ental[k](f[T](i, 1)) * (f[C1 + k](i, 0) - f[C1 + k](i, 1)) / H;
        sumMidNext   += ental[k](f[T](i, 0)) * (f[C1 + k](i, 0) - f[C1 + k](i - 1, 0)) / L;
        sumMidPrev   += ental[k](f[T](i - 2, 0)) * (f[C1 + k](i - 1, 0) - f[C1 + k](i - 2, 0)) / L;
    }
    ans[0]  = f[J](i - 2, 0) * compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
            - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, 0) * compd('y', {J}, {}, i - 1, 0) //0.2
            + 1.0 / Pr / (L * L) * (-f[J](i - 1, 0) * (f[MU](i, 0) + f[MU](i - 2, 0)) + f[MU](i - 2, 0) * f[J](i - 2, 0)) //1
            + (1.0 - alpha) / ymu(0) / Pr * compdd('y', {Y, MU}, {}, 'y', {J}, {}, i - 1, 0) //2
            + (4.0 / 3.0 - 1.0 / Pr) * compdd('x', {MU, U}, {}, 'x', {U}, {}, i - 1, 0) //3
            + (1.0 - 1.0 / Pr) * compdd('x', {MU, W}, {}, 'x', {W}, {}, i - 1, 0) //4
            + (1.0 - 1.0 / Pr) * compdd('x', {MU, V}, {}, 'x', {V}, {}, i - 1, 0) //5
            + (1.0 - alpha) * (1.0 - 1.0 / Pr) / ymu(0) * compdd('y', {Y, MU, U}, {}, 'y', {U}, {}, i - 1, 0) + alpha * (1.0 - 1.0 / Pr) / ymu(0) * compdd('y', {Y, MU, U}, {}, 'y', {U}, {}, i, 0) //6
            + (1.0 - alpha) * (4.0 / 3.0 - 1.0 / Pr) / ymu(0) * compdd('y', {Y, MU, V}, {}, 'y', {V}, {}, i - 1, 0) + alpha * (4.0 / 3.0 - 1.0 / Pr) / ymu(0) * compdd('y', {Y, MU, V}, {}, 'y', {V}, {}, i, 0) //7
            + (1.0 - alpha) * (1.0 - 1.0 / Pr) / ymu(0) * compdd('y', {Y, MU, W}, {}, 'y', {W}, {}, i - 1, 0) + alpha * (1.0 - 1.0 / Pr) / ymu(0) * compdd('y', {Y, MU, W}, {}, 'y', {W}, {}, i, 0) //8
            + compdd('x', {MU, V}, {}, 'y', {U}, {}, i - 1, 0) //9
            + 1.0 / ymu(0) * compdd('y', {MU, Y, U}, {}, 'x', {V}, {}, i - 1, 0) //10
            - 2.0 / 3.0 * nu * compd('x', {MU, U, V}, {Y}, i - 1, 0) //11
            - (1.0 - alpha) * 2.0 / 3.0 * nu / ymu(0) * compd('y', {MU, V, V}, {}, i - 1, 0) - alpha * 2.0 / 3.0 * nu / ymu(0) * compd('y', {MU, V, V}, {}, i, 0) //12
            - (1.0 - alpha) * nu / ymu(0) * compd('y', {MU, W, W}, {}, i - 1, 0) - alpha * nu / ymu(0) * compd('y', {MU, W, W}, {}, i, 0) //13
            - 2.0 / 3.0 * compdd('x', {MU, U}, {}, 'y', {V}, {}, i - 1, 0) //14
            - 2.0 / 3.0 / ymu(0) * compdd('y', {MU, V, Y}, {}, 'x', {U}, {}, i - 1, 0) //15
            - (1.0 - 1.0 / Le) / (L * Pr) * (f[MU](i, 0) * sumMidNext - f[MU](i, 0) * sumMidPrev)//16
            - (1.0 - alpha) * (1.0 - 1.0 / Le) / (H * Pr * ymu(0)) * (comp({Y, MU}, {}, i - 1, 1) * sumLeftPrev - comp({Y, MU}, {}, i - 1, 1) * sumRightPrev)//17
            - alpha * (1.0 - 1.0 / Le) / (H * Pr * ymu(0)) * (comp({Y, MU}, {}, i, 1) * sumLeftNext - comp({Y, MU}, {}, i, 1) * sumRightNext)//17
            - (1.0 - alpha) * comp({RHO, Q}, {}, i - 1, 0) - alpha * comp({RHO, Q}, {}, i, 0) //18
            ;
    for (uint64_t j = 1; j < n - 1; ++j) {
        M(j, j - 1) = -alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                    - alpha / Pr / ymu(j) / (H * H) * comp({Y, MU}, {}, i, j - 1) //2
                    ;
        M(j, j) = compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                - 1.0 / Pr / (L * L) * f[MU](i, j) //1
                + alpha / Pr / ymu(j) / (H * H) * (comp({Y, MU}, {}, i, j + 1) + comp({Y, MU}, {}, i, j - 1)) //2
                ;
        M(j, j + 1) = alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                    - alpha / Pr / ymu(j) / (H * H) * comp({Y, MU}, {}, i, j + 1) //2
                    ;
        sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
        for (uint64_t k = 0; k < ental.size(); ++k) {
            sumLeftPrev  += ental[k](f[T](i - 1, j + 1)) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) / H;
            sumRightPrev += ental[k](f[T](i - 1, j - 1)) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) / H;
            sumLeftNext  += ental[k](f[T](i, j + 1)) * (f[C1 + k](i, j + 1) - f[C1 + k](i, j)) / H;
            sumRightNext += ental[k](f[T](i, j - 1)) * (f[C1 + k](i, j) - f[C1 + k](i, j - 1)) / H;
            sumMidNext   += ental[k](f[T](i, j)) * (f[C1 + k](i, j) - f[C1 + k](i - 1, j)) / L;
            sumMidPrev   += ental[k](f[T](i - 2, j)) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 2, j)) / L;
        }
        ans[j]  = f[J](i - 2, j) * compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, j) * compd('y', {J}, {}, i - 1, j) //0.2
                + 1.0 / Pr / (L * L) * (-f[J](i - 1, j) * (f[MU](i, j) + f[MU](i - 2, j)) + f[MU](i - 2, j) * f[J](i - 2, j)) //1
                + (1.0 - alpha) / ymu(j) / Pr * compdd('y', {Y, MU}, {}, 'y', {J}, {}, i - 1, j) //2
                + (4.0 / 3.0 - 1.0 / Pr) * compdd('x', {MU, U}, {}, 'x', {U}, {}, i - 1, j) //3
                + (1.0 - 1.0 / Pr) * compdd('x', {MU, W}, {}, 'x', {W}, {}, i - 1, j) //4
                + (1.0 - 1.0 / Pr) * compdd('x', {MU, V}, {}, 'x', {V}, {}, i - 1, j) //5
                + (1.0 - alpha) * (1.0 - 1.0 / Pr) / ymu(j) * compdd('y', {Y, MU, U}, {}, 'y', {U}, {}, i - 1, j) + alpha * (1.0 - 1.0 / Pr) / ymu(j) * compdd('y', {Y, MU, U}, {}, 'y', {U}, {}, i, j) //6
                + (1.0 - alpha) * (4.0 / 3.0 - 1.0 / Pr) / ymu(j) * compdd('y', {Y, MU, V}, {}, 'y', {V}, {}, i - 1, j) + alpha * (4.0 / 3.0 - 1.0 / Pr) / ymu(j) * compdd('y', {Y, MU, V}, {}, 'y', {V}, {}, i, j) //7
                + (1.0 - alpha) * (1.0 - 1.0 / Pr) / ymu(j) * compdd('y', {Y, MU, W}, {}, 'y', {W}, {}, i - 1, j) + alpha * (1.0 - 1.0 / Pr) / ymu(j) * compdd('y', {Y, MU, W}, {}, 'y', {W}, {}, i, j) //8
                + compdd('x', {MU, V}, {}, 'y', {U}, {}, i - 1, j) //9
                + 1.0 / ymu(j) * compdd('y', {MU, Y, U}, {}, 'x', {V}, {}, i - 1, j) //10
                - 2.0 / 3.0 * nu * compd('x', {MU, U, V}, {Y}, i - 1, j) //11
                - (1.0 - alpha) * 2.0 / 3.0 * nu / ymu(j) * compd('y', {MU, V, V}, {}, i - 1, j) - alpha * 2.0 / 3.0 * nu / ymu(j) * compd('y', {MU, V, V}, {}, i, j) //12
                - (1.0 - alpha) * nu / ymu(j) * compd('y', {MU, W, W}, {}, i - 1, j) - alpha * nu / ymu(j) * compd('y', {MU, W, W}, {}, i, j) //13
                - 2.0 / 3.0 * compdd('x', {MU, U}, {}, 'y', {V}, {}, i - 1, j) //14
                - 2.0 / 3.0 / ymu(j) * compdd('y', {MU, V, Y}, {}, 'x', {U}, {}, i - 1, j) //15
                - (1.0 - 1.0 / Le) / (L * Pr) * (f[MU](i, j) * sumMidNext - f[MU](i, j) * sumMidPrev)//16
                - (1.0 - alpha) * (1.0 - 1.0 / Le) / (H * Pr * ymu(j)) * (comp({Y, MU}, {}, i - 1, j + 1) * sumLeftPrev - comp({Y, MU}, {}, i - 1, j - 1) * sumRightPrev)//17
                - alpha * (1.0 - 1.0 / Le) / (H * Pr * ymu(j)) * (comp({Y, MU}, {}, i, j + 1) * sumLeftNext - comp({Y, MU}, {}, i, j - 1) * sumRightNext)//17
                - (1.0 - alpha) * comp({RHO, Q}, {}, i - 1, j) - alpha * comp({RHO, Q}, {}, i, j) //18
                ;
    }
    M(n - 1, n - 2) = -alpha * comp({RHO, V}, {}, i, n - 1) / (2 * H) //0.2
                    - alpha / Pr / ymu(n - 1) / (H * H) * comp({Y, MU}, {}, i, n - 2) //2
                    ;
    M(n - 1, n - 1) = compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
                    - 1.0 / Pr / (L * L) * f[MU](i, n - 1) //1
                    + alpha / Pr / ymu(n - 1) / (H * H) * (comp({Y, MU}, {}, i, n - 1) + comp({Y, MU}, {}, i, n - 2)) //2
                    ;
    sumLeftPrev = sumRightPrev = sumLeftNext = sumRightNext = 0;
    for (uint64_t k = 0; k < ental.size(); ++k) {
        sumLeftPrev  += ental[k](f[T](0, n - 1)) * (f[C1 + k](0, n - 1) - f[C1 + k](i - 1, n - 1)) / H;
        sumRightPrev += ental[k](f[T](i - 1, n - 2)) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 1, n - 2)) / H;
        sumLeftNext  += ental[k](f[T](0, n - 1)) * (f[C1 + k](0, n - 1) - f[C1 + k](i, n - 1)) / H;
        sumRightNext += ental[k](f[T](i, n - 2)) * (f[C1 + k](i, n - 1) - f[C1 + k](i, n - 2)) / H;
        sumMidNext   += ental[k](f[T](i, n - 1)) * (f[C1 + k](i, n - 1) - f[C1 + k](i - 1, n - 1)) / L;
        sumMidPrev   += ental[k](f[T](i - 2, n - 1)) * (f[C1 + k](i - 1, n - 1) - f[C1 + k](i - 2, n - 1)) / L;
    }
    ans[n - 1]  = f[J](i - 2, n - 1) * compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
                - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, n - 1) * compd('y', {J}, {}, i - 1, n - 1) //0.2
                + 1.0 / Pr / (L * L) * (-f[J](i - 1, n - 1) * (f[MU](i, n - 1) + f[MU](i - 2, n - 1)) + f[MU](i - 2, n - 1) * f[J](i - 2, n - 1)) //1
                + (1.0 - alpha) / ymu(n - 1) / Pr * compdd('y', {Y, MU}, {}, 'y', {J}, {}, i - 1, n - 1) //2
                + (4.0 / 3.0 - 1.0 / Pr) * compdd('x', {MU, U}, {}, 'x', {U}, {}, i - 1, n - 1) //3
                + (1.0 - 1.0 / Pr) * compdd('x', {MU, W}, {}, 'x', {W}, {}, i - 1, n - 1) //4
                + (1.0 - 1.0 / Pr) * compdd('x', {MU, V}, {}, 'x', {V}, {}, i - 1, n - 1) //5
                + (1.0 - alpha) * (1.0 - 1.0 / Pr) / ymu(n - 1) * compdd('y', {Y, MU, U}, {}, 'y', {U}, {}, i - 1, n - 1) + alpha * (1.0 - 1.0 / Pr) / ymu(n - 1) * compdd('y', {Y, MU, U}, {}, 'y', {U}, {}, i, n - 1) //6
                + (1.0 - alpha) * (4.0 / 3.0 - 1.0 / Pr) / ymu(n - 1) * compdd('y', {Y, MU, V}, {}, 'y', {V}, {}, i - 1, n - 1) + alpha * (4.0 / 3.0 - 1.0 / Pr) / ymu(n - 1) * compdd('y', {Y, MU, V}, {}, 'y', {V}, {}, i, n - 1) //7
                + (1.0 - alpha) * (1.0 - 1.0 / Pr) / ymu(n - 1) * compdd('y', {Y, MU, W}, {}, 'y', {W}, {}, i - 1, n - 1) + alpha * (1.0 - 1.0 / Pr) / ymu(n - 1) * compdd('y', {Y, MU, W}, {}, 'y', {W}, {}, i, n - 1) //8
                + compdd('x', {MU, V}, {}, 'y', {U}, {}, i - 1, n - 1) //9
                + 1.0 / ymu(n - 1) * compdd('y', {MU, Y, U}, {}, 'x', {V}, {}, i - 1, n - 1) //10
                - 2.0 / 3.0 * nu * compd('x', {MU, U, V}, {Y}, i - 1, n - 1) //11
                - (1.0 - alpha) * 2.0 / 3.0 * nu / ymu(n - 1) * compd('y', {MU, V, V}, {}, i - 1, n - 1) - alpha * 2.0 / 3.0 * nu / ymu(n - 1) * compd('y', {MU, V, V}, {}, i, n - 1) //12
                - (1.0 - alpha) * nu / ymu(n - 1) * compd('y', {MU, W, W}, {}, i - 1, n - 1) - alpha * nu / ymu(n - 1) * compd('y', {MU, W, W}, {}, i, n - 1) //13
                - 2.0 / 3.0 * compdd('x', {MU, U}, {}, 'y', {V}, {}, i - 1, n - 1) //14
                - 2.0 / 3.0 / ymu(n - 1) * compdd('y', {MU, V, Y}, {}, 'x', {U}, {}, i - 1, n - 1) //15
                - (1.0 - 1.0 / Le) / (L * Pr) * (f[MU](i, n - 1) * sumMidNext - f[MU](i, n - 1) * sumMidPrev)//16
                - (1.0 - alpha) * (1.0 - 1.0 / Le) / (H * Pr * ymu(n - 1)) * (comp({Y, MU}, {}, i - 1, n - 1) * sumLeftPrev - comp({Y, MU}, {}, i - 1, n - 2) * sumRightPrev)//17
                - alpha * (1.0 - 1.0 / Le) / (H * Pr * ymu(n - 1)) * (comp({Y, MU}, {}, i, n - 1) * sumLeftNext - comp({Y, MU}, {}, i, n - 2) * sumRightNext)//17
                - (1.0 - alpha) * comp({RHO, Q}, {}, i - 1, n - 1) - alpha * comp({RHO, Q}, {}, i, n - 1) //18
                ;
    if (alpha != 0.0) {
        ans = GaussSolveSLAE(M, ans, "J");
    }
    for (uint64_t j = 0; j < ans.size(); ++j) {
        f[J](i, j) = ans[j];
    }

    //T (V)
    for (uint64_t j = 0; j < f[T].size().m - 1; ++j) {
        auto func = [&] (double T) -> double {
            double ans = 0;
            for (uint64_t k = 0; k < ental.size(); ++k) {
                //std::cout << "ental " << k << ": " << ental[k](T) << "\n";
                ans += ental[k](T) * f[C1 + k](i, j);
            }
            double up = f[U](i, j);
            double vp = f[V](i, j);
            double wp = f[W](i, j);
            ans += 1.0 / 2.0 * (up*up + vp*vp + wp*wp);
            return -ans + f[J](i, j);
        };
        double tmpT = f[T](i, j);
        f[T](i, j) = NewtonFindT(func, f[T](0, 0), 0.01, DiffConfig::POINTS2_ORDER1_WAY2);
        if (std::isnan(f[T](i, j))) {
            std::cout << "T: NaN encountered at point (" << i << ", " << j << ")! Exit\n"
                      << "\tU: " << f[U](i, j) / valProfs[U].in << "\n"
                      << "\tV: " << f[V](i, j) / valProfs[V].in << "\n"
                      << "\tW: " << f[W](i, j) / valProfs[W].in << "\n"
                      << "\tJ: " << f[J](i, j) / valProfs[J].in << "\n"
                      << "\tT: " << tmpT << "\n"
                      << "\tf(T): " << func(f[T](i, j)) << "\n";
            std::cout.flush();
            exit(1);
            f[T](i, j) = f[T](i - 1, j);
        }
    }

    //rho (V)
    double tmpRho = 0;
    // M(0, 0) = comp({U, Y}, {}, i, 0) / (2 * L) - alpha * comp({V, Y}, {}, i, 0) / H;
    // M(0, 1) = alpha * comp({V, Y}, {}, i, 1) / H;
    // ans[0]  = comp({RHO, U, Y}, {}, i - 2, 0) / (2 * L)
    //         - (1.0 - alpha) * compd('y', {RHO, V, Y}, {}, i - 1, 0);
    M(0, 0) = 1;
    M(0, 1) = -1;
    ans[0]  = 0;
    for (uint64_t j = 1; j < n - 1; ++j) {
        M(j, j - 1) = -alpha * comp({V, Y}, {}, i, j - 1) / (2 * H);
        M(j, j) = comp({U, Y}, {}, i, j) / (2 * L);
        M(j, j + 1) = alpha * comp({V, Y}, {}, i, j + 1) / (2 * H);
        ans[j]  = comp({RHO, U, Y}, {}, i - 2, j) / (2 * L)
                - (1.0 - alpha) * compd('y', {RHO, V, Y}, {}, i - 1, j);
    }
    // M(n - 1, n - 2) = -alpha * comp({V, Y}, {}, i, n - 2) / H;
    // M(n - 1, n - 1) = comp({U, Y}, {}, i, n - 1) / (2 * L) + alpha * comp({V, Y}, {}, i, n - 1) / H;
    // ans[n - 1]  = comp({RHO, U, Y}, {}, i - 2, n - 1) / (2 * L)
    //             - (1.0 - alpha) * compd('y', {RHO, V, Y}, {}, i - 1, n - 1);
    M(n - 1, n - 2) = 0;
    M(n - 1, n - 1) = 1;
    ans[n - 1]  = f[RHO](0, n - 1);
    if (alpha != 0.0) {
        ans = GaussSolveSLAE(M, ans, "RHO");
    }
    for (uint64_t j = 0; j < ans.size(); ++j) {
        f[RHO](i, j) = ans[j];
    }

    //P (V)
    for (uint64_t j = 0; j < n; ++j) {
        double muSum = 0;
        f[P](i, j) = f[RHO](i, j) * f[T](i, j) * ChemicalSystem::R;
        for (uint64_t k = C1; k < f.size(); ++k) {
            muSum += f[k](i, j);
        }
        f[P](i, j) *= muSum;
    }

    //Ci (V)
    for (uint64_t k = 0; k < chemCount; ++k) {
        //ChemicalSystem *sys = (ChemicalSystem *)info.task;
        //auto ode = sys->getODE();
        //1
        // M(0, 0) = compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
        //         - 1.0 / (L * L * Sc) * f[MU](i, 0) //1
        //         + alpha / ymu(0) / (H * H * Sc) * (comp({Y, MU}, {}, i, 1) + comp({Y, MU}, {}, 0, 1)) //2
        //         ;
        // M(0, 1) = alpha * comp({RHO, V}, {}, i, 0) / (2 * H) //0.2
        //         - alpha / ymu(0) / (H * H * Sc) * ymu(1) * f[MU](i, 1) //2
        //         ;
        // ans[0]  = f[C1 + k](i - 2, 0) * compAlpha({RHO, U}, {}, i, 0) / (2 * L) //0.1
        //         - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, 0) * compd('y', {C1 + k}, {}, i - 1, 0) //0.2
        //         + 1.0 / (L * L * Sc) * (-(f[MU](i, 0) + f[MU](i - 2, 0)) * f[C1 + k](i - 1, 0) + f[MU](i - 2, 0) * f[C1 + k](i - 2, 0)) //1
        //         + (1.0 - alpha) / ymu(0) / Sc * compdd('y', {Y, MU}, {}, 'y', {C1 + k}, {}, i - 1, 0) //2
        //         ;
        M(0, 0) = 1;
        M(0, 1) = -1;
        ans[0]  = 0;
        for (uint64_t j = 1; j < n - 1; ++j) {
            M(j, j - 1) = -alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                        - alpha / ymu(j) / (H * H * Sc) * ymu(j - 1) * f[MU](i, j - 1) //2
                        ;
            M(j, j) = compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                    - 1.0 / (L * L * Sc) * f[MU](i, j) //1
                    + alpha / ymu(j) / (H * H * Sc) * (comp({Y, MU}, {}, i, j + 1) + comp({Y, MU}, {}, i, j - 1)) //2
                    ;
            M(j, j + 1) = alpha * comp({RHO, V}, {}, i, j) / (2 * H) //0.2
                        - alpha / ymu(j) / (H * H * Sc) * ymu(j + 1) * f[MU](i, j + 1) //2
                        ;
            ans[j]  = f[C1 + k](i - 2, j) * compAlpha({RHO, U}, {}, i, j) / (2 * L) //0.1
                    - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, j) * compd('y', {C1 + k}, {}, i - 1, j) //0.2
                    + 1.0 / (L * L * Sc) * (-(f[MU](i, j) + f[MU](i - 2, j)) * f[C1 + k](i - 1, j) + f[MU](i - 2, j) * f[C1 + k](i - 2, j)) //1
                    + (1.0 - alpha) / ymu(j) / Sc * compdd('y', {Y, MU}, {}, 'y', {C1 + k}, {}, i - 1, j) //2
                    ;
            // for (uint64_t l = 0; l < chemCount; ++l) {
            //     oldConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i - 1, j));
            //     newConc[l + 1] = chemProfs[l].toPhysical(f[C1 + l](i, j));
            // }
            // oldConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i - 1, j));
            // oldConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i - 1, j));
            // newConc[chemCount + 1] = valProfs[RHO].toPhysical(f[RHO](i, j));
            // newConc[chemCount + 2] = valProfs[T].toPhysical(f[T](i, j));
            // ans[j]  = f[C1 + k](i - 1, j) / L
            //         + (1.0 - alpha) * coeff * (comp({MU, RHO, U}, {}, i, j + 1) * ymu2(i, j + 1) * (f[C1 + k](i - 1, j + 1) - f[C1 + k](i - 1, j)) / H - comp({MU, RHO, U}, {}, i, j - 1) * ymu2(i, j - 1) * (f[C1 + k](i - 1, j) - f[C1 + k](i - 1, j - 1)) / H)
            //         //+ alpha * ode[k](newConc) / comp({RHO, U}, {}, i, j) / (valProfs[RHO].in * valProfs[U].in) + (1.0 - alpha) * ode[k](oldConc) / comp({RHO, U}, {}, i - 1, j) / (valProfs[RHO].in * valProfs[U].in)
            //         //+ alpha * chemProfs[k].toNormal(ode[k](newConc) / comp({RHO, U}, {}, i, j) / (valProfs[RHO].in * valProfs[U].in)) + (1.0 - alpha) * chemProfs[k].toNormal(ode[k](oldConc) / comp({RHO, U}, {}, i - 1, j) / (valProfs[RHO].in * valProfs[U].in))
            //         ;
            
        }
        // M(n - 1, n - 2) = - alpha * comp({RHO, V}, {}, i, n - 1) / (2 * H) //0.2
        //                 - alpha / ymu(n - 1) / (H * H * Sc) * ymu(n - 2) * f[MU](i, n - 2) //2
        //                 ;
        // M(n - 1, n - 1) = compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
        //                 - 1.0 / (L * L * Sc) * f[MU](i, n - 1) //1
        //                 + alpha / ymu(n - 1) / (H * H * Sc) * (comp({Y, MU}, {}, i, n - 1) + comp({Y, MU}, {}, i, n - 2)) //2
        //                 ;
        // ans[n - 1]  = f[C1 + k](i - 2, n - 1) * compAlpha({RHO, U}, {}, i, n - 1) / (2 * L) //0.1
        //             - (1.0 - alpha) * comp({RHO, V}, {}, i - 1, n - 1) * compd('y', {C1 + k}, {}, i - 1, n - 1) //0.2
        //             + 1.0 / (L * L * Sc) * (-(f[MU](i, n - 1) + f[MU](i - 2, n - 1)) * f[C1 + k](i - 1, n - 1) + f[MU](i - 2, n - 1) * f[C1 + k](i - 2, n - 1)) //1
        //             + (1.0 - alpha) / ymu(n - 1) / Sc * compdd('y', {Y, MU}, {}, 'y', {C1 + k}, {}, i - 1, n - 1) //2
        //             ;
        M(n - 1, n - 2) = 0;
        M(n - 1, n - 1) = 1;
        ans[n - 1]  = f[C1 + k](0, n - 1);
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
            //std::cout << "ode conc\n";
            double an = 0;
            for (uint64_t k = 0; k < chemCount; ++k) {
                //f[C1 + k](i, j) += chemProfs[k].toNormal(ode[k](conc) / 1'000'000);
                //f[C1 + k](i, j) += chemProfs[k].toNormal(ode[k](conc) / comp({RHO, U}, {}, i, j) / (valProfs[RHO].in * valProfs[U].in));
                //an += ode[k](conc) / 1000000;
                //std::cout << ode[k](conc) / 1000000 << " ";
            }
        }
    }

    //MU (V)
    for (uint64_t j = 0; j < n; ++j) {
        f[MU](i, j) = turbulence[usingTurb](f, valProfs, chemProfs, R0, i, j);
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
void ExplNonExpl (ReportInfo &info, double theta, double R0, std::vector<Matrix<double>> &f, double x_h, double y_h, double t_h, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs) {
    if (theta < 0 || theta > 1) {
        throw std::logic_error("ExplNonExpl: \"theta\" must be in range [0, 1]");
    }
    for (uint64_t i = 2; i < f[0].size().n; ++i) {
        for (int j = 0; j < f.size(); ++j) {
            for (int k = 0; k < f[0].size().m; ++k) {
                f[j](i, k) = f[j](i - 1, k);
            }
        }
        //используем итерацию для вычисления следующего временного слоя
        for (int j = 0; j < 10; ++j) {
            ExplNonExplIteration(info, theta, R0, f, x_h, y_h, t_h, i, ental, valProfs, chemProfs);
        }
        std::cout << "\n////////////////////////////////////////////////////////\n";
        for (int j = 0; j < f.size(); ++j) {
            if (j >= C1) {
                std::cout << "\n" << i << " before C" << j - C1 + 1 << ":\n";
            } else {
                std::cout << "\n" << i << " before " << argToStr(ARGS(j)) << ":\n";
            }
            for (uint64_t k = 0; k < f[j].size().m; ++k) {
                if (j < C1) {
                    if (j != Y && j != MU) {
                        std::cout << std::setw(13) << f[j](i - 1, k) / valProfs[j].in << " ";
                    } else {
                        std::cout << std::setw(13) << f[j](i - 1, k) << " ";
                    }
                } else {
                    std::cout << std::setw(13) << f[j](i - 1, k) / chemProfs[j - C1].in << " ";
                }
            }
            if (j >= C1) {
                std::cout << "\n" << i << " after C" << j - C1 + 1 << ":\n";
            } else {
                std::cout << "\n" << i << " after " << argToStr(ARGS(j)) << ":\n";
            }
            for (uint64_t k = 0; k < f[j].size().m; ++k) {
                if (j < C1) {
                    if (j != Y && j != MU) {
                        std::cout << std::setw(13) << f[j](i, k) / valProfs[j].in << " ";
                    } else {
                        std::cout << std::setw(13) << f[j](i, k) << " ";
                    }
                } else {
                    std::cout << std::setw(13) << f[j](i, k) / chemProfs[j - C1].in << " ";
                }
            }
            std::cout << "\n";
        }
    }
}

std::vector<Matrix<double>> SolveIBVP (ReportInfo &info, std::vector<Matrix<double>> &u, double alpha, double R0, double xh, double yh, double th, Method method, ApproxLevel approx, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs) {
    //решения задачи в зависимости от метода и метода аппроксимации
    switch (method) {
        case Method::EXPLICIT:
            ExplNonExpl(info, 0, R0, u, xh, yh, th, ental, valProfs, chemProfs); //theta = 0, получаем явную схему
            break;
        case Method::NOT_EXPLICIT:
            ExplNonExpl(info, 1.0, R0, u, xh, yh, th, ental, valProfs, chemProfs); //theta = 1, получаем неявную схему
            break;
        case Method::KRANK_NICOLAS:
            ExplNonExpl(info, alpha, R0, u, xh, yh, th, ental, valProfs, chemProfs); //theta = 0.5, получаем схему Кранка-Николоса
            break;
        default:
            break;
    }
    return u;
}