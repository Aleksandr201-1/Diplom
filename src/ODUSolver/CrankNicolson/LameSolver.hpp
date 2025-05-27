
/*

*/
#ifndef LAME_SOLVER_HPP
#define LAME_SOLVER_HPP

#include <NumericMethods/Gauss.hpp>
#include <NumericMethods/Integral.hpp>
#include <NumericMethods/Differentiation.hpp>
#include <Math/Matrix.hpp>
#include <General/Enum.hpp>
#include <General/Underlying.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <PDFReporter/ReportGenerator.hpp>
#include <iomanip>

enum ARGS {
    U = 0,
    W,
    V,
    J,
    T,
    P,
    RHO,
    Q,
    Y,
    MU,
    C1,
    TOTAL_COUNT,
    ERROR
};

int strToArgs (const std::string &str);

std::string argToStr (ARGS arg);

//методы: явный, неявный и Кранка-Николаса
enum class Method {
    EXPLICIT,
    NOT_EXPLICIT,
    KRANK_NICOLAS
};

enum class SIDE {
    LEFT,
    RIGHT
};

//аппроксимация
enum class ApproxLevel {
    POINT2_ORDER1, //двухточечная первого порядка
    POINT2_ORDER2, //двухточечная второго порядка
    POINT3_ORDER2, //трёхточечная второго порядка
    NONE
};

struct ValueProfile {
    double in;
    std::vector<std::pair<double, double>> profile;
    double toPhysical (double val) const;
    //double toPhysical (uint64_t i) const;
    double toNormal (double val) const;
    double getPhysical (double x) const;

    ValueProfile ();
    ValueProfile (uint64_t size);
    ValueProfile (double in, uint64_t size);
};

struct ChemicalProfile {
    double in, out;
    std::vector<std::pair<double, double>> profile;
    double toPhysical (double val) const;
    //double toPhysical (uint64_t i) const;
    double toNormal (double val) const;
    double getPhysical (double x) const;
    double getNormal (double val) const;

    ChemicalProfile ();
    ChemicalProfile (uint64_t size);
    ChemicalProfile (double in, double out, uint64_t size);
};

inline int usingTurb = 3;

inline std::vector<std::function<double (const std::vector<Matrix<double>> &, const std::vector<ValueProfile> &, const std::vector<ChemicalProfile> &, double R, uint64_t, uint64_t)>> turbulence = {
    [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
        auto size = f[U].size();
        double Um = f[U](i, 0), Ub = f[U](i, size.m - 1);
        double RHOm = f[RHO](i, 0), RHOb = f[RHO](i, size.m - 1);
        return 0.0135 * f[RHO](i, j) * Um * std::sqrt(0.07 * (1.0 - Ub / Um) * (1.0 - Ub / Um)) * 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
    },
    [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
        double dU = 0;
        if (j == 0) {
            dU = (f[U](i, j) - f[U](i, j + 1)) / (f[Y](i, j) - f[Y](i, j + 1));
        } else {
            dU = (f[U](i, j) - f[U](i, j - 1)) / (f[Y](i, j) - f[Y](i, j - 1));
        }
        auto size = f[U].size();
        double l = 0.011 * f[Y](0, size.m - 1);
        return f[RHO](i, j) * l * l * std::abs(dU);
    },
    [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
        auto size = f[U].size();
        double U0 = f[U](0, 0), Un = f[U](0, size.m - 1);
        uint64_t u02 = 1, u05 = 1, u08 = 1;
        double m = Un / U0, n = 0.5 - 0.3 * std::sqrt(m / (1.0 - m));
        //double yStar = f[Y](i, u05), dyHalf = f[Y](i, u02) - f[Y](i, u08);
        double yStar = 0.2 * 300, dyHalf = (0.2 - 0.8) * 300;
        double dUm = (f[U](i, j) - Un) / (U0 - Un);
        return 0.014 * f[RHO](i, j) * f[U](i, 0) * dyHalf * std::pow(dUm, n) * (std::abs(1.0 - m) * dyHalf / yStar + 1.0 + dyHalf / yStar);
    },
    [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
        auto size = f[U].size();
        uint64_t b = size.m - 1, m = 0;
        while (f[U](0, b) == f[U](0, size.m - 1)) {
            --b;
        }
        while (f[U](0, m) == f[U](0, 0)) {
            ++m;
        }
        // while (std::abs(f[U](i, b) - f[U](0, b)) / vals[U].in > 0.01 && b < size.m - 1) {
        //     ++b;
        // }
        // while (std::abs(f[U](i, m) - f[U](0, m)) / vals[U].in > 0.01 && m > 0) {
        //     --m;
        // }
        ++b;
        --m;
        double Ub = vals[U].toNormal(f[U](i, b)), Um = vals[U].toNormal(f[U](i, m));
        double RHOb = vals[RHO].toNormal(f[RHO](i, b)), RHOm = vals[RHO].toNormal(f[RHO](i, m));
        double yHalf = 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
        if (j >= m && j <= b) {
            return vals[MU].in * (f[T](i, j) + 110.0) / (293.15 + 110.0) * std::pow(f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in + vals[MU].in / (R * f[RHO](i, j) * f[U](i, j)) + 0.135 * vals[RHO].toNormal(f[RHO](i, j)) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
        } else {
            return vals[MU].in * (f[T](i, j) + 110.0) / (293.15 + 110.0) * std::pow(f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in; 
        }
    },
    [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
        auto size = f[U].size();
        uint64_t b = size.m - 1, m = 0;
        while (f[U](0, b) == f[U](0, size.m - 1)) {
            --b;
        }
        while (f[U](0, m) == f[U](0, 0)) {
            ++m;
        }
        // while (std::abs(f[U](i, b) - f[U](0, b)) > 0.1 && b < size.m - 1) {
        //     ++b;
        // }
        // while (std::abs(f[U](i, m) - f[U](0, m)) > 0.1 && m > 0) {
        //     --m;
        // }
        ++b;
        --m;
        if (j > m && j < b) {
            double Ub = f[U](i, b), Um = f[U](i, m);
            double RHOb = f[RHO](i, b), RHOm = f[RHO](i, m);
            double yHalf = 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
            return vals[MU].in * (293.15 + 110.0) / (110.0 + vals[T].in * f[T](i, j)) * std::pow(vals[T].in * f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in + 0.0135 * f[RHO](i, j) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
        } else {
            return vals[MU].in * (293.15 + 110.0) / (110.0 + vals[T].in * f[T](i, j)) * std::pow(vals[T].in * f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in; 
        }
    }
};

/**
    * \brief Функция решения краевой задачи. Принимает следующие аргументы
    *
    *   
    * \param task Структура, содержащая задачу.
    * \param timeLimit Верхний предел времени.
    * \param xh Размер шага для X.
    * \param th Размер шага для T.
    * \param method Метод решения краевой задачи.
    * \param approx Метод аппроксимации производной.
    * 
    * \return Таблица функции U(x, t)
    *
**/
std::vector<Matrix<double>> SolveIBVP (ReportInfo &info, std::vector<Matrix<double>> &u, double alpha, double R0, double xh, double yh, double th, Method method, ApproxLevel approx, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs);


#endif