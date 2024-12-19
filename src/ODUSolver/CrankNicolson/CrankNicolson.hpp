#ifndef CRANK_NICOLSON_HPP
#define CRANK_NICOLSON_HPP

#include <NumericMethods/Gauss.hpp>
#include <NumericMethods/RUN.hpp>
#include <NumericMethods/SI.hpp>
#include <NumericMethods/LU.hpp>
#include <NumericMethods/CubeSpline.hpp>
#include <NumericMethods/Integral.hpp>
//#include <General/Interpolation.hpp>
#include <NumericMethods/Differentiation.hpp>
#include <Math/Matrix.hpp>
#include <General/Enum.hpp>
#include <General/Underlying.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <PDFReporter/ReportGenerator.hpp>

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

enum class TurbulentMode {
    MODE1 = 0,
    MODE2,
    MODE3,
    MODE4,
    MODE5,
    TOTAL_COUNT
};

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

    ChemicalProfile ();
    ChemicalProfile (uint64_t size);
    ChemicalProfile (double in, double out, uint64_t size);
};

uint64_t findDeltaU (const std::vector<Matrix<double>> &f, uint64_t i, double value);

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
            uint64_t u02 = findDeltaU(f, i, 0.2), u05 = findDeltaU(f, i, 0.5), u08 = findDeltaU(f, i, 0.8);
            double m = Un / U0, n = 0.5 - 0.3 * std::sqrt(m / (1.0 - m));
            //double yStar = f[Y](i, u05), dyHalf = f[Y](i, u02) - f[Y](i, u08);
            double yStar = 0.2 * 300, dyHalf = (0.2 - 0.8) * 300;
            double dUm = (f[U](i, j) - Un) / (U0 - Un);
            return 0.014 * f[RHO](i, j) * f[U](i, 0) * dyHalf * std::pow(dUm, n) * (std::abs(1.0 - m) * dyHalf / yStar + 1.0 + dyHalf / yStar);
        },
        [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
            auto size = f[U].size();
            uint64_t b = size.m - 1, m = 0;
            while (std::abs(f[U](i, b) - f[U](i, size.m - 1)) / vals[U].in < 0.01 && b > 1) {
                --b;
            }
            while (std::abs(f[U](i, m) - f[U](i, 0)) / vals[U].in < 0.01 && m < size.m - 1) {
                ++m;
            }
            if (b > 1) {
                ++b;
            }
            if (m < size.m - 1) {
                --m;
            }
            auto fy = [&] (uint64_t idx) -> double {
                return f[Y](i, idx);
            };
            auto fL0 = [&] (uint64_t idx) -> double {
                return vals[U].toNormal(f[U](i, idx)) * vals[W].toNormal(f[W](i, idx)) * fy(idx) * fy(idx);
            };
            auto fG0 = [&] (uint64_t idx) -> double {
                return (vals[P].toNormal(f[P](i, idx)) * vals[RHO].toNormal(f[RHO](i, idx)) * vals[U].toNormal(f[U](i, idx)) * vals[U].toNormal(f[U](i, idx))) * fy(idx);
            };
            double L0 = 2 * std::acos(-1.0) * vals[RHO].toNormal(f[RHO](i, j)) * IntegralTrapeze(fy, fL0, size.m), G0 = 2 * std::acos(-1.0) * IntegralTrapeze(fy, fG0, size.m);
            double ae = 2 * L0 / (G0 * R);
            //std::cout << "ae: " << ae << "\n";
            double k = 0.26 + 0.6 * ae, l = k * f[Y](i, j);
            double Ub = vals[U].toNormal(f[U](i, b)), Um = vals[U].toNormal(f[U](i, m));
            double RHOb = vals[RHO].toNormal(f[RHO](i, b)), RHOm = vals[RHO].toNormal(f[RHO](i, m));
            double yHalf = 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
            if (j > m && j < b) {
                return vals[MU].in * (293.15 + 110.0) / (110.0 + f[T](i, j)) * std::pow(f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in + vals[MU].in / (R * f[RHO](i, j) * f[U](i, j)) + 0.0135 * vals[RHO].toNormal(f[RHO](i, j)) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
                //return vals[MU].in / (R * f[RHO](i, j) * f[U](i, j)) + vals[MU].in / (R * f[RHO](i, j) * f[U](i, j)) + 0.0135 * vals[RHO].toNormal(f[RHO](i, j)) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
            } else {
                //return vals[MU].in / (R * f[RHO](i, j) * f[U](i, j));
                return vals[MU].in * (293.15 + 110.0) / (110.0 + f[T](i, j)) * std::pow(f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in; 
            }
            // auto size = f[U].size();
            // double Ub = vals[U].toPhysical(f[U](i, size.m - 1)), Um = vals[U].toPhysical(f[U](i, 0));
            // double RHOb = vals[RHO].toPhysical(f[RHO](i, size.m - 1)), RHOm = vals[RHO].toPhysical(f[RHO](i, 0));
            // double yHalf = 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
            // return 0.0135 * vals[RHO].toPhysical(f[RHO](i, j)) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf / 100;
            // auto size = f[U].size();
            // double Ub = f[U](i, size.m - 1), Um = f[U](i, 0);
            // double RHOb = f[RHO](i, size.m - 1), RHOm = f[RHO](i, 0);
            // double yHalf = 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
            // return 0.0135 * f[RHO](i, j) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
        },
        [] (const std::vector<Matrix<double>> &f, const std::vector<ValueProfile> &vals, const std::vector<ChemicalProfile> &chems, double R, uint64_t i, uint64_t j) -> double {
            auto size = f[U].size();
            uint64_t b = size.m - 1, m = 0;
            while (std::abs(f[U](i, b) - f[U](i, size.m - 1)) < 0.01 && b > 1) {
                --b;
                //std::cout << "B: " << b << "\n";
            }
            while (std::abs(f[U](i, m) - f[U](i, 0)) < 0.01 && m < size.m - 1) {
                ++m;
                //std::cout << "M: " << m << "\n";
            }
            if (b > 1) {
                ++b;
            }
            if (m < size.m - 1) {
                --m;
            }
            auto fy = [&] (uint64_t idx) -> double {
                return f[Y](i, idx);
            };
            auto fL0 = [&] (uint64_t idx) -> double {
                return f[U](i, idx) * f[W](i, idx) * fy(idx) * fy(idx);
            };
            auto fG0 = [&] (uint64_t idx) -> double {
                return (f[P](i, idx) * f[RHO](i, idx) * f[U](i, idx) * f[U](i, idx)) * fy(idx);
            };
            double L0 = IntegralTrapeze(fy, fL0, size.m), G0 = IntegralTrapeze(fy, fG0, size.m);
            double ae = 2 * L0 / (G0 * R);
            //std::cout << "ae: " << ae << "\n";
            double k = 0.26 + 0.6 * ae, l = k * f[Y](i, j);
            if (j > m && j < b) {
                double Ub = f[U](i, b), Um = f[U](i, m);
                double RHOb = f[RHO](i, b), RHOm = f[RHO](i, m);
                double yHalf = 1.0 / 2.0 * (RHOb * Ub + RHOm * Um);
                return vals[MU].in * (293.15 + 110.0) / (110.0 + vals[T].in * f[T](i, j)) * std::pow(vals[T].in * f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in + 0.0135 * f[RHO](i, j) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
                //return vals[MU].in / (vals[RHO].in * vals[U].in * R * f[RHO](i, j) * f[U](i, j)) + 0.0135 * f[RHO](i, j) * Um * std::sqrt(0.07 * std::pow(1.0 - Ub / Um, 2.0)) * yHalf;
            } else {
                //return vals[MU].in / (vals[RHO].in * vals[U].in * R * f[RHO](i, j) * f[U](i, j));
                return vals[MU].in * (293.15 + 110.0) / (110.0 + vals[T].in * f[T](i, j)) * std::pow(vals[T].in * f[T](i, j) / 293.15, 3.0 / 2.0) / R / vals[U].in / vals[RHO].in; 
            }
            //0.18e-4 / (0.2072*0.7363 * 0.7361*790 * 2.001)
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
std::vector<Matrix<double>> SolveIBVP (ReportInfo &info, std::vector<Matrix<double>> &u, double alpha, double R0, double xh, double th, Method method, ApproxLevel approx, const std::vector<std::function<double (double)>> &ental, const std::vector<ValueProfile> &valProfs, const std::vector<ChemicalProfile> &chemProfs);

#endif