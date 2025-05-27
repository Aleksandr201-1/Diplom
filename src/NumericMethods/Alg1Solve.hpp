#ifndef ALG_1_SOLVE
#define ALG_1_SOLVE

#include <NumericMethods/Differentiation.hpp>
#include <exception>

// bool checkNewton (const std::function<double (double)> &f, double x) {
//     return f(x) * derivative(f, x, 0.01, DiffConfig::POINTS3_ORDER2_WAY1) > 0.0;
// }

double NewtonFind (const std::function<double (double)> &f, double x, double approx, DiffConfig conf) {
    // if (!checkNewton(f, x)) {
    //     throw std::runtime_error("NewtonFind: Incorrect begining point");
    // }
    uint64_t count = 1;
    double x2 = x - f(x) / derivative(f, x, approx, conf);
    while (std::abs(x2 - x) > approx) {
        x = x2;
        x2 = x - f(x) / derivative(f, x, approx, conf);
        ++count;
        if (count > 100) {
            return x2;
        }
    }
    return x2;
}

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

// bool checkSI (double a, double b, double (*function)(double), double q) {
//     double e = findEpsillon(), x = a;
//     while (x + e < b) {
//         x += e;
//         if (derivative(function, x) > q) {
//             return false;
//         }
//         if (function(x) < a || function(x) > b) {
//             return false;
//         }
//     }
//     return true;
// }

// std::pair<double, uint64_t> SIFind (double x, double a, double b, double approx, double (*fi)(double)) {
//     double q = std::max(std::abs(derivative(fi, a)), std::abs(derivative(fi, b)));
//     if (!checkSI(a, b, fi, q)) {
//         throw std::runtime_error("SI: Incorrect function fi(x). fi'(x) > q, x in (a, b) or fi(x) not in (a, b)");
//     }
//     uint64_t count = 1;
//     double x2 = fi(x);
//     while ((q / (1 - q)) * (std::abs(x2 - x)) > approx) {
//         x = x2;
//         x2 = fi(x2);
//         ++count;
//         if (count > 100) {
//             throw std::runtime_error("SI: the maximum number of iterations has been reached");
//         }
//     }
//     return {fi(x2), count + 1};
// }

#endif