#include "General.hpp"

uint64_t getOrder (const std::string &task) {
    uint64_t ans = 0;
    for (uint64_t i = 0; i < task.size(); ++i) {
        if (task[i] == 'y') {
            ++i;
            uint64_t tmp = 0;
            while (i < task.size() && task[i] == '\'') {
                ++tmp;
                ++i;
            }
            ans = std::max(ans, tmp);
        }
    }
    return ans;
}

void printVector(const std::vector<double> &vec) {
    for (double el : vec) {
        std::cout << el << " ";
    }
    std::cout << "\n";
}

double findEpsillon () {
    return std::pow(__DBL_EPSILON__, 1.0 / 3) * 10;
}

bool isEqual(double x, double y) {
    return std::fabs(x - y) < std::numeric_limits<double>::epsilon();
}

double derivative (const std::function<double(double)> &f, double x, uint64_t degree) {
    static double eps = findEpsillon();
    if (degree == 1) {
        return (f(x + eps) - f(x - eps)) / (2 * eps);
    }
    return (derivative(f, x + eps, degree - 1) - derivative(f, x - eps, degree - 1)) / (2 * eps);
}

double derivative (const std::function<double(std::vector<double> &)> &f, const std::vector<double> &X, double h, uint64_t idx) {
    std::vector<double> args(X);
    double a, b, c, d;
    args[idx] += 2 * h;
    a = f(args);
    args[idx] -= h;
    b = f(args);
    args[idx] -= 2 * h;
    c = f(args);
    args[idx] -= h;
    d = f(args);
    return (-a + 8 * b - 8 * c + d) / 12 * h;
}

std::string toString (double val, uint64_t precision) {
    return std::to_string(val).substr(0, std::to_string(val).find(".") + precision + 1);
}

std::string readLine () {
    std::string str;
    while (str.empty()) {
        std::getline(std::cin, str);
        if (!str.empty() && str[0] == '#') {
            str = "";
        }
    }
    return str;
}

double stringFix (const std::string &str) {
    double val = 0;
    std::string ans = "";
    for (uint64_t i = 0; i < str.size(); ++i) {
        if (str[i] == 'y') {
            std::string valStr;
            while (str[i] != '(') {
                ans += str[i];
                ++i;
            }
            ++i;
            while (str[i] != ')') {
                valStr += str[i];
                ++i;
            }
            val = FunctionalTree(valStr).calculate();
        } else {
            ans += str[i];
        }
    }
    return val;
}

Task getTaskInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn) {
    uint64_t idx = 0, size = 0;
    Task task;
    std::string tmp = "y";
    std::vector<std::string> args; //x y y' y'' ...
    args.push_back("x");
    for (uint64_t i = 0; i <= order; ++i) {
        args.push_back(tmp);
        tmp += "'";
    }
    for (uint64_t i = 0; i < order - 1; ++i) {
        auto func = [=] (const std::vector<double> &args) -> double {
            return args[i + 2];
        };
        task.odu_system.push_back(func);
    }

    idx = system[0].find('=');
    size = system[0].size();
    FunctionalTree y_order1(system[0].substr(0, idx), args), y_order2(system[0].substr(idx + 1, size - idx), args);
    FunctionalTree coeff = y_order1.getCoeff(order + 1);
    //std::cout << "General: " << y_order.toString(Style::LATEX) << "\n";
    auto func = [=] (const std::vector<double> &args) -> double {
        return (y_order2(args) - y_order1(args)) / coeff(args);
    };
    task.odu_system.push_back(func);

    for (uint64_t i = 1; i <= order; ++i) {
        idx = system[i].find('=');
        size = system[i].size();
        tmp = system[i].substr(0, idx);
        task.X0 = stringFix(tmp);
        if (task.X0 != X0) {
            throw std::logic_error("getSysInfo: y(X0) != X0");
        }
        //task.a = FunctionalTree(system[i].substr(idx + 1, size - idx), {}).func(0);
        task.Y.push_back(FunctionalTree(system[i].substr(idx + 1, size - idx), std::vector<std::string>()).calculate());
    }
    task.order = order;
    task.X0 = X0;
    task.Xn = Xn;

    return task;
}

Task getSysInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn) {
    uint64_t idx = 0, size = 0;
    Task task;
    std::vector<std::string> args; //x y y' y'' ...
    args.push_back("x");
    for (uint64_t i = 0; i < order; ++i) {
        uint64_t pos = system[i].find('\'');
        args.push_back(system[i].substr(0, pos));
        std::cout << "var: " << args.back() << "\n";
    }
    // for (uint64_t i = 0; i < order - 1; ++i) {
    //     auto func = [=] (const std::vector<double> &args) -> double {
    //         return args[i + 2];
    //     };
    //     task.odu_system.push_back(func);
    // }

    for (uint64_t i = 0; i < order; ++i) {
        idx = system[i].find('=');
        size = system[i].size();
        FunctionalTree y_n(system[i].substr(idx + 1, size - idx), args);
        task.odu_system.push_back(y_n);
    }

    for (uint64_t i = order; i < order * 2; ++i) {
        idx = system[i].find('=');
        size = system[i].size();
        //tmp = system[i].substr(0, idx);
        task.X0 = stringFix(system[i].substr(0, idx));
        if (task.X0 != X0) {
            throw std::logic_error("getSysInfo: y(X0) != X0");
        }
        //task.a = FunctionalTree(system[i].substr(idx + 1, size - idx), {}).func(0);
        task.Y.push_back(FunctionalTree(system[i].substr(idx + 1, size - idx), std::vector<std::string>()).calculate());
    }
    task.order = order;
    task.X0 = X0;
    task.Xn = Xn;

    return task;
}