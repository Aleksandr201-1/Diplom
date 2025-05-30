#include "KoshiTask.hpp"

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
            val = FuncMaker(valStr).calculate();
        } else {
            ans += str[i];
        }
    }
    return val;
}

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

KoshiTask::KoshiTask () : X0(0), Xn(0) {}

KoshiTask::~KoshiTask () {}

void KoshiTask::setTaskInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn) {
    ode_system.clear();
    Y.clear();
    uint64_t idx = 0, size = 0;
    std::string tmp = "y";
    std::vector<std::string> args; //x y y' y'' ...
    args.push_back("x");
    for (uint64_t i = 0; i < order; ++i) {
        args.push_back(tmp);
        tmp += "'";
    }
    for (uint64_t i = 0; i < order - 1; ++i) {
        auto func = [=] (const std::vector<double> &args) -> double {
            return args[i + 2];
        };
        ode_system.push_back(func);
    }

    idx = system[0].find('=');
    size = system[0].size();
    // FuncMaker y_order1(system[0].substr(0, idx), args), y_order2(system[0].substr(idx + 1, size - idx), args);
    // FuncMaker coeff = y_order1.getCoeff(order + 1);
    // auto func = [=] (const std::vector<double> &args) -> double {
    //     return (y_order2(args) - y_order1(args)) / coeff(args);
    // };
    FuncMaker y_order2(system[0].substr(idx + 1, size - idx), args);
    //FuncMaker coeff = y_order1.getCoeff(order + 1);
    auto func = [=] (const std::vector<double> &args) -> double {
        return y_order2(args);
    };
    ode_system.push_back(func);

    for (uint64_t i = 1; i <= order; ++i) {
        idx = system[i].find('=');
        size = system[i].size();
        tmp = system[i].substr(0, idx);
        this->X0 = stringFix(tmp);
        if (this->X0 != X0) {
            throw std::logic_error("getSysInfo: y(X0) != X0");
        }
        //task.a = FuncMaker(system[i].substr(idx + 1, size - idx), {}).func(0);
        Y.push_back(FuncMaker(system[i].substr(idx + 1, size - idx)).calculate());
    }
    this->X0 = X0;
    this->Xn = Xn;
}

void KoshiTask::setTaskInfo(const std::string &ode, uint64_t order, const std::vector<double> &Y0, double X0, double Xn) {
    ode_system.clear();
    Y.clear();
    uint64_t idx = 0, size = 0;
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
        ode_system.push_back(func);
    }
    idx = ode.find('=');
    size = ode.size();
    // FuncMaker y_order1(system[0].substr(0, idx), args), y_order2(system[0].substr(idx + 1, size - idx), args);
    // FuncMaker coeff = y_order1.getCoeff(order + 1);
    // auto func = [=] (const std::vector<double> &args) -> double {
    //     return (y_order2(args) - y_order1(args)) / coeff(args);
    // };
    FuncMaker y_order2(ode.substr(idx + 1, size - idx), args);
    //FuncMaker coeff = y_order1.getCoeff(order + 1);
    auto func = [=] (const std::vector<double> &args) -> double {
        return y_order2(args);
    };
    ode_system.push_back(func);
    Y = Y0;
    this->X0 = X0;
    this->Xn = Xn;
}

void KoshiTask::setSystemInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn) {
    ode_system.clear();
    Y.clear();
    uint64_t idx = 0, size = 0;
    std::vector<std::string> args; //x y y' y'' ...
    args.push_back("x");
    for (uint64_t i = 0; i < order; ++i) {
        uint64_t pos = system[i].find('\'');
        args.push_back(system[i].substr(0, pos));
        std::cout << "var: " << args.back() << "\n";
    }

    for (uint64_t i = 0; i < order; ++i) {
        idx = system[i].find('=');
        size = system[i].size();
        FuncMaker y_n(system[i].substr(idx + 1, size - idx), args);
        ode_system.push_back(y_n);
    }

    for (uint64_t i = order; i < order * 2; ++i) {
        idx = system[i].find('=');
        size = system[i].size();
        //tmp = system[i].substr(0, idx);
        this->X0 = stringFix(system[i].substr(0, idx));
        if (this->X0 != X0) {
            throw std::logic_error("getSysInfo: y(X0) != X0");
        }
        //task.a = FuncMaker(system[i].substr(idx + 1, size - idx), {}).func(0);
        Y.push_back(FuncMaker(system[i].substr(idx + 1, size - idx)).calculate());
    }
    this->X0 = X0;
    this->Xn = Xn;
    //std::cout << "size: " << ode_system.size() << "\n";
}

std::tuple<double, double> KoshiTask::getBorders () const {
    return std::make_tuple(X0, Xn);
}