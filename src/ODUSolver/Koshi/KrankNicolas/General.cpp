#include "General.hpp"

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

std::string toString (double val, uint64_t precision) {
    return std::to_string(val).substr(0, std::to_string(val).find(".") + precision + 1);
}

double stringFix (std::string &str) {
    double val = 0;
    std::string tmp;
    std::swap(tmp, str);
    for (uint64_t i = 0; i < tmp.size(); ++i) {
        if (tmp[i] == 'u') {
            std::string valStr;
            while (tmp[i] != '(') {
                str += tmp[i];
                ++i;
            }
            ++i;
            while (tmp[i] != ')') {
                valStr += tmp[i];
                ++i;
            }
            val = std::atof(valStr.c_str());
        } else {
            str += tmp[i];
        }
    }
    return val;
}

Task getTaskInfo(const std::vector<std::string> &system) {
    uint64_t idx = 0, size = 0;
    std::string tmp;
    Task task;

    task.X.resize(2);
    task.coeff = std::vector<double>(4, 0);

    idx = system[0].find('=');
    size = system[0].size();
    task.trees.push_back(std::move(FunctionalTree(system[0].substr(idx + 1, size - idx), {"uxx", "ux", "u", "x", "t"}))); //записываем правую часть уравнения
    task.coeff[0] = task.trees[0].getCoeff(0).func(0);//коэффициент при uxx
    task.coeff[1] = task.trees[0].getCoeff(1).func(0);//коэффициент при ux
    task.coeff[2] = task.trees[0].getCoeff(2).func(0);//коэффициент при u

    for (uint64_t i = 1; i < 3; ++i) { //записываем краевые условия
        FunctionalTree tmpTree;
        idx = system[i].find('=');
        size = system[i].size();
        tmp = system[i].substr(0, idx);
        task.X[i - 1] = stringFix(tmp);

        tmpTree.reset(tmp, {"ux", "u"});
        task.ux.push_back({tmpTree.getCoeff(0).func(0), tmpTree.getCoeff(1).func(0)});
        task.trees.push_back(std::move(FunctionalTree(system[i].substr(idx + 1, size - idx), {"t"})));
    }
    //записываем условие в начальный момент времени
    idx = system[3].find('=');
    size = system[3].size();
    tmp = system[3].substr(0, idx);
    task.trees.push_back(std::move(FunctionalTree(system[3].substr(idx + 1, size - idx), {"x"})));

    return task;
}


//функция чтения строки целиком и игнорирования строк, начинающихся с '#'
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