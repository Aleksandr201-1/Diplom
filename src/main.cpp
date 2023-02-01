#include <iostream>
#include <iomanip>
//#include <io.h>
//#include <fcntl.h>
//#include "gnuplot-iostream.h"
#include "NotTough.hpp"
#include "General/LSM.hpp"

// void plot (const std::vector<std::string> &func, double a, double b) {
//     static std::vector<std::string> colors = {"red", "green", "blue"};
//     Gnuplot gp;

//     gp << "set xlabel \"X\"\n";
//     gp << "set ylabel \"Y\"\n";
//     gp << "set xzeroaxis lw 1\n";
//     gp << "set yzeroaxis lw 1\n";
//     gp << "set xrange[" << a << ":" << b << "]\n";
//     gp << "set grid\n";
//     gp << "set title \"Plot\" font \"Helvetica Bold, 20\"\n";
    
//     gp << "plot ";
//     for (uint64_t i = 0; i < func.size(); ++i) {
//         gp << func[i] << " title \"";
//         gp << (i == 0 ? "our func" : "real func");
//         gp << "\"  lc rgb \"" << colors[i % colors.size()] << "\"";
//         if (i == func.size() - 1) {
//             gp << "\n";
//         } else {
//             gp << ",";
//         }
//     }
// }

std::vector<double> getY (const FunctionalTree &func, const std::vector<double> &X) {
    std::vector<double> Y;
    for (uint64_t i = 0; i < X.size(); ++i) {
        Y.push_back(func(X[i]));
    }
    return Y;
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

void printTable (const std::vector<double> &X, const std::vector<double> &Yn, const std::vector<double> &Ya) {
    const uint64_t WIDTH = 12;
    std::cout << "Таблица:\n";
    std::cout << " ____________________________________________\n";
    std::cout << "| " << std::setw(WIDTH) << "X";
    std::cout << " | " << std::setw(WIDTH) << "Y nu";
    std::cout << " | " << std::setw(WIDTH) << "Y an";
    std::cout << " |\n";
    std::cout << "|______________|______________|______________|\n";
    for (uint64_t i = 0; i < X.size(); ++i) {
        std::cout << "| " << std::setw(WIDTH) << X[i];
        std::cout << " | " << std::setw(WIDTH) << Yn[i];
        std::cout << " | " << std::setw(WIDTH) << Ya[i];
        std::cout << " |\n";
        std::cout << "|______________|______________|______________|\n";
    }
}

int main () {
    //_setmode(_fileno(stdout), _O_U16TEXT);
    //_setmode(_fileno(stdin),  _O_U16TEXT);
    //_setmode(_fileno(stderr), _O_U16TEXT);

    Matrix<double> b;
    b = createButcherTable(SolveMethod::FALBERG, 2, 1);
    std::cout << "Таблица Бутчера " << b.size().n << "x" << b.size().m << ":\n" << b << "\n";
    std::cout << "Введите задачу:\n";
    std::vector<std::string> system;
    uint64_t order = 0;
    std::string str;
    double h, X0, Xn;
    std::pair<std::vector<double>, std::vector<double>> res1, res2;
    Task task;
    FunctionalTree check;

    system.push_back(readLine());
    order = getOrder(system[0]);
    std::cout << "Порядок задачи: " << order << "\n";
    for (uint64_t i = 0; i < order; ++i) {
        //str = readLine();
        system.push_back(readLine());
    }
    std::cout << "Введите размер шага:\n";
    std::cin >> h;
    std::cout << "Введите границы:\n";
    std::cin >> X0 >> Xn;
    std::cout << "Введите функцию для сравнения:\n";
    str = readLine();
    check.reset(str, {"x"});
    task = getTaskInfo(system, order, X0, Xn);
    std::cout << "Func test: " << task.f[0]({1,2,3,4}) << "\n";
    std::cout << "Tree test: " << task.trees[0]({1,2,3,4}) << "\n";

    std::vector<std::function<double(const std::vector<double> &)>> sus = {
        // [&] (const std::vector<double> &args) -> double {
        //     return args[2];
        // },
        // [&] (const std::vector<double> &args) -> double {
        //     static auto c2 = task.trees[0].getCoeff(order + 1);
        //     return -task.trees[0](args) / c2(args[0]);
        // }
        [&] (const std::vector<double> &args) -> double {
            return -10000 * args[1] - 4999 * args[2];
        },
        [&] (const std::vector<double> &args) -> double {
            return -10001 * args[1] - 5000 * args[2];
        }
    };
    std::cout << "Жёсткость: " << ToughCoeff(sus, task) << "\n";

    std::cout << "=====Рунге=====\n";
    // res1 = RungeKutta4(task, h);
    // res2 = RungeKutta4(task, 2 * h);
    //res2 = RungeKutta4(task, h);
    std::cout << "rk\n";
    res1 = RungeKutta6(task, b, h);
    //res1 = Falberg(task, b, h);
    //res1 = NonExpl(task, b, h);
    std::cout << "X: ";
    printVector(res1.first);
    std::cout << "Y numeric1: ";
    printVector(res1.second);
    //std::cout << "Y numeric2: ";
    //printVector(res2.second);
    std::cout << "Y analitic: ";
    auto analitic = getY(check, res1.first);
    printVector(analitic);

    printTable(res1.first, res1.second, analitic);

    auto func = LeastSquareMethod(res1.first, res1.second, 3);
    std::cout << LSMToText(func) << "\n";
    //plot({LSMToText(func), str}, res1.first[0], res1.first.back());
    //plot({LSMToText(func), str}, res1.first[0], res1.first.back());

    return 0;
}