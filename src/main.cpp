#include <iostream>
//#include <io.h>
//#include <fcntl.h>
//#include "gnuplot-iostream.h"
#include "NotTough.hpp"

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

int main () {
    //_setmode(_fileno(stdout), _O_U16TEXT);
    //_setmode(_fileno(stdin),  _O_U16TEXT);
    //_setmode(_fileno(stderr), _O_U16TEXT);

    Matrix<double> b;
    b = createButcherTable(4, 1);
    std::cout << "Таблица Бутчера:\n" << b << "\n";
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
    task = getTaskInfo(system, order, X0, Xn);

    std::cout << "=====Рунге=====\n";
    // res1 = RungeKutta4(task, h);
    // res2 = RungeKutta4(task, 2 * h);
    res1 = RungeKutta6(task, b, h);
    res2 = RungeKutta6(task, b, 2 * h);
    std::cout << "X: ";
    printVector(res1.first);
    std::cout << "Y: ";
    printVector(res1.second);
    //std::cout << "Погрешность: " << RungeRomberg(res1.second[2], res2.second[1]) << "\n";

    auto func = LeastSquareMethod(res1.first, res1.second, 3);
    std::cout << LSMToText(func) << "\n";
    //plot({LSMToText(func), str}, res1.first[0], res1.first.back());
    //plot({LSMToText(func), str}, res1.first[0], res1.first.back());

    return 0;
}