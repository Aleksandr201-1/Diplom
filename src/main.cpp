#include <iostream>
#include <iomanip>
//#include <io.h>
//#include <fcntl.h>
//#include "gnuplot-iostream.h"
#include "NotTough.hpp"
//#include "General/LSM.hpp"
#include "RightPartGen.hpp"
#include "PDF-reporter/ReportGenerator.hpp"

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

// std::vector<double> getY (const FunctionalTree &func, const std::vector<double> &X) {
//     std::vector<double> Y;
//     for (uint64_t i = 0; i < X.size(); ++i) {
//         Y.push_back(func(X[i]));
//     }
//     return Y;
// }

// void printTable (const std::vector<double> &X, const std::vector<double> &Yn, const std::vector<double> &Ya) {
//     const uint64_t WIDTH = 12;
//     std::cout << "Таблица:\n";
//     std::cout << " ____________________________________________\n";
//     std::cout << "| " << std::setw(WIDTH) << "X";
//     std::cout << " | " << std::setw(WIDTH) << "Y nu";
//     std::cout << " | " << std::setw(WIDTH) << "Y an";
//     std::cout << " |\n";
//     std::cout << "|______________|______________|______________|\n";
//     for (uint64_t i = 0; i < X.size(); ++i) {
//         std::cout << "| " << std::setw(WIDTH) << X[i];
//         std::cout << " | " << std::setw(WIDTH) << Yn[i];
//         std::cout << " | " << std::setw(WIDTH) << Ya[i];
//         std::cout << " |\n";
//         std::cout << "|______________|______________|______________|\n";
//     }
// }

void help (const std::string &name) {
    std:: cout << "Usage: " << name << " [KEYS] [OPTIONS]\n"
                  "\t-m, --method\t[METHOD]\tРешать задачу методом [METHOD]\n"
                  "\t\tПоддерживаемые методы:\n\t\tRunge-Kutta, Cheskino, Dorman-Prince, Falberg, Gauss, Lobatto, L Stable Diagonal, Merson, Rado\n"
                  "\t-o, --order\t[ORDER]\t\tРешать задачу методом [METHOD] с порядком точности [ORDER]\n"
                  "\t-w, --way\t[WAY]\t\tРешать задачу методом [METHOD] с порядком точности [ORDER] и способом реализации [WAY]\n"
                  "\t-i, --iter\t[ALGORYTHM]\tИспользовать алгоритм [ALGORYTHM] для решения системы уравнений в жёстких методах\n"
                  "\t\tПоддерживаемые алгоритмы:\n\t\tSI, Zeidel, Newton\n"
                  "\t-a, --approx\t[NUM]\t\tПродолжать процесс итераций в жёстких методах пока не будет достигнута точность [NUM]\n";
}

int main (int argc, char* argv[]) {
    //ChemicalReaction::setAtomList({"H", "O"});
    //ChemicalReaction::setSubstanceList({"H2", "O2", "H2O", "OH", "H", "O"});
    //ChemicalReaction("H + OH => H2O");
    ReportInfo info;
    //uint64_t methodPos = 0;
    std::vector<std::string> args(argv, argv + argc);
    for (auto el : args) {
        std::cout << el << "\n";
    }
    for (uint64_t i = 1; i < args.size(); ++i) {
        if (args[i] == "--help" || args[i] == "-h") {
            help(args[0]);
            exit(0);
        } else if (args[i] == "--method" || args[i] == "-m") {
            if (args.size() > i + 1) {
                info.method = stringToSolveMethod(args[i + 1]);
                //methodPos = i + 1;
                ++i;
            }
        } else if (args[i] == "--order" || args[i] == "-o") {
            if (args.size() > i + 1) {
                info.order = std::stoull(args[i + 1]);
                ++i;
            }
        } else if (args[i] == "--way" || args[i] == "-w") {
            if (args.size() > i + 1) {
                info.way = std::stoull(args[i + 1]);
                ++i;
            }
        } else if (args[i] == "--iter" || args[i] == "-i") {
            if (args.size() > i + 1) {
                //info.way = std::stoull(args[i + 1]);
                if (args[i + 1] == "Newton") {
                    info.algo = IterationAlgo::NEWTON;
                } else if (args[i + 1] == "Zeidel") {
                    info.algo = IterationAlgo::ZEIDEL;
                } else if (args[i + 1] == "SI") {
                    info.algo = IterationAlgo::SIMPLE_ITERATION;
                }
                ++i;
            }
        } else if (args[i] == "--approx" || args[i] == "-a") {
            if (args.size() > i + 1) {
                info.approx = std::stod(args[i + 1]);
                ++i;
            }
        }
    }
    // info.method = SolveMethod::RUNGE_KUTTA;
    // info.order = 4;
    // info.way = 1;
    // info.method = SolveMethod::FALBERG;
    // info.order = 2;
    // info.way = 2;
    // info.method = SolveMethod::GAUSS;
    // info.order = 4;
    // info.way = 1;
    info.butcher = createButcherTable(info.method, info.order, info.way);

    std::cout << "Введите задачу:\n";
    std::vector<std::string> system;
    uint64_t order = 0;
    double X0, Xn;
    Task task;
    FunctionalTree check;

    info.input_task.push_back(readLine());
    order = getOrder(info.input_task[0]);
    for (uint64_t i = 0; i < order; ++i) {
        //system.push_back(readLine());
        info.input_task.push_back(readLine());
    }

    std::cout << "Введите размер шага:\n";
    std::cin >> info.h;
    std::cout << "Введите границы:\n";
    std::cin >> X0 >> Xn;
    std::cout << "Введите функцию для сравнения:\n";
    check.reset(readLine(), {"x"});
    info.task = getTaskInfo(info.input_task, order, X0, Xn);
    info.analitic = check;
    std::cout << "analitic: " << info.analitic << "\n";

    // std::vector<std::function<double(const std::vector<double> &)>> sus = {
    //     [&] (const std::vector<double> &args) -> double {
    //         return -10000 * args[1] - 4999 * args[2];
    //     },
    //     [&] (const std::vector<double> &args) -> double {
    //         return -10001 * args[1] - 5000 * args[2];
    //     }
    // };
    info.tough_coeff = ToughCoeff(info.task);

    std::cout << "=====Рунге=====\n";
    //info.solution = RungeKutta6(info.task, info.butcher, info.h);
    //info.solution = Falberg(info.task, info.butcher, info.h);
    //info.solution = NonExplZeidel(info.task, info.butcher, info.h);
    info.solution = ODESolve(info.method, info.task, info.butcher, info.h, info.algo, info.approx);

    
    std::ofstream file("./report/report.tex");
    //auto &out = std::cout;
    auto &out = file;
    generateReport(info, ReportType::TEX, out);
    file.close();
    //plot({LSMToText(func), str}, res1.first[0], res1.first.back());

    return 0;
}