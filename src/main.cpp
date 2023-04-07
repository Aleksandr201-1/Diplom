#include <iostream>
#include <iomanip>
//#include "gnuplot-iostream.h"
// #include "ODUSolver/RungeKutta.hpp"
// #include "ChemicalGenerator/RightPartGen.hpp"
// #include "PDFReporter/ReportGenerator.hpp"
#include <ODUSolver/RungeKutta.hpp>
#include <ChemicalGenerator/RightPartGen.hpp>
#include <PDFReporter/ReportGenerator.hpp>

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

//схема CROS
int main (int argc, char* argv[]) {
    ReportInfo info;
    ChemicalSystem sys;
    sys.initFromFile("./test/ChemicTest/bufermm.txt");
    sys.setPressure(100'000);
    sys.setTemperature(1000);
    //sys.setDensity(0.95);
    sys.addReaction("H2 + O2 => 2OH");
    sys.setReactionParameters(0, 1.7 * std::pow(10, 7), 0, 24044);
    sys.addReaction("H + O2 => OH + O");
    sys.setReactionParameters(1, 1.987 * std::pow(10, 8), 0, 8456);
    sys.addReaction("H2 + OH => H2O + H");
    sys.setReactionParameters(2, 1.024 * std::pow(10, 2), 1.6, 1660);
    sys.addReaction("H2 + O => OH + H");
    sys.setReactionParameters(3, 5.119 * std::pow(10, -2), 2.67, 3163);
    sys.addReaction("2OH => H2O + O");
    sys.setReactionParameters(4, 1.506 * std::pow(10, 3), 1.14, 50);
    sys.addReaction("H + OH => H2O");
    sys.setReactionParameters(5, 2.212 * std::pow(10, 16), -2.0, 0);
    //sys.addReaction("2H => H2");
    //sys.setReactionParameters(6, 9.791 * std::pow(10, 10), -0.6, 0);

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
        } else if (args[i] == "--system" || args[i] == "-s") {
            info.multigraph = true;
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

    if (info.multigraph) {
        std::cin >> order;
        for (uint64_t i = 0; i < order * 2; ++i) {
            info.input_task.push_back(readLine());
        }
    } else {
        info.input_task.push_back(readLine());
        order = getOrder(info.input_task[0]);
        for (uint64_t i = 0; i < order; ++i) {
            //system.push_back(readLine());
            info.input_task.push_back(readLine());
        }
    }

    std::cout << "Введите размер шага:\n";
    std::cin >> info.h;
    std::cout << "Введите границы:\n";
    std::cin >> X0 >> Xn;
    std::cout << "Введите функцию для сравнения:\n";
    uint64_t countOfAnalitic = info.multigraph ? order : 1;
    for (uint64_t i = 0; i < countOfAnalitic; ++i) {
        check.reset(readLine(), {"x"});
        info.analitic.push_back(check);
        std::cout << "analitic " << i + 1 << ": " << info.analitic[i] << "\n";
    }

    if (info.multigraph) {
        info.task = getSysInfo(info.input_task, order, X0, Xn);
    } else {
        info.task = getTaskInfo(info.input_task, order, X0, Xn);
    }
    // std::vector<std::function<double(const std::vector<double> &)>> sus = {
    //     [&] (const std::vector<double> &args) -> double {
    //         return -10000 * args[1] - 4999 * args[2];
    //     },
    //     [&] (const std::vector<double> &args) -> double {
    //         return -10001 * args[1] - 5000 * args[2];
    //     }
    // };

    // info.task.odu_system = sys.rightPartGen();
    // info.task.order = 6;
    // info.task.X0 = 0;
    // info.h = 0.01;
    // info.task.Xn = 3;
    // info.task.Y = {0.5, 0, 0, 0.5, 0, 0};

    // std::cout << "Chem size: " << info.task.odu_system.size() << "\n";
    // std::cout << "Chemic func exmpl: " << info.task.odu_system[0]({0, 0.5, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0}) << "\n";
    // return 0;

    // double ru = 10;
    // // auto f1 = [=] (const std::vector<double> &X) -> double {
    // //     return -(ru + 2)*X[1] + ru*X[2]*X[2];
    // // };
    // // auto f2 = [=] (const std::vector<double> &X) -> double {
    // //     return X[1] - X[2] - X[2]*X[2];
    // // };
    // auto f1 = [=] (const std::vector<double> &X) -> double {
    //     return -(ru + 2)*X[2] + ru*X[1]*X[1];
    // };
    // auto f2 = [=] (const std::vector<double> &X) -> double {
    //     return X[2] - X[1] - X[1]*X[1];
    // };
    // info.task.odu_system = {f2, f1};
    // info.task.order = 2;
    // info.task.X0 = 0;
    // info.task.Xn = 3;
    // info.task.Y = {1, 1};
    // info.analitic = FunctionalTree("exp(-x)", std::vector<std::string>{"x"});
    //return 0;

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