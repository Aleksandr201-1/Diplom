#include <iostream>
#include <iomanip>
#include <ODUSolver/Koshi/KoshiSolver.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <PDFReporter/ReportGenerator.hpp>
#include <Input/InputHandle.hpp>
#include <chrono>

using duration_t = std::chrono::milliseconds;

void printTable1 (const std::vector<std::vector<float128_t>> &Y, std::ostream &out) {
    for (uint64_t i = 0; i < Y[0].size(); ++i) {
        for (uint64_t j = 0; j < Y.size(); ++j) {
            out << Y[j][i] << "\t";
        }
        out << "\n";
    }
}

//схема CROS
int main (int argc, char* argv[]) {
    std::chrono::time_point <std::chrono::system_clock> startT, endT;
    uint64_t time = 0;
    ReportInfo info;
    std::string fileInput = "";

    std::vector<std::string> args(argv, argv + argc);
    // for (auto el : args) {
    //     std::cout << el << "\n";
    // }
    argsHandler(args, info);
    info.butcher = createButcherTable(info.method, info.order, info.way);

    std::cout << "Введите задачу ";
    uint64_t order = 0;
    ChemicalSystem sys;
    KoshiTask koshi;
    if (info.type == TaskType::KOSHI) {
        std::cout << "(Коши):\n";
        std::vector<std::string> system;
        //KoshiTask koshi;
        float128_t X0, Xn;
        FunctionalTree check;
        info.input_task.push_back(readLine());
        order = getOrder(info.input_task[0]);
        std::cout << "Порядок: " << order << "\n";
        for (uint64_t i = 0; i < order; ++i) {
            info.input_task.push_back(readLine());
        }
        std::cout << "Введите размер шага:\n";
        std::cin >> info.h;
        std::cout << "Введите границы интегрирования:\n";
        std::cin >> X0 >> Xn;
        std::cout << "Введите функцию для сравнения:\n";
        check.reset(readLine(), {"x"});
        info.analitic.push_back(check);
        koshi.setTaskInfo(info.input_task, order, X0, Xn);
        info.task = &koshi;
    } else if (info.type == TaskType::KOSHI_SYSTEM) {
        std::cout << "(Систему Коши):\n";
        std::vector<std::string> system;
        float128_t X0, Xn;
        FunctionalTree check;
        std::cin >> order;
        std::cout << "order: " << order << "\n";
        for (uint64_t i = 0; i < order * 2; ++i) {
            info.input_task.push_back(readLine());
        }
        std::cout << "Введите размер шага:\n";
        std::cin >> info.h;
        std::cout << "Введите границы интегрирования:\n";
        std::cin >> X0 >> Xn;
        std::cout << "Введите функции для сравнения:\n";
        for (uint64_t i = 0; i < order; ++i) {
            check.reset(readLine(), {"x"});
            info.analitic.push_back(check);
        }
        koshi.setSystemInfo(info.input_task, order, X0, Xn);
        info.task = &koshi;
    } else if (info.type == TaskType::CHEMICAL) {
        std::cout << "(Химической кинетики):\n";
        std::string filename, additive;
        float128_t T, P;
        uint64_t additiveCount;
        std::cout << "Введите название файла с данными о веществах:\n";
        std::cin >> filename;
        std::cout << "Введите температуру и давление:\n";
        std::cin >> T >> P;
        std::cout << "Введите количество добавок и добавки:\n";
        std::cin >> additiveCount;
        //std::cout << "filename: " << filename << "\n";
        sys.initFromFile(filename);
        sys.setTemperature(T);
        sys.setPressure(P);
        for (uint64_t i = 0; i < additiveCount; ++i) {
            std::vector<float128_t> add(sys.getSubstanceList().size());
            std::cin >> additive;
            for (uint64_t j = 0; j < add.size(); ++j) {
                std::cin >> add[j];
            }
            sys.addAdditive(additive, add);
        }
        sys.printInfo(std::cout);

        std::cout << "Введите количество реакций и реакции:\n";
        std::cin >> order;
        for (uint64_t i = 0; i < order; ++i) {
            float128_t A, n, E;
            std::string reaction = readLine();
            std::cin >> A >> n >> E;
            sys.addReaction(reaction, A, n, E);
        }

        std::cout << "Введите количество начальных концентраций и концентрации:\n";
        std::cin >> additiveCount;
        auto substances = sys.getSubstanceList();
        std::vector<float128_t> initGamma(substances.size(), 0);
        for (uint64_t i = 0; i < additiveCount; ++i) {
            float128_t concentration;
            std::string sub, tmp;
            std::cin >> sub >> tmp >> concentration;
            auto it = std::find(substances.cbegin(), substances.cend(), sub);
            uint64_t idx = std::distance(substances.cbegin(), it);
            if (idx == substances.size()) {
                std::cerr << "Некорректное название вещества\n";
                return 1;
            }
            initGamma[idx] = concentration;
        }
        std::cout << "Введите начальный размер шага:\n";
        std::cin >> info.h;
        sys.setConcentrations(initGamma);
        //sys.setConcentrations({0, 0, 0, 0, 0, 0});
        sys.rightPartGen();
        std::cout << "info:\n";
        sys.printInfo(std::cout);
        info.task = &sys;
    } else {
        std::cerr << "Некорректный тип задачи\n";
        return 1;
    }

    //chemic
    std::vector<std::string> colors = {
        "blue",
        "purple",
        //"red",
        "black",
        "magenta",
        "green",
        "orange"
    };

    std::vector<std::string> substances = sys.getSubstanceList();
    // info.task.odu_system = sys.rightPartGen();
    // info.task.order = info.task.odu_system.size();
    // info.task.X0 = 0;
    // info.task.Xn = 2e-6;
    // info.h = (info.task.Xn - info.task.X0) / 50;
    // info.task.Y = sys.getY0();
    // printVector(info.task.Y);
    switch (info.type) {
        case TaskType::CHEMICAL:
            for (uint64_t i = 0; i < substances.size(); ++i) {
                info.graph_info.push_back({colors[i], substances[i]});
            }
            info.table = sys.getTable();
            break;
        case TaskType::KOSHI:
        case TaskType::KOSHI_SYSTEM:
            for (uint64_t i = 0; i < order; ++i) {
                info.graph_info.push_back({colors[i], std::to_string(i + 1)});
            }
            break;
        default:
            break;
    }

    //info.tough_coeff = ToughCoeff(info.task);
    info.tough_coeff = 0;

    std::cout << "=========Расчёт========\n";
    startT = std::chrono::system_clock::now();
    switch (info.type) {
        case TaskType::KOSHI:
        case TaskType::KOSHI_SYSTEM:
            info.solution = KoshiSolver(info.method, *static_cast<KoshiTask*>(info.task), info.butcher, info.h, info.algo, info.approx);
            //info.solution = RungeKutta(koshi, info.butcher, info.h);
            break;
        case TaskType::CHEMICAL:
            //info.solution = KoshiSolver(info.method, *static_cast<KoshiTask*>(info.task), info.butcher, info.h, info.algo, info.approx);
            info.solution = ChemicalSolver(info.method, *static_cast<ChemicalSystem*>(info.task), info.butcher, info.h, info.algo, info.approx, ReactionType::ISOTERM_CONST_P);
            //info.solution = NonExpl(*static_cast<ChemicalSystem*>(info.task), info.butcher, info.h, info.algo, info.approx);
            break;
        default:
            exit(1);
    }
    //info.solution = KoshiSolver(info.method, *static_cast<KoshiTask*>(&info.task), info.butcher, info.h, info.algo, info.approx);
    endT = std::chrono::system_clock::now();
    std::cout << "=====Конец расчёта=====\n";
    time += std::chrono::duration_cast<duration_t>(endT - startT).count();
    std::cout << "Time: " << time << "\n";
    info.workTime = time;
    
    std::ofstream file("./report/report.tex");
    auto &out = file;
    generateReport(info, ReportType::TEX, out);
    file.close();
    //printTable1(info.solution, std::cout);

    return 0;
}