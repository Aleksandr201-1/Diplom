#include <iostream>
#include <iomanip>
#include <ODUSolver/Koshi/KoshiSolver.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <PDFReporter/ReportGenerator.hpp>
#include <Input/InputHandle.hpp>
#include <chrono>

using duration_t = std::chrono::milliseconds;

void printTable1 (const std::vector<std::vector<double>> &Y, std::ostream &out) {
    for (uint64_t i = 0; i < Y[0].size(); ++i) {
        for (uint64_t j = 0; j < Y.size(); ++j) {
            out << Y[j][i] << "\t";
        }
        out << "\n";
    }
}

int main (int argc, char* argv[]) {
    std::chrono::time_point <std::chrono::system_clock> startT, endT;
    uint64_t time = 0;
    ReportInfo info;
    std::stringstream str_buf;

    std::vector<std::string> args(argv, argv + argc);
    argsHandler(args, info);
    info.butcher = createButcherTable(info.method, info.order, info.way);

    if (info.fileInput != "") {
        std::ifstream fileInput(info.fileInput);
        std::string str;
        if (!fileInput.good()) {
            throw std::logic_error("main: cant open file with name \"" + info.fileInput + "\"");
        }
        while (!fileInput.eof()) {
            str = readLine(fileInput);
            str_buf << str << '\n';
        }
    }
    auto &in = (info.fileInput != "") ? str_buf : std::cin;
    std::ofstream file("./report/report.tex");
    auto &out = file;

    std::cout << "Введите задачу ";
    uint64_t order = 0;
    ChemicalSystem sys;
    KoshiTask koshi;
    if (info.type == TaskType::KOSHI) {
        std::cout << "(Коши):\n";
        std::vector<std::string> system;
        double X0, Xn;
        FuncMaker check;
        info.input_task.push_back(readLine(in));
        order = getOrder(info.input_task[0]);
        std::cout << "Порядок: " << order << "\n";
        for (uint64_t i = 0; i < order; ++i) {
            info.input_task.push_back(readLine(in));
        }
        std::cout << "Введите размер шага:\n";
        in >> info.h_min;
        std::cout << "Введите границы интегрирования:\n";
        in >> X0 >> Xn;
        std::cout << "Введите функцию для сравнения:\n";
        check.reset(readLine(in), {"x"});
        info.analitic.push_back(check);
        koshi.setTaskInfo(info.input_task, order, X0, Xn);
        info.task = &koshi;
    } else if (info.type == TaskType::KOSHI_SYSTEM) {
        std::cout << "(Систему Коши):\n";
        std::vector<std::string> system;
        double X0, Xn;
        FuncMaker check;
        in >> order;
        std::cout << "order: " << order << "\n";
        for (uint64_t i = 0; i < order * 2; ++i) {
            info.input_task.push_back(readLine(in));
        }
        std::cout << "Введите размер шага:\n";
        in >> info.h_min;
        std::cout << "Введите границы интегрирования:\n";
        in >> X0 >> Xn;
        std::cout << "Введите функции для сравнения:\n";
        for (uint64_t i = 0; i < order; ++i) {
            check.reset(readLine(in), {"x"});
            info.analitic.push_back(check);
        }
        koshi.setSystemInfo(info.input_task, order, X0, Xn);
        info.task = &koshi;
    } else if (info.type == TaskType::CHEMICAL) {
        std::cout << "(Химической кинетики):\n";
        std::string filename, additive;
        double T, P;
        uint64_t additiveCount;
        std::cout << "Введите название файла с данными о веществах:\n";
        in >> filename;
        std::cout << "Введите температуру и давление:\n";
        in >> T >> P;
        std::cout << "Введите количество добавок и добавки:\n";
        in >> additiveCount;
        //std::cout << "filename: " << filename << "\n";
        sys.initFromFile(filename);
        sys.setTemperature(T);
        sys.setPressure(P);
        for (uint64_t i = 0; i < additiveCount; ++i) {
            std::vector<double> add(sys.getSubstanceList().size());
            in >> additive;
            for (uint64_t j = 0; j < add.size(); ++j) {
                in >> add[j];
            }
            sys.addAdditive(additive, add);
        }
        sys.printInfo(std::cout);

        std::cout << "Введите количество реакций и реакции:\n";
        in >> order;
        for (uint64_t i = 0; i < order; ++i) {
            double A, n, E;
            std::string reaction = readLine(in);
            //std::getline(in, reaction);
            in >> A >> n >> E;
            sys.addReaction(reaction, A, n, E);
        }
        std::cout << "Введите режим ввода концентраций:\n";
        std::string mode;
        in >> mode;

        std::cout << "Введите количество начальных концентраций и концентрации:\n";
        in >> additiveCount;
        auto substances = sys.getSubstanceList();
        std::vector<double> initGamma(substances.size(), 0);
        for (uint64_t i = 0; i < additiveCount; ++i) {
            double concentration;
            std::string sub, tmp;
            in >> sub >> tmp >> concentration;
            auto it = std::find(substances.cbegin(), substances.cend(), sub);
            uint64_t idx = std::distance(substances.cbegin(), it);
            if (idx == substances.size()) {
                std::cerr << "Некорректное название вещества\n";
                return 1;
            }
            initGamma[idx] = concentration;
        }
        std::cout << "Введите начальный размер шага:\n";
        in >> info.h_min >> info.h_max >> info.h_last;
        if (mode == "percent") {
            sys.setConcentrations(initGamma, ConcentrationMode::PERCENT);
        } else if (mode == "molar") {
            sys.setConcentrations(initGamma, ConcentrationMode::MOLAR_MASS);
        }
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
        "red",
        "black",
        "magenta",
        "green",
        "orange",
        "cyan",
        "pink",
        "yellow",
        "blue"
    };

    std::vector<std::string> substances = sys.getSubstanceList();
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
            info.solution = KoshiSolver(info.method, *static_cast<KoshiTask*>(info.task), info.butcher, info.h_min, info.algo, info.approx);
            break;
        case TaskType::CHEMICAL:
            info.solution = ChemicalSolver(info.method, *static_cast<ChemicalSystem*>(info.task), info.butcher, info.h_min, info.h_max, info.h_last, info.algo, info.approx, ReactionType::ADIABAT_CONST_RHO);
            break;
        default:
            exit(1);
    }
    endT = std::chrono::system_clock::now();
    std::cout << "=====Конец расчёта=====\n";
    time += std::chrono::duration_cast<duration_t>(endT - startT).count();
    std::cout << "Time: " << time << "\n";
    info.workTime = time;
    
    generateReport(info, ReportType::TEX, out);
    file.close();
    //printTable1(info.solution, std::cout);

    return 0;
}