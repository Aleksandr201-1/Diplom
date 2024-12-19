#include <iostream>
#include <iomanip>
#include <ODUSolver/Koshi/KoshiSolver.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <PDFReporter/ReportGenerator.hpp>
#include <Input/InputHandle.hpp>
#include <GUI/Plot2D.hpp>
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

    std::vector<std::string> args;
    std::string oneArg;
    std::ifstream fileArgs("init.txt");
    while (fileArgs >> oneArg) {
        args.push_back(oneArg);
    }
    fileArgs.close();
    std::cout << "Args:\n";
    for (auto el : args) {
        std::cout << el << "|\n";
    }
    std::cout << "\n";

    //std::vector<std::string> args(argv, argv + argc);
    argsHandler(args, info);
    info.type = TaskType::CHEMICAL;
    // info.fileInput = "../../test/ChemicTest/test2-1.txt";
    // info.order = 4;
    // info.way = 1;
    // info.method = SolveMethod::RUNGE_KUTTA;
    // info.approx = 0.001;
    // info.algo = IterationAlgo::ZEIDEL;
    info.butcher = createButcherTable(info.method, info.order, info.way);
    std::cout << "Butcher matrix:\n" << info.butcher << "\n";

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
        //sys.printInfo(std::cout);

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
        std::string mode, type;
        in >> mode >> type;
        info.react = stringToReactionType(type);
        std::cout << "Тип реакции: " << reactionTypeToString(info.react) << "\n";

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
        std::cout << "Введите начальный размер шага и конечное значение по времени:\n";
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
    uint64_t dropNum;
    in >> dropNum;
    std::cout << "Размер выходных сеток: " << dropNum << "\n";

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

    std::cout << "========Расчёт========\n";
    startT = std::chrono::system_clock::now();
    switch (info.type) {
        case TaskType::KOSHI:
        case TaskType::KOSHI_SYSTEM:
            info.solution = KoshiSolver(info.method, *static_cast<KoshiTask*>(info.task), info.butcher, info.h_min, info.algo, info.approx);
            break;
        case TaskType::CHEMICAL:
            info.solution = ChemicalSolver(info.method, *static_cast<ChemicalSystem*>(info.task), info.butcher, info.h_min, info.h_max, info.h_last, info.algo, info.approx, info.react);
            break;
        default:
            exit(1);
    }
    //info.tough_coeff = ToughCoeff(info.task->getODE(), info.task->getY0());
    //std::cout << "Tough coeff: " << info.tough_coeff << "\n";
    endT = std::chrono::system_clock::now();
    std::cout << "=====Конец расчёта=====\n";
    time += std::chrono::duration_cast<duration_t>(endT - startT).count();
    std::cout << "Time: " << time << "\n";
    info.workTime = time;
    //info.solution = drop(info.solution, dropNum);
    //generateReport(info, ReportType::TEX, out);
    file.close();
    //printTable1(info.solution, std::cout);
    //exit(0);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 6;
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Plot for... idk? maybe U?", 7U, settings);
    window.setFramerateLimit(20);
    uint64_t solutionSize = 11;//info.solution.size();
    std::ifstream chemLog("chem_log.txt"), toughOut("tough_output.txt");
    std::string str;
    std::getline(chemLog, str);
    std::vector<std::vector<float>> solution(solutionSize);
    std::vector<float> toughVec;
    while (!chemLog.eof()) {
        int step;
        chemLog >> step;
        for (uint64_t j = 0; j < solutionSize; ++j) {
            float num;
            chemLog >> num;
            std::cout << num << "0";
            solution[j].push_back(num);
        }
        std::cout << "\n";
    }
    while (!toughOut.eof()) {
        float num;
        toughOut >> num;
        toughVec.push_back(num);
    }
    chemLog.close();
    toughOut.close();
    for (uint64_t j = 0; j < solutionSize; ++j) {
        //solution[j].erase(solution[j].begin(), solution[j].begin() + 550);
    }
    //toughVec.erase(toughVec.begin(), toughVec.begin() + 550);
    std::cout << "\n";
    std::cin >> str;
    Plot2D plot, plotTemp, plotRho, plotTough;
    plot.setYName(L"nu/m");
    plot.setXName(L"t");
    plot.setPosition(0, 0);
    for (uint64_t i = solutionSize - 4; i < solutionSize - 2; ++i) {
        plot.addFunc(solution[0], solution[i], Line(Line::Type::LINE, sf::Color::Black, 2.f));
    }
    auto list = sys.getSubstanceList();
    plot.setColors({sf::Color::Red, sf::Color::Green, sf::Color::Blue, sf::Color::Black, sf::Color::Yellow, sf::Color::Magenta, sf::Color::Cyan, sf::Color(124, 145, 243), sf::Color(145, 123, 243), sf::Color(243, 145, 132), sf::Color(154, 234, 128)});
    plot.setNames(std::vector<sf::String>(list.begin(), list.end()));
    auto plotLeg = plot.generateLegend();
    plotLeg->setPosition({plot.getSize().x, 0.f});
    plotLeg->update();

    plotTemp.setYName(L"T, K");
    plotTemp.setXName(L"t, s");
    plotTemp.addFunc(solution[0], solution[solutionSize - 2], Line(Line::Type::LINE, sf::Color::Black, 2.f));
    plotTemp.setColor(0, sf::Color::Black);
    plotTemp.setName(0, sf::String(L"Плотность"));
    plotTemp.setPosition({plot.getSize().x + plotLeg->getSize().x, 0.f});
    auto tempLeg = plotTemp.generateLegend();
    tempLeg->setPosition({plot.getSize().x + plotLeg->getSize().x + plotTemp.getSize().x, 0.f});
    tempLeg->update();

    plotRho.setYName(L"rho");
    plotRho.setXName(L"t, s");
    plotRho.addFunc(solution[0], solution[solutionSize - 1], Line(Line::Type::LINE, sf::Color::Black, 2.f));
    plotRho.setColor(0, sf::Color::Black);
    plotRho.setName(0, sf::String(L"Температура"));
    plotRho.setPosition({0.f, plot.getSize().y});
    auto rhoLeg = plotRho.generateLegend();
    rhoLeg->setPosition({plotRho.getSize().x, plot.getSize().y});
    rhoLeg->update();

    plotTough.setYName(L"");
    plotTough.setXName(L"t, s");
    plotTough.addFunc(solution[0], toughVec, Line(Line::Type::LINE, sf::Color::Black, 2.f));
    plotTough.setColor(0, sf::Color::Black);
    plotTough.setName(0, sf::String(L"Жёсткость"));
    plotTough.setPosition({plot.getSize().x + plotLeg->getSize().x, plot.getSize().y});
    auto toughLeg = plotTough.generateLegend();
    toughLeg->setPosition({plot.getSize().x + plotLeg->getSize().x + plotTough.getSize().x, plot.getSize().y});
    toughLeg->update();

    // plot2.setColors({sf::Color::Black, sf::Color::Green, sf::Color::Blue});
    // //plot2.addLegends({"Plot for " + argToStr(toDraw) + ", 1", "Plot for " + argToStr(toDraw) + ", " + std::to_string(n / 2), "Plot for " + argToStr(toDraw) + ", " + std::to_string(n)});
    // plot2.setXName(L"x/R0");
    // plot2.setYName(L"f");
    // plot2.setPosition(0, 550);

    while (window.isOpen()) {
        sf::Event event;

        while (window.pollEvent(event)) {
            sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
            sf::Vector2f mousePos = window.mapPixelToCoords(pixelPos);
            if (event.type == sf::Event::Closed || sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
                window.close();
            }
            plot.handleEvent(event, mousePos);
            plotTemp.handleEvent(event, mousePos);
            plotRho.handleEvent(event, mousePos);
            plotTough.handleEvent(event, mousePos);
        }
        plot.update();
        plotTemp.update();
        plotRho.update();
        plotTough.update();

        window.clear(sf::Color(127, 127, 127));
        window.draw(plot);
        window.draw(plotTemp);
        window.draw(plotRho);
        window.draw(plotTough);
        window.draw(*plotLeg);
        window.draw(*tempLeg);
        window.draw(*rhoLeg);
        window.draw(*toughLeg);
        window.display();
    }

    return 0;
}