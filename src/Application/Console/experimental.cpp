//#include <ODUSolver/CrankNicolson/CrankNicolson.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <ODUSolver/CrankNicolson/LameSolver.hpp>
#include <Input/InputHandle.hpp>
#include <GUI/Plot2D.hpp>
#include <GUI/Plot3D.hpp>
#include <GUI/HeatMap.hpp>
#include <GUI/Selector.hpp>
#include <GUI/TextField.hpp>
#include <GUI/PageHolder.hpp>
#include <GUI/ScrollingWindow.hpp>
#include <GUI/Box.hpp>
#include <NumericMethods/Differentiation.hpp>

using duration_t = std::chrono::milliseconds;

void getContextWindow (PageHolder &holder, uint64_t wWidth, uint64_t wHeight) {
    auto window = std::make_shared<ScrollingWindow>();
    auto mainBox = std::make_shared<Box>(40, Box::Alignment::LEFT, Box::Orientation::VERTICAL);
    auto profileChoiceBox = std::make_shared<Box>(20, Box::Alignment::LEFT, Box::Orientation::VERTICAL);
    auto sectionModBox = std::make_shared<Box>(20, Box::Alignment::LEFT, Box::Orientation::VERTICAL);
    auto saveModBox = std::make_shared<Box>(20, Box::Alignment::LEFT, Box::Orientation::VERTICAL);


}

int main (int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.precision(10);
    std::cout.setf(std::ios_base::fixed);

    ReportInfo info;
    uint64_t order = 0;
    std::vector<std::string> args(argv, argv + argc);
    
    std::stringstream str_buf;
    info.type = TaskType::CHEMICAL;
    //info.fileInput = "../../test/ChemicTest/test4.txt";
    info.fileInput = "./test7.txt";
    info.order = 4;
    info.way = 1;
    info.algo = IterationAlgo::ZEIDEL;
    info.method = SolveMethod::RUNGE_KUTTA;
    info.butcher = createButcherTable(info.method, info.order, info.way);
    info.approx = 0.01;
    argsHandler(args, info);
    auto time = std::chrono::system_clock::now();
    std::time_t date = std::chrono::system_clock::to_time_t(time);
    std::string date_str(std::ctime(&date));
    date_str.pop_back();
    std::cout << "Data: " << date_str << "\n";

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

    ChemicalSystem sys;
    std::cout << "(Химической кинетики):\n";
    std::string filename, additive;
    double Temp, Press;
    uint64_t additiveCount, n, m, initProfileCount, profileSize;
    double R0, x_h, xi_h;
    double nu = 1.0;
    uint64_t calcEndX, calcEndY;
    std::cout << "Введите размеры сетки:\n";
    in >> R0 >> calcEndX >> calcEndY >> n >> m;
    x_h = R0 * calcEndX / n;
    xi_h = R0 * calcEndY / m;
    in >> profileSize;
    std::cout << "x_h = " << x_h << "\nxi_h = " << xi_h << "\n";
    std::cout << "Введите название файла с данными о веществах:\n";
    in >> filename;
    sys.initFromFile(filename);
    std::vector<ValueProfile> valProf(ARGS::TOTAL_COUNT - 1, ValueProfile(profileSize));
    std::vector<ChemicalProfile> chemProf(sys.getSubstanceList().size(), ChemicalProfile(profileSize));
    in >> initProfileCount;
    for (uint64_t i = 0; i < initProfileCount; ++i) {
        std::string valName;
        in >> valName;
        double inV;
        in >> inV;
        valProf[strToArgs(valName)].in = inV;
    }
    in >> initProfileCount;
    std::vector<std::string> valueNames(initProfileCount);
    for (uint64_t i = 0; i < initProfileCount; ++i) {
        in >> valueNames[i];
    }
    for (uint64_t i = 0; i < valProf.size(); ++i) {
        for (uint64_t j = 0; j < profileSize; ++j) {
            valProf[i].profile[j].first = j * xi_h / (profileSize - 1) * m;
        }
    }
    for (uint64_t i = 0; i < profileSize; ++i) {
        for (uint64_t j = 0; j < initProfileCount; ++j) {
            double val;
            in >> val;
            valProf[strToArgs(valueNames[j])].profile[i].second = val;
        }
    }
    valueNames.resize(sys.getSubstanceList().size());
    for (uint64_t i = 0; i < sys.getSubstanceList().size(); ++i) {
        in >> valueNames[i];
    }
    for (uint64_t i = 0; i < profileSize; ++i) {
        for (uint64_t j = 0; j < chemProf.size(); ++j) {
            double val;
            in >> val;
            chemProf[j].profile[i].second = val;
            chemProf[j].profile[i].first = i * xi_h / (profileSize - 1) * m;
        }
    }
    in >> usingTurb;
    std::cout << "Turbulence mode: " << usingTurb << "\n";
    std::cout << "Введите количество добавок и добавки:\n";
    in >> additiveCount;
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
        in >> A >> n >> E;
        sys.addReaction(reaction, A, n, E);
    }
    std::cout << "Введите режим ввода концентраций:\n";
    std::string mode;
    in >> mode;

    std::cout << "Введите количество начальных концентраций и концентрации:\n";
    in >> additiveCount;
    auto substances = sys.getSubstanceList();
    std::vector<double> insideGamma(substances.size(), 0), outsideGamma(substances.size(), 0);
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
        insideGamma[idx] = concentration * 1000;
    }
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
        outsideGamma[idx] = concentration * 1000;
    }
    double sum1 = 0, sum2 = 0;
    double gammaSum = 0;
    for (int i = 0; i < insideGamma.size(); ++i) {
        sum1 += insideGamma[i] * sys.getSubstanceMasses().find(sys.getSubstanceList()[i])->second;
        sum2 += outsideGamma[i] * sys.getSubstanceMasses().find(sys.getSubstanceList()[i])->second;
        gammaSum += insideGamma[i];
    }
    //нормироквка (отключена) 
    for (int i = 0; i < insideGamma.size(); ++i) {
        //insideGamma[i]  /= sum1;
        //outsideGamma[i] /= sum2;
    }
    std::cout << "Init gamma sum: " << sum1 << " " << sum2 << " " << 1.0/gammaSum << "\n";
    sum1 = sum2 = gammaSum = 0;
    for (int i = 0; i < insideGamma.size(); ++i) {
        sum1 += insideGamma[i] * sys.getSubstanceMasses().find(sys.getSubstanceList()[i])->second;
        sum2 += outsideGamma[i] * sys.getSubstanceMasses().find(sys.getSubstanceList()[i])->second;
        gammaSum += insideGamma[i];
    }
    std::cout << "Init gamma sum: " << sum1 << " " << sum2 << " " << 1.0/gammaSum << "\n";
    for (uint64_t i = 0; i < sys.getSubstanceMasses().size(); ++i) {
        chemProf[i].in  = insideGamma[i];
        chemProf[i].out = outsideGamma[i];
    }
    std::cout << "Введите начальный размер шага:\n";
    in >> info.h_min >> info.h_max >> info.h_last;
    sys.setPressure(valProf[P].in);
    sys.setTemperature(valProf[T].in);
    sys.rightPartGen();
    std::cout << "info:\n";
    auto ental = sys.getEnthalpy();
    info.task = &sys;

    std::cout << "Введите название выходного файла:\n";
    std::string outputFile;
    in >> outputFile;
    int choice = 0;
    in >> choice;

    std::vector<Matrix<double>> f(TOTAL_COUNT + sys.getSubstanceMasses().size() - 1, Matrix<double>(n, m));
    valProf[V].in = valProf[W].in = valProf[U].in;
    //переход к физическим величинам
    for (uint64_t i = 0; i < C1; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](0, j) = valProf[i].getPhysical((double)j * xi_h);
        }
    }
    for (uint64_t i = C1; i < f.size(); ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](0, j) = chemProf[i - C1].getPhysical((double)j * xi_h);
        }
    }
    double Xi_0 = 0;
    std::vector<Matrix<double>> tmpF = f;

    //RHO
    auto pointGamma = insideGamma;
    for (uint64_t i = 0; i < m; ++i) {
        double muSum = 0;
        std::cout << "press for rho " << i << ":\n";
        for (uint64_t j = 0; j < pointGamma.size(); ++j) {
            pointGamma[j] = f[C1 + j](0, i);
            muSum += pointGamma[j];
            std::cout << f[C1 + j](0, i) << " "; 
        }

        sys.setTemperature(f[T](0, i));
        sys.setPressure(f[P](0, i));
        sys.setConcentrations(pointGamma, ConcentrationMode::MOLAR_MASS);
        std::cout << "T: " << f[T](0, i) << " P: " << f[P](0, i) << " muSum: " << muSum << " rho: " << sys.getRho() << "\n";
        f[RHO](0, i) = f[P](0, i) / 8.3144 / f[T](0, i) / muSum;
    }
    valProf[RHO].in = f[RHO](0, 0);
    valProf[P].in = valProf[U].in * valProf[U].in * valProf[RHO].in;
    //MU
    for (uint64_t j = 0; j < m; ++j) {
        f[MU](0, j) = turbulence[usingTurb](f, valProf, chemProf, R0, 0, j);
    }
    //J
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[J](i, j) = 0;
            for (uint64_t k = 0; k < ental.size(); ++k) {
                f[J](i, j) += ental[k](f[T](i, j)) * f[C1 + k](i, j);
            }
            f[J](i, j) += 1.0 / 2.0 * (f[U](i, j) * f[U](i, j) + f[V](i, j) * f[V](i, j) + f[W](i, j) * f[W](i, j));
        }
    }
    valProf[J].in = valProf[U].in * valProf[U].in;
    valProf[Q].in = 1;
    valProf[V].in = valProf[W].in = valProf[U].in;

    std::cout << "physical values:\n";
    for (uint64_t i = 0; i < f.size(); ++i) {
        std::string str = argToStr(ARGS(i));
        if (str == "Error") {
            str = "C" + std::to_string(i - C1 + 1);
        }
        std::cout << "physical " << str << ":\n";
        for (uint64_t j = 0; j < m; ++j) {
            std::cout << f[i](0, j) << " ";
        }
        std::cout << "\n";
    }
    //заполнение доп. слоя
    for (uint64_t i = 0; i < f.size(); ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](1, j) = f[i](0, j);
        }
    }

    //нормировка (отключена)
    // for (uint64_t i = 0; i < C1; ++i) {
    //     for (uint64_t j = 0; j < m; ++j) {
    //         f[i](0, j) = valProf[i].toNormal(f[i](0, j));
    //     }
    // }
    // for (uint64_t i = C1; i < f.size(); ++i) {
    //     for (uint64_t j = 0; j < m; ++j) {
    //         f[i](0, j) = chemProf[i - C1].toNormal(f[i](0, j));
    //     }
    // }
    // for (uint64_t j = 0; j < m; ++j) {
    //     f[MU](0, j) *= valProf[MU].in;
    // }

    std::cout << "initial profiles:\n";
    //for vals
    for (uint64_t i = 0; i < valProf.size(); ++i) {
        std::string str = argToStr(ARGS(i));
        if (str == "Error") {
            str = "C" + std::to_string(i - C1 + 1);
        }
        std::cout << "profile " << str << ":\n";
        std::cout << valProf[i].in << "\n";
        for (uint64_t j = 0; j < profileSize; ++j) {
            std::cout << "(" << valProf[i].profile[j].first << ", " << valProf[i].profile[j].second << ") ";
        }
        std::cout << "\n";
    }
    //for chem
    for (uint64_t i = 0; i < chemProf.size(); ++i) {
        std::string str = argToStr(ARGS(i + TOTAL_COUNT));
        if (str == "Error") {
            str = "C" + std::to_string(i - C1 + TOTAL_COUNT);
        }
        std::cout << "profile " << str << ":\n";
        std::cout << chemProf[i].in << " " << chemProf[i].out << "\n";
        for (uint64_t j = 0; j < profileSize; ++j) {
            std::cout << "(" << chemProf[i].profile[j].first << ", " << chemProf[i].profile[j].second << ") ";
        }
        std::cout << "\n";
    }
    std::cout << "initial values:\n";
    for (uint64_t i = 0; i < f.size(); ++i) {
        std::string str = argToStr(ARGS(i));
        if (str == "Error") {
            str = "C" + std::to_string(i - C1 + 1);
        }
        std::cout << "initial " << str << ":\n";
        for (uint64_t j = 0; j < m; ++j) {
            std::cout << f[i](0, j) << " ";
        }
        std::cout << "\n";
    }

    std::vector<Matrix<double>> ans(n, Matrix<double>(m));
    std::cout << "solve method: kranc-nichalas\n";
    if (choice == 0) {
        auto startT = std::chrono::system_clock::now();
        ans = SolveIBVP(info, f, 0.7, R0, x_h, xi_h, 1.e-6, Method::KRANK_NICOLAS, ApproxLevel::NONE, ental, valProf, chemProf);
        auto endT = std::chrono::system_clock::now();
        std::cout << "=====Конец расчёта=====\n";
        auto time = std::chrono::duration_cast<duration_t>(endT - startT).count();
        std::cout << "Time: " << time << "\n";
        std::cout.flush();
        std::ofstream output(outputFile);
        uint64_t size = ans.size();
        //output.write(reinterpret_cast<const char *>(&size), sizeof(size));
        output << ans.size() << " " << ans[0].size().n << " " << ans[0].size().m << "\n";
        for (uint64_t i = 0; i < ans.size(); ++i) {
            //output << ans[i];
            for (uint64_t j = 0; j < ans[0].size().n; ++j) {
                for (uint64_t k = 0; k < ans[0].size().m; ++k) {
                    output << ans[i](j, k) << " ";
                }
                output << "\n";
            }
            output << "\n";
        }
        output.flush();
        output.close();
        //exit(0);
    }
    std::cout << "solved\n";
    setlocale(0, "");
    sf::ContextSettings settings;
    settings.antialiasingLevel = 6;
    std::cout << "Plot drawing\n";
    if (choice != 0) {
        std::ifstream input(outputFile);
        uint64_t ansSize, nSize, mSize;
        input >> ansSize >> nSize >> mSize;
        //input.read(reinterpret_cast<char *>(&ansSize), sizeof(ansSize));
        //input >> ansSize;
        ans.resize(ansSize, Matrix<double>(nSize, mSize));
        for (uint64_t i = 0; i < ans.size(); ++i) {
            //input >> ans[i];
            for (uint64_t j = 0; j < ans[0].size().n; ++j) {
                for (uint64_t k = 0; k < ans[0].size().m; ++k) {
                    input >> ans[i](j, k);
                }
            }
        }
        input.close();
        std::cout << "mat:\n" << ans[C1] << "\n";
    }
    std::vector<std::vector<sf::Vector2f>> points = {{}, {}, {}};
    ARGS toDraw = U;
    n = ans[toDraw].size().n;
    for (int i = 0; i < ans[toDraw].size().m; ++i) {
        points[0].push_back({ans[Y](0, i), ans[toDraw](0, i)});
        points[1].push_back({ans[Y](n / 2 - 1, i), ans[toDraw](n / 2 - 1, i)});
        points[2].push_back({ans[Y](n - 1, i), ans[toDraw](n - 1, i)});
    }
    std::vector<sf::Vector3f> plane;
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            //float x = minX + j * stepX;
            //float y = z / 10 * func1(x / 10) + (1 - z / 10) * func2(x / 10);
            //float y = std::cos(std::sqrt(std::pow(x, 2) + std::pow(z, 2)));
            //float y = std::pow(x, 2) / 40 - std::pow(z, 2) / 40;
            //plane.push_back({ans[Y](i, j) / R0, f[toDraw](i, j), j * R0 / 4});
            //plane.push_back({i, f[toDraw](i, j), j});
        }
    }

    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Plot for... idk? maybe U?", 7U, settings);
    window.setFramerateLimit(20);
    Plot2D plot(points), plot2(points);
    Line line1 = Line(Line::Type::LINE, sf::Color::Black, 2.f);
    Line line2 = Line(Line::Type::LINE, sf::Color::Green, 2.f);
    Line line3 = Line(Line::Type::LINE, sf::Color::Blue,  2.f);
    Selector selector;
    TextField field, animSpeedField;
    for (uint64_t i = 0; i < f.size(); ++i) {
        sf::String str;
        if (i > C1) {
            str = "C" + std::to_string(i - C1 + 1);
        } else {
            str = argToStr(ARGS(i));
        }
        selector.addVar(str);
    }
    plot.setColors({sf::Color::Black, sf::Color::Green, sf::Color::Blue});
    plot.setLines({line1, line2, line3});
    plot.setNames({"Plot for " + argToStr(toDraw) + ", 1", "Plot for " + argToStr(toDraw) + ", " + std::to_string(n / 2), "Plot for " + argToStr(toDraw) + ", " + std::to_string(n)});
    auto leg1 = plot.generateLegend(), leg2 = plot2.generateLegend();
    plot.setXName(L"Y/R0");
    plot.setYName(L"U/U0");
    plot.setSize(1200, 500);
    plot.setPosition(0, 0);
    plot2.setColors({sf::Color::Blue});
    plot2.setLines({line1});
    //plot2.addLegends({"Plot for " + argToStr(toDraw) + ", 1", "Plot for " + argToStr(toDraw) + ", " + std::to_string(n / 2), "Plot for " + argToStr(toDraw) + ", " + std::to_string(n)});
    plot2.setXName(L"x/R0");
    plot2.setYName(L"f");
    plot2.setPosition(0, 500);
    plot2.setSize(1200, 500);
    leg1->setPosition(plot.getPosition() + sf::Vector2f(plot.getSize().x, 0));
    leg1->update();
    leg2->setPosition(plot2.getPosition() + sf::Vector2f(plot2.getSize().x, 0));
    leg2->update();
    selector.setPosition(leg1->getPosition() + sf::Vector2f(0, leg1->getSize().y));
    field.setPosition(selector.getPosition() + sf::Vector2f(selector.getSize().x, 0));
    animSpeedField.setPosition(leg2->getPosition() + sf::Vector2f(0, leg2->getSize().y));
    uint64_t anim_idx = -1;
    float timePerFrame = 0.02f;
    auto changePlotSig = [&] () {
        toDraw = ARGS(selector.getCurrIdx());
        n = field.getInputLine().toInt64();
        points.clear();
        points.resize(3);
        for (int i = 0; i < ans[toDraw].size().m; ++i) {
            points[0].push_back({ans[Y](0, i), ans[toDraw](0, i)});
            points[1].push_back({ans[Y](n / 2 - 1, i), ans[toDraw](n / 2 - 1, i)});
            points[2].push_back({ans[Y](n - 1, i), ans[toDraw](n - 1, i)});
        }
        plot.setFuncs(points, {line1, line2, line3}, {"Plot for " + argToStr(toDraw) + ", 1", "Plot for " + argToStr(toDraw) + ", " + std::to_string(n / 2), "Plot for " + argToStr(toDraw) + ", " + std::to_string(n)});
        leg1 = plot.generateLegend();
        leg1->setPosition(plot.getPosition() + sf::Vector2f(plot.getSize().x, 0));
        leg1->update();
        anim_idx = -1;
    };
    selector.setSignal(Interactable::Signal::LOST_FOCUS, changePlotSig);
    field.setSignal(Interactable::Signal::LOST_FOCUS, changePlotSig);
    animSpeedField.setSignal(Interactable::Signal::LOST_FOCUS, [&] () {
        timePerFrame = field.getInputLine().toFloat64();
    });
    sf::Clock clock;

    while (window.isOpen()) {
        sf::Event event;

        while (window.pollEvent(event)) {
            sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
            sf::Vector2f mousePos = window.mapPixelToCoords(pixelPos);
            if (event.type == sf::Event::Closed || sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
                window.close();
            }
            if (event.type == sf::Event::KeyPressed && sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
                anim_idx = 0;
            }
            plot.handleEvent(event, mousePos);
            selector.handleEvent(event, mousePos);
            field.handleEvent(event, mousePos);
            //plot2.handleEvent(event, mousePos);
            animSpeedField.handleEvent(event, mousePos);
        }
        if (anim_idx != -1 && anim_idx < n && clock.getElapsedTime().asSeconds() > timePerFrame) {
            points.clear();
            points.resize(1);
            for (int i = 0; i < ans[toDraw].size().m; ++i) {
                points[0].push_back({ans[Y](anim_idx, i), ans[toDraw](anim_idx, i)});
            }
            plot2.setFuncs(points, {line1}, {"Plot for " + std::to_string(anim_idx)});
            ++anim_idx;
            leg2 = plot2.generateLegend();
            leg2->setPosition(plot2.getPosition() + sf::Vector2f(plot2.getSize().x, 0));
            leg2->update();
            clock.restart();
        }
        plot.update();
        selector.update();
        field.update();
        plot2.update();
        animSpeedField.update();

        window.clear(sf::Color(127, 127, 127));
        window.draw(plot);
        window.draw(plot2);
        window.draw(*leg1);
        window.draw(*leg2);
        window.draw(selector);
        window.draw(field);
        window.draw(animSpeedField);
        //window.draw(map);
        window.display();
    }
    return 0;
}