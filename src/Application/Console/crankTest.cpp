#include <ODUSolver/CrankNicolson/CrankNicolson.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <Input/InputHandle.hpp>
#include <GUI/Plot2D.hpp>
#include <GUI/Plot3D.hpp>
#include <GUI/HeatMap.hpp>

//Численное моделирование течения струй газа с неравновесным химическим процессом

//Newton для основного текста
//Must-Type для математических формул
//Курсив запрещён

std::vector<std::pair<double, double>> Lagrange (const std::vector<std::pair<double, double>> &points) {
    std::vector<std::pair<double, double>> polynom(points);
    for (uint64_t i = 0; i < points.size(); ++i) {
        double w_el = 1;
        for (uint64_t j = 0; j < i; ++j) {
            w_el *= (points[i].first - points[j].first);
        }
        for (uint64_t j = i + 1; j < points.size(); ++j) {
            w_el *= (points[i].first - points[j].first);
        }
        polynom[i].second = w_el;
    }
    for (uint64_t i = 0; i < points.size(); ++i) {
        polynom[i].second = points[i].second / polynom[i].second;
    }
    return polynom;
}

double LagrangeFunc (const std::vector<std::pair<double, double>> &polynom, double x) {
    double ans = 0;
    for (uint64_t i = 0; i < polynom.size(); ++i) {
        double tmp = 1;
        for (uint64_t j = 0; j < i; ++j) {
            tmp *= (x - polynom[j].first);
        }
        for (uint64_t j = i + 1; j < polynom.size(); ++j) {
            tmp *= (x - polynom[j].first);
        }
        tmp *= polynom[i].second;
        ans += tmp;
    }
    return ans;
}

void initializer (std::vector<Matrix<double>> &f, std::vector<std::pair<double, double>> &points, double xn, int i) {
    auto polynom = Lagrange(points);
    double step = xn / f[i].size().n;
    for (uint64_t idx = 0; idx < f[i].size().n; ++idx) {
        f[i](0, idx) = LagrangeFunc(polynom, idx * step);
    }
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
    info.fileInput = "./test5.txt";
    info.order = 4;
    info.way = 1;
    info.algo = IterationAlgo::ZEIDEL;
    info.method = SolveMethod::RUNGE_KUTTA;
    info.butcher = createButcherTable(info.method, info.order, info.way);
    info.approx = 0.01;
    argsHandler(args, info);

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
    //???
    //page №25 for V
    for (uint64_t i = 0; i < initProfileCount; ++i) {
        std::string valName;
        in >> valName;
        double inV;
        in >> inV;
        valProf[strToArgs(valName)].in = inV;
        //profiles[strToArgs(valName)].out = outV;
    }
    in >> initProfileCount;
    std::vector<std::string> valueNames(initProfileCount);
    for (uint64_t i = 0; i < initProfileCount; ++i) {
        in >> valueNames[i];
    }
    for (uint64_t i = 0; i < valProf.size(); ++i) {
        for (uint64_t j = 0; j < profileSize; ++j) {
            valProf[i].profile[j].first = j * xi_h;
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
    for (uint64_t i = 0; i < chemProf.size(); ++i) {
        for (uint64_t j = 0; j < profileSize; ++j) {
            chemProf[i].profile[j].first = j * xi_h;
        }
    }
    for (uint64_t i = 0; i < profileSize; ++i) {
        for (uint64_t j = 0; j < chemProf.size(); ++j) {
            double val;
            in >> val;
            //chemProf[strToArgs(valueNames[j])].profile[i].second = val;
            chemProf[j].profile[i].second = val;
        }
    }
    in >> usingTurb;
    std::cout << "Turbulence mode: " << usingTurb << "\n";
    std::cout << "Введите количество добавок и добавки:\n";
    in >> additiveCount;
    //std::cout << "filename: " << filename << "\n";
    // sys.setTemperature(Temp);
    // sys.setPressure(Press);
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
    for (int i = 0; i < insideGamma.size(); ++i) {
        insideGamma[i] /= sum1;
        outsideGamma[i] /= sum2;
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
    //sys.printInfo(std::cout);
    auto ental = sys.getEnthalpy();
    info.task = &sys;

    std::vector<Matrix<double>> f(TOTAL_COUNT + sys.getSubstanceMasses().size() - 1, Matrix<double>(n, m));
    valProf[V].in = valProf[W].in = valProf[U].in;
    //переход к физическим величинам
    for (uint64_t i = 0; i < C1; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](0, j) = valProf[i].getPhysical((double)j * xi_h * profileSize / m);
        }
    }
    for (uint64_t i = C1; i < f.size(); ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](0, j) = chemProf[i - C1].getPhysical((double)j * xi_h * profileSize / m);
        }
    }
    double Xi_0 = 0;

    //RHO
    auto pointGamma = insideGamma;
    for (uint64_t i = 0; i < m; ++i) {
        std::cout << "press for rho " << i << ":\n";
        for (uint64_t j = 0; j < pointGamma.size(); ++j) {
            pointGamma[j] = f[C1 + j](0, i);
            std::cout << f[C1 + j](0, i) << " "; 
        }

        sys.setTemperature(f[T](0, i));
        sys.setPressure(f[P](0, i));
        sys.setConcentrations(pointGamma, ConcentrationMode::MOLAR_MASS);
        std::cout << "T: " << f[T](0, i) << " P: " << f[P](0, i) << " rho: " << sys.getRho() << "\n";
        f[RHO](0, i) = sys.getRho();
    }
    valProf[RHO].in = f[RHO](0, 0);
    valProf[P].in = valProf[U].in * valProf[U].in * valProf[RHO].in;
    //Y
    double nu = 1.0;
    f[Y](0, 0) = 0;
    for (uint64_t j = 0; j < m; ++j) {
        double xi = j * xi_h + Xi_0;
        f[Y](0, j) = std::pow(std::pow(xi, nu + 1) * (nu + 1) / (valProf[RHO].toNormal(f[RHO](0, j)) * valProf[U].toNormal(f[U](0, j))), 1.0 / (nu + 1));
        //f[Y](0, j) = std::pow(std::pow(xi, nu + 1) * (nu + 1) / (f[RHO](0, j) * f[U](0, j)), 1.0 / (nu + 1));
    }
    // f[Y](0, 0) = 0;
    // for (uint64_t j = 1; j < n; ++j) {
    //     auto yf = [&] (uint64_t idx) -> double {
    //         return std::pow(xi_h * idx + Xi_0, nu) / valProf[RHO].toNormal(f[RHO](0, idx - 1)) / valProf[U].toNormal(f[U](0, idx - 1));
    //     };
    //     auto xf = [&] (uint64_t idx) -> double {
    //         return idx * xi_h + Xi_0;
    //     };
    //     f[Y](0, j) = std::pow((nu + 1) * IntegralTrapeze(yf, xf, j), 1.0 / (nu + 1));
    // }
    // f[Y](0, 0) = 0;
    // //f[Y](0, 1) = std::pow(std::pow(xi_h, nu) / (valProf[RHO].toNormal(f[RHO](0, 1)) * valProf[U].toNormal(f[U](0, 1))) * xi_h, 1.0 / (nu + 1));
    // f[Y](0, 1) = std::pow(std::pow(xi_h, nu) / (f[RHO](0, 1) * f[U](0, 1)) * xi_h, 1.0 / (nu + 1));
    // for (uint64_t j = 2; j < m; ++j) {
    //     double xi = xi_h * (j - 1) + Xi_0;
    //     //f[Y](0, j) = std::pow(xi / f[Y](0, j - 1), nu) / (valProf[RHO].toNormal(f[RHO](0, j - 1)) * valProf[U].toNormal(f[U](0, j - 1))) * 2*xi_h + f[Y](0, j - 2);
    //     f[Y](0, j) = std::pow(xi / f[Y](0, j - 1), nu) / (f[RHO](0, j - 1) * f[U](0, j - 1)) * 2*xi_h + f[Y](0, j - 2);
    // }
//Шлихтинг толстая книжка про вязкость
    //MU
    for (uint64_t j = 0; j < m; ++j) {
        f[MU](0, j) = turbulence[usingTurb](f, valProf, chemProf, R0, 0, j);
    }
    //valProf[MU].in = 1;//valProf[U].in * valProf[RHO].in * R0;
    //valProf[MU].in = 1.0 / 0.18e-4;

    //J
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[J](i, j) = 0;
            for (uint64_t k = 0; k < ental.size(); ++k) {
                f[J](i, j) += ental[k](f[T](i, j)) * f[C1 + k](i, j);
            }
            //f[J](i, j) /= (f[U](i, j) * f[U](i, j));
            f[J](i, j) += 1.0 / 2.0 * (f[U](i, j) * f[U](i, j) + f[V](i, j) * f[V](i, j) + f[W](i, j) * f[W](i, j));
        }
    }
    valProf[J].in = valProf[U].in * valProf[U].in;
    // for (uint64_t i = 0; i < valProf[J].profile.size(); ++i) {
    //     valProf[J].profile[i].second = f[J](0, i) / valProf[J].in;
    //     valProf[RHO].profile[i].second = f[RHO](0, i) / valProf[RHO].in;
    //     valProf[MU].profile[i].second = f[MU](0, i);
    // }

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

    //нормировка
    for (uint64_t i = 0; i < C1; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](0, j) = valProf[i].toNormal(f[i](0, j));
        }
    }
    for (uint64_t i = C1; i < f.size(); ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            f[i](0, j) = chemProf[i - C1].toNormal(f[i](0, j));
        }
    }
    for (uint64_t j = 0; j < m; ++j) {
        f[MU](0, j) *= valProf[MU].in;
    }

    std::cout << "initial profiles:\n";
    //for vals
    for (uint64_t i = 0; i < valProf.size(); ++i) {
        std::string str = argToStr(ARGS(i));
        if (str == "Error") {
            str = "C" + std::to_string(i - C1 + 1);
        }
        std::cout << "profile " << str << ":\n";
        std::cout << valProf[i].in << "\n";// << profiles[i].out << "\n";
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
            //f[P](0, i) = sys.getPressure();
            //f[MU](0, i) = 0.4;
        }
        std::cout << "\n";
    }
//     //Y
//     f[Y](0, 0) = 0.1 * R0;
//     for (uint64_t i = 1; i < m; ++i) {
//         double xi = i * (R0 / 100) + 0;
//         //f[Y](0, i) = f[Y](0, 0) + 0.1 * R0 / 20 * i;
//         f[Y](0, i) = 1.0 / std::pow(f[Y](0, i - 1), 1.0) * xi / (f[RHO](0, i - 1) * f[U](0, i - 1));
//         f[Y](0, i) = f[Y](0, i) * R0 / m + f[Y](0, i - 1);
//     }
    // auto yf = [&] (uint64_t k) -> double {
    //     double xi = k * (R0 / 20) + 0;
    //     return 1.0 / std::pow(f[Y](0, k), 1.0) * std::pow(xi, 1.0) / (f[RHO](0, k) * f[U](0, k));
    // };
    // for (uint64_t j = 1; j < f[Y].size().m; ++j) {
    //     f[Y](0, j) = f[Y](0, j - 1) + (R0 / 20) * yf(j - 1);
    // }
    std::vector<Matrix<double>> ans(n, Matrix<double>(m));
    std::cout << "solve method: kranc-nichalas\n";
    ans = SolveIBVP(info, f, 0.4, R0, x_h, xi_h, Method::KRANK_NICOLAS, ApproxLevel::NONE, ental, valProf, chemProf);
    std::cout << "solved\n";
    sys.setConcentrations(insideGamma, ConcentrationMode::MOLAR_MASS);
    info.solution = ChemicalSolver(info.method, *static_cast<ChemicalSystem*>(info.task), info.butcher, info.h_min, info.h_max, info.h_last, info.algo, info.approx, ReactionType::ADIABAT_CONST_RHO);
    setlocale(0, "");
    sf::ContextSettings settings;
    settings.antialiasingLevel = 6;
    uint64_t chemCount = f.size() - C1;
    std::vector<std::vector<sf::Vector2f>> float_ans(chemCount);
    // for (uint64_t i = 0; i < chemCount; ++i) {
    //     for (uint64_t j = 0; j < f[C1 + i].size().n; ++j) {
    //         float_ans[i].push_back({info.h * j, f[C1 + i](j, m / 2)});
    //     }
    // }
    for (uint64_t i = 0; i < chemCount; ++i) {
        for (uint64_t j = 0; j < info.solution[i].size(); ++j) {
            if (j % 10 == 0) {
                float_ans[i].push_back({info.solution[0][j], info.solution[i + 1][j]});
            }
        }
    }
    exit(0);
    std::cout << "Plot drawing\n";
    int a;
    std::cin >> a;

    std::vector<std::vector<sf::Vector2f>> points = {{}, {}, {}}, points2 = {{}};
    HeatMap map;
    Matrix<float> matrix(ans[U].size().m * 2, ans[U].size().n, 0);
    for (int i = 0; i < ans[U].size().n; ++i) {
        for (int j = 0; j < ans[U].size().m; ++j) {
            matrix(j + ans[U].size().m / 2, i)     = ans[U](i, j);
            matrix(ans[U].size().m / 2 - 1 - j, i) = ans[U](i, j);
        }
    }
    map.setMap(matrix);
    ARGS toDraw = U;
    for (int i = 0; i < ans[toDraw].size().m; ++i) {
        points[0].push_back({ans[Y](0, i) / R0, ans[toDraw](0, i)});
        points[1].push_back({ans[Y](n / 2 - 1, i) / R0, ans[toDraw](n / 2 - 1, i)});
        points[2].push_back({ans[Y](n - 1, i) / R0, ans[toDraw](n - 1, i)});
        //points2[0].push_back({i * R0 / 4, ans[toDraw](i, m / 2 - 1)});
    }
    int chemSysSize = info.solution.size();
    for (int i = 0; i < info.solution[chemSysSize - 1].size(); ++i) {
        if (i % 10 == 0) {
        points2[0].push_back({info.solution[0][i], info.solution[chemSysSize - 1][i]});
        }
    }
    std::vector<sf::Vector3f> plane;
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            //float x = minX + j * stepX;
            //float y = z / 10 * func1(x / 10) + (1 - z / 10) * func2(x / 10);
            //float y = std::cos(std::sqrt(std::pow(x, 2) + std::pow(z, 2)));
            //float y = std::pow(x, 2) / 40 - std::pow(z, 2) / 40;
            //plane.push_back({ans[Y](i, j) / R0, f[toDraw](i, j), j * R0 / 4});
            plane.push_back({i, f[toDraw](i, j), j});
        }
    }
    exit (0);
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Plot for... idk? maybe U?", 7U, settings);
    window.setFramerateLimit(20);
    Plot2D plot(points), plot2(points2);
    Plot3D suffering(plane, n, m);
    suffering.setPosition(900, 550);
    suffering.setLineColor(sf::Color::Red);
    suffering.setPlaneColor(sf::Color::Green);
    plot.setColors({sf::Color::Black, sf::Color::Green, sf::Color::Blue});
    //plot.addLegends({"Plot for " + argToStr(toDraw) + ", 1", "Plot for " + argToStr(toDraw) + ", " + std::to_string(n / 2), "Plot for " + argToStr(toDraw) + ", " + std::to_string(n)});
    plot.setXName(L"y/R0");
    plot.setYName(L"f");
    plot.setPosition(0, 0);
    plot2.setColors({sf::Color::Black, sf::Color::Green, sf::Color::Blue});
    //plot2.addLegends({"Plot for " + argToStr(toDraw) + ", 1", "Plot for " + argToStr(toDraw) + ", " + std::to_string(n / 2), "Plot for " + argToStr(toDraw) + ", " + std::to_string(n)});
    plot2.setXName(L"x/R0");
    plot2.setYName(L"f");
    plot2.setPosition(0, 550);
    //plot.scale(2, 2);

    Plot2D chemPlot(float_ans);
    chemPlot.setColors({sf::Color::Red, sf::Color::Green, sf::Color::Blue, sf::Color::Black, sf::Color::Yellow, sf::Color::Magenta, sf::Color::Cyan, sf::Color(124, 145, 243), sf::Color(145, 123, 243), sf::Color(243, 145, 132), sf::Color(154, 234, 128)});
    //chemPlot.addLegends({"CO", "CO2", "H2", "OH", "H2O", "O2", "N2", "NO", "H", "O", "N"});
    chemPlot.setXName("T(us)");
    chemPlot.setYName("rho");
    chemPlot.setPosition(900, 0);
    while (window.isOpen()) {
        sf::Event event;

        while (window.pollEvent(event)) {
            sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
            sf::Vector2f mousePos = window.mapPixelToCoords(pixelPos);
            if (event.type == sf::Event::Closed || sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
                window.close();
            }
            plot.handleEvent(event, mousePos);
            plot2.handleEvent(event, mousePos);
            chemPlot.handleEvent(event, mousePos);
            suffering.handleEvent(event, mousePos);
        }
        plot.update();
        plot2.update();
        chemPlot.update();
        suffering.update();

        window.clear(sf::Color(127, 127, 127));
        window.draw(plot);
        window.draw(plot2);
        window.draw(chemPlot);
        window.draw(suffering);
        window.draw(map);
        window.display();
    }
    return 0;
}