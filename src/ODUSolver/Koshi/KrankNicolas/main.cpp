#include "KrancNicolas.hpp"

int main () {
    std::cout.precision(PRECISION);
    std::cout.setf(std::ios_base::fixed);

    const double TIME_LIMIT = 10; //ограничение по времени от 0 до 10 секунд
    std::vector<std::string> system;
    std::string str;
    double xh, th;
    std::pair<std::vector<double>, std::vector<double>> res1, res2;
    Task task;
    FunctionalTree check;

    std::cout << "Введите начально-краевую задачу:\n";

    str = readLine();
    system.push_back(str); //считываем уравнение
    for (uint64_t i = 0; i < 3; ++i) { // считываем начальные условия
        str = readLine();
        system.push_back(str);
    }
    std::cout << "Введите размер шага для \"x\" и для \"t\":\n";
    std::cin >> xh >> th;
    std::cout << "Введите функцию для сравнения:\n";
    str = readLine();
    task = getTaskInfo(system); //считываем аналитическое решение

    uint64_t xSize = uint64_t((task.X[1] - task.X[2]) / xh), tSize = uint64_t(TIME_LIMIT / th);
    uint64_t S = xSize * tSize; //вычисляем площадь таблицы U
    std::cout << "Размер таблицы функции U: " << xSize << "x" << tSize << "\n";
    for (uint64_t i = 0; i < 3; ++i) {
        ApproxLevel approx = static_cast<ApproxLevel>(i);
        switch (approx) {
            case ApproxLevel::POINT2_ORDER1:
                std::cout << "Point 2 Order 1\n";
                break;
            case ApproxLevel::POINT2_ORDER2:
                std::cout << "Point 2 Order 2\n";
                break;
            case ApproxLevel::POINT3_ORDER2:
                std::cout << "Point 3 Order 2\n";
                break;
            default:
                std::cout << "Error\n";
                break;
        }
        auto ans1 = SolveIBVP(task, TIME_LIMIT, xh, th, Method::EXPLICIT, approx);
        auto ans2 = SolveIBVP(task, TIME_LIMIT, xh, th, Method::NOT_EXPLICIT, approx);
        auto ans3 = SolveIBVP(task, TIME_LIMIT, xh, th, Method::KRANK_NICOLAS, approx);
        double mistakeSum1 = 0, mistakeSum2 = 0, mistakeSum3 = 0;
        double X0 = task.X[0];
        check.reset(str, {"x", "t"});
        for (uint64_t i = 0; i < ans1.size(); ++i) {
            for (uint64_t j = 0; j < ans1[i].size(); ++j) {
                mistakeSum1 += std::abs(ans1[i][j] - check({X0 + j * xh, i * th}));
                mistakeSum2 += std::abs(ans2[i][j] - check({X0 + j * xh, i * th}));
                mistakeSum3 += std::abs(ans3[i][j] - check({X0 + j * xh, i * th}));
            }
        }
        std::cout << "Средняя погрешность явного метода: " << mistakeSum1 / S << "\n";
        std::cout << "Средняя погрешность неявного метода: " << mistakeSum2 / S << "\n";
        std::cout << "Средняя погрешность метода Кранка-Николаса: " << mistakeSum3 / S << "\n";
    }

    return 0;
}