#include "FuncMaker.hpp"
#include <chrono>

using duration_t = std::chrono::microseconds;

double f (const std::vector<double> &X) {
    double x = X[0], y = X[1];
    //return -x - y + std::sin(-x) - std::acos(-1) * std::acos(-1);
    return -y + std::sin(-x) * std::sin(-x) + std::sin(x*x);
    //return 4 - 5;
}

double test_func (double x) {
    return 11 - x*x*x*(3*x - 8);
}

double useAsFunctor (const std::function<double(const std::vector<double> &)> &func, const std::vector<double> &X) {
    return func(X);
}

int main () {
    std::chrono::time_point <std::chrono::system_clock> startt, endt;
    uint64_t time = 0;
    std::vector<std::string> var_list = {"x", "y"};
    FunctionalTree tree("-y + sin(-x)**2 + sin(y**2)", var_list);
    FunctionalTree tr;
    std::ifstream file("out.txt");
    file >> tr;
    file.close();
    //FunctionalTree tree("4 - 5", {"x", "y"});
    double miss = 0;

    for (uint64_t i = 0; i < 1000; ++i) {
        miss += std::abs(tree({1.0*i, 2}) - tr({1.0*i, 2}));
    }

    if (miss == 0.0) {
        std::cout << "Good\n" << tree << "\n" << tr << "\n";
    } else {
        std::cout << "Bad\n";
        tree.printFunc();
        std::cout << "\n----------\n";
        tree.printTree();
        std::cout << "\n";
    }
    tr.reset("(11 - x**3 * (3*x - 8))", {"x"});
    std::cout << tr(-1) << "\n" << tr << "\n";

    startt = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < 100000; ++i) {
        double val = tr(0);
    }
    endt = std::chrono::system_clock::now();
    time += std::chrono::duration_cast <duration_t>(endt - startt).count();
    std::cout << "FuncMaker time: " << time << "mcs\n";

    auto lambda_func = [] (double x) -> double {
        return 11 - x*x*x*(3*x - 8);
    };
    startt = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < 100000; ++i) {
        double val = lambda_func(0);
    }
    endt = std::chrono::system_clock::now();
    time += std::chrono::duration_cast <duration_t>(endt - startt).count();
    std::cout << "Lambda func time: " << time << "mcs\n";

    startt = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < 100000; ++i) {
        double val = test_func(0);
    }
    endt = std::chrono::system_clock::now();
    time += std::chrono::duration_cast <duration_t>(endt - startt).count();
    std::cout << "Normal func time: " << time << "mcs\n";

    std::cout << "Use as Functor: " << useAsFunctor(tr, {0}) << "\n";

    return 0;
}