#include "FuncMaker.hpp"

double f (const std::vector<double> &X) {
    double x = X[0], y = X[1];
    //return -x - y + std::sin(-x) - std::acos(-1) * std::acos(-1);
    return -y + std::sin(-x) * std::sin(-x) + std::sin(x*x);
    //return 4 - 5;
}

int main () {
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
    return 0;
}