#include "RungeKutta.hpp"

std::pair<std::vector<double>, std::vector<double>> RungeKutta4 (const Task &task, double h) {
    double X0 = task.X0, Xn = task.Xn, a = task.a, b = task.b;
    const std::vector<FunctionalTree> &trees = task.trees;
    uint64_t order = task.n;

    std::vector<double> X, Y, Z;
    std::vector<std::vector<double>> Yi(task.Y.size());
    //std::cout << X0 << " " << Xn << "\n";
    for (double i = X0; i <= Xn; i += h) {
        X.push_back(i);
        //std::cout << i << "\n";
    }
    for (uint64_t i = 0; i < task.Y.size(); ++i) {
        Yi[i].push_back(task.Y[i]);
    }
    Y.push_back(a);
    Z.push_back(b);
    auto g = [&] (double x, double y, double z) -> double {
        static auto c2 = trees[0].getCoeff(1);
        return -trees[0]({x, 0, z, y}) / c2(x);
    };
    std::vector<std::vector<double>> K(order, std::vector<double>(4));
    //std::vector<double> delta(order);
    for (uint64_t i = 1; i < X.size(); ++i) {
        // double delta;
        // for (uint64_t j = 0; j < order; ++j) {
        //     K[j][0] = h * Yi[j][i - 1];
        //     delta = 1.0 / 6.0 * (K[j][0] + 2 * K[j][1] + 2 * K[j][2] + K[j][3]);
        //     Yi[j].push_back(Yi[j][i - 1] + delta);
        // }
        double delta_y, delta_z;
        double K1 = h * Z[i - 1];
        double L1 = h * g(X[i - 1], Y[i - 1], Z[i - 1]);
        double K2 = h * (Z[i - 1] + 1.0 / 2.0 * L1);
        double L2 = h * g(X[i - 1] + 1.0 / 2.0 * h, Y[i - 1] + 1.0 / 2.0 * K1, Z[i - 1] + 1.0 / 2.0 * L1);
        double K3 = h * (Z[i - 1] + 1.0 / 2.0 * L2);
        double L3 = h * g(X[i - 1] + 1.0 / 2.0 * h, Y[i - 1] + 1.0 / 2.0 * K2, Z[i - 1] + 1.0 / 2.0 * L2);
        double K4 = h * (Z[i - 1] + L3);
        double L4 = h * g(X[i - 1] + h, Y[i - 1] + K3, Z[i - 1] + L3);
        delta_y = 1.0 / 6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
        delta_z = 1.0 / 6.0 * (L1 + 2 * L2 + 2 * L3 + L4);
        Y.push_back(Y[i - 1] + delta_y);
        Z.push_back(Z[i - 1] + delta_z);
    }
    return std::make_pair(X, Y);
}

std::pair<std::vector<double>, std::vector<double>> RungeKutta6 (const Task &task, double h) {
    double X0 = task.X0, Xn = task.Xn, a = task.a, b = task.b;
    const std::vector<FunctionalTree> &trees = task.trees;

    std::vector<double> X, Y, Z;
    for (double i = X0; i <= Xn; i += h) {
        X.push_back(i);
    }
    Y.push_back(a);
    Z.push_back(b);
    auto g = [&] (double x, double y, double z) -> double {
        static auto c2 = trees[0].getCoeff(1);
        return -trees[0]({x, 0, z, y}) / c2(x);
    };
    for (uint64_t i = 1; i < X.size(); ++i) {
        double delta_y, delta_z;
        double K1 = h * Z[i - 1];
        double L1 = h * g(X[i - 1], Y[i - 1], Z[i - 1]);
        double K2 = h * (Z[i - 1] + 1.0 / 2.0 * L1);
        double L2 = h * g(X[i - 1] + 1.0 / 2.0 * h, Y[i - 1] + 1.0 / 2.0 * K1, Z[i - 1] + 1.0 / 2.0 * L1);
        double K3 = h * (Z[i - 1] + 1.0 / 2.0 * L2);
        double L3 = h * g(X[i - 1] + 1.0 / 2.0 * h, Y[i - 1] + 1.0 / 2.0 * K2, Z[i - 1] + 1.0 / 2.0 * L2);
        double K4 = h * (Z[i - 1] + L3);
        double L4 = h * g(X[i - 1] + h, Y[i - 1] + K3, Z[i - 1] + L3);
        delta_y = 1.0 / 6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
        delta_z = 1.0 / 6.0 * (L1 + 2 * L2 + 2 * L3 + L4);
        Y.push_back(Y[i - 1] + delta_y);
        Z.push_back(Z[i - 1] + delta_z);
    }
    return std::make_pair(X, Y);
}