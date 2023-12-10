#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "gnuplot-iostream.h"

void plot (const std::vector<std::string> &func, const std::vector<double> points) {
    static std::vector<std::string> colors = {"red", "green", "blue"};
    Gnuplot gp;

    gp << "set xlabel \"X\"\n";
    gp << "set ylabel \"Y\"\n";
    gp << "set xzeroaxis lw 1\n";
    gp << "set yzeroaxis lw 1\n";
    gp << "set grid\n";
    gp << "set title \"Plot\" font \"Helvetica Bold, 20\"\n";
    
    gp << "plot ";
    for (uint64_t i = 0; i < func.size(); ++i) {
        gp << func[i] << " title \"polynom " << i + 1 << "\"  lc rgb \"" << colors[i % colors.size()] << "\"";
        if (i == func.size() - 1) {
            if (points.size()) {
                gp << ", '-' " << " title \"our func\n";
                for (uint64_t j = 0; j < points.size(); j += 2) {
                    gp << points[j] << " " << points[j + 1] << "\n";
                }
                gp << "e\n";
            }
            gp << "\n";
        } else {
            gp << ",";
        }
    }
}

void plotCube () {
    Gnuplot gp;
    gp << "set contour\n";
    // gp << "splot sin(x)*cos(y)\n";
    gp << "splot sin(sqrt(x**2 + y**2))\n";
}

int main () {
    plotCube();
    return 0;
}