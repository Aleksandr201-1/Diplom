#include <iostream>
#include <General/FuncMaker.hpp>

int main () {
    std::string line;
    std::cout << "> ";
    while (std::getline(std::cin, line)) {
        std::cout << "> " << FunctionalTree(line).calculate() << "\n\n";
        std::cout << "> ";
    }
    return 0;
}