#include <iostream>
#include <ODUSolver/Chemical/ChemicalTask.hpp>

int main () {
    //Перумов "Что-то там горение чего-то"
    std::cout << "Its testing time\n";
    //system("ls");
    ChemicalSystem sys;
    sys.initFromFile("./test/ChemicTest/bufermm.txt");

    sys.setPressure(101'325);
    sys.setTemperature(2300);
    sys.addAdditive("M", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0});

    sys.addReaction("H2 + O2 <==> 2OH", 1.7 * std::pow(10, 7), 0, 24044);
    sys.addReaction("H + O2 <==> OH + O", 1.987 * std::pow(10, 8), 0, 8456);
    sys.addReaction("H2 + OH <==> H2O + H", 1.024 * std::pow(10, 2), 1.6, 1660);
    sys.addReaction("H2 + O <==> OH + H", 5.119 * std::pow(10, -2), 2.67, 3163);
    sys.addReaction("2OH <==> H2O + O", 1.506 * std::pow(10, 3), 1.14, 50);
    sys.addReaction("H + OH + M <==> H2O + M", 2.212 * std::pow(10, 10), -2.0, 0);
    sys.addReaction("2H + M <==> H2 + M", 9.791 * std::pow(10, 7), -0.6, 0);

    sys.printInfo(std::cout);

    sys.setConcentrations({0.5, 0.0, 0.0, 0.5, 0.0, 0.0});
    sys.rightPartGen();
    std::cout << "Testing: right part generation\n";
    auto rightPart = sys.getODE();
    auto gamma = sys.getY0();
    gamma.insert(gamma.begin(), 0);
    gamma.insert(gamma.end(), 0.09);
    gamma.insert(gamma.end(), 2300);
    for (auto el : gamma) {
        std::cout << el << "\n";
    }
    for (uint64_t i = 0; i < rightPart.size(); ++i) {
        std::cout << "\tFunc " << i + 1 << ": " << rightPart[i](gamma) << "\n";
    }
    std::cout << "Looks like we tested all other the place\n";
    return 0;
}