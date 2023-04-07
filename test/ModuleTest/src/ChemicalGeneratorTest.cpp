#include <iostream>
#include <ChemicalGenerator/RightPartGen.hpp>

int main () {
    std::cout << "Its testing time\n";
    //system("ls");
    ChemicalSystem sys;
    sys.initFromFile("./test/ChemicTest/bufermm.txt");
    sys.setPressure(100'000);
    sys.setTemperature(1000);
    //sys.setDensity(0.95);
    std::cout << "Testing: Chemical Generator\n";
    sys.printInfo(std::cout);
    std::cout << "Testing: Adding reactions\n";
    sys.addReaction("H2 + O2 => 2OH");
    sys.setReactionParameters(0, 1.7 * std::pow(10, 7), 0, 24044);
    sys.addReaction("H + O2 => OH + O");
    sys.setReactionParameters(1, 1.987 * std::pow(10, 8), 0, 8456);
    sys.addReaction("H2 + OH => H2O + H");
    sys.setReactionParameters(2, 1.024 * std::pow(10, 2), 1.6, 1660);
    sys.addReaction("H2 + O => OH + H");
    sys.setReactionParameters(3, 5.119 * std::pow(10, -2), 2.67, 3163);
    sys.addReaction("2OH => H2O + O");
    sys.setReactionParameters(4, 1.506 * std::pow(10, 3), 1.14, 50);
    sys.addReaction("H + OH => H2O");
    sys.setReactionParameters(5, 2.212 * std::pow(10, 16), -2.0, 0);
    //sys.addReaction("2H => H2");
    //sys.setReactionParameters(6, 9.791 * std::pow(10, 10), -0.6, 0);
    sys.printInfo(std::cout);
    std::cout << "Testing: right part generation\n";
    auto rightPart = sys.rightPartGen();
    std::vector<double> X = {0, 0.5, 0, 0, 0.5, 0};
    for (uint64_t i = 0; i < rightPart.size(); ++i) {
        std::cout << "\tFunc " << i + 1 << ": " << rightPart[i](X) << "\n";
    }
    std::cout << "Looks like we tested all other the place\n";
    return 0;
}