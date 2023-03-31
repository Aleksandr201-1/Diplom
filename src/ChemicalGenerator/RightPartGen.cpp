#include "RightPartGen.hpp"
#include <iostream>

// std::string AtomToString (Atom atom) {
//     switch (atom) {
//         case Atom::H:
//             return "H";
//         case Atom::O:
//             return "O";
//         default:
//             return "";
//     }
// }

// std::string SubstanceToString (Substance substance) {
//     switch (substance) {
//         case Substance::H2:
//             return "H2";
//         case Substance::O2:
//             return "O2";
//         case Substance::H2O:
//             return "H2O";
//         case Substance::OH:
//             return "OH";
//         case Substance::H:
//             return "H";
//         case Substance::O:
//             return "O";
//         default:
//             return "";
//     }
// }

//std::vector<std::string> ChemicalReaction::atoms = {"H", "O"};
//std::vector<std::string> ChemicalReaction::substances = {"H2", "O2", "H2O", "OH", "H", "O"};

std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> tt;

std::string ChemicalReaction::readAtom (const std::string &sub, uint64_t &i) const {
    std::string ans;
    ans += sub[i];
    ++i;
    if (i < sub.size() && sub[i] >= 'a' && sub[i] <= 'z') {
        ans += sub[i];
        ++i;
    }
    return ans;
}

std::string ChemicalReaction::readSubstance (const std::string &sub, uint64_t &i) const {
    std::string ans;
    while (i < sub.size() && sub[i] != '+') {
        ans += sub[i];
        ++i;
    }
    return ans;
}

uint64_t ChemicalReaction::readCoeff (const std::string &sub, uint64_t &i) const {
    uint64_t coeff = 1;
    std::string num;
    while (i < sub.size() && sub[i] >= '0' && sub[i] <= '9') {
        num += sub[i];
        ++i;
    }
    if (!num.empty()) {
        coeff = std::stoull(num);
    }
    return coeff;
}

std::vector<uint64_t> ChemicalReaction::getSubstanceContent (const std::string &sub) const {
    std::vector<uint64_t> ans(atoms.size(), 0);
    uint64_t i = 0;
    while (i < sub.size()) {
    //for (uint64_t i = 0; i < sub.size(); ++i) {
        std::string atom = readAtom(sub, i);
        uint64_t coeff = readCoeff(sub, i);
        auto it = std::find(atoms.cbegin(), atoms.cend(), atom);
        uint64_t idx = std::distance(atoms.cbegin(), it);
        if (idx == ans.size()) {
            throw std::logic_error("function \"getSubstanceContent\": atom \"" + atom + "\" not presented in list");
        }
        ans[idx] += coeff;
    }
    return ans;
}

void ChemicalReaction::initTable () {
    //for (uint64_t i = 0; i < atoms.size(); ++i) {
        for (uint64_t i = 0; i < substances.size(); ++i) {
            auto coeffs = getSubstanceContent(substances[i]);
            for (uint64_t j = 0; j < atoms.size(); ++j) {
                tt[atoms[j]][substances[i]] = coeffs[j];
            }
        }
    //}
}

bool ChemicalReaction::checkForCorrect () const {
    std::vector<uint64_t> in(atoms.size(), 0), out(atoms.size(), 0);
    for (uint64_t i = 0; i < atoms.size(); ++i) {
        for (uint64_t j = 0; j < input.size(); ++j) {
            in[i] += input[j].second * tt[atoms[i]][substances[input[j].first]];
        }
        for (uint64_t j = 0; j < output.size(); ++j) {
            out[i] += output[j].second * tt[atoms[i]][substances[output[j].first]];
        }
    }
    return in == out;
}

// std::unordered_map<std::string, Atom> atomTable = {
//     {"H", Atom::H},
//     {"O", Atom::O}
// };
// std::unordered_map<std::string, Substance> substanceTable = {
//     {"H2",  Substance::H2},
//     {"O2",  Substance::O2},
//     {"H2O", Substance::H2O},
//     {"OH",  Substance::OH},
//     {"H",   Substance::H},
//     {"O",   Substance::O}
// };
std::unordered_map<std::string, double> atomMas = {
    {"H", 1.00784},
    {"O", 15.99903}
};

//std::vector<std::pair<Atom, double>> atomMas;

//Matrix<uint64_t> tableOfSubstance(to_underlying(Atom::COUNT_OF_ATOMS), to_underlying(Substance::COUNT_OF_SUBSTANCE));

ChemicalReaction::ChemicalReaction (const std::string &str, const std::vector<std::string> &atoms, const std::vector<std::string> &substances) : atoms(atoms), substances(substances) {
    //initTable();
    std::cout << "  ";
    //std::pair<std::string, std::unordered_map<std::string, uint64_t>> &elt = tt;
    for (auto &elt : tt) {
        for (auto &ell : elt.second) {
            std::cout << ell.first << "    ";
        }
        break;
    }
    
    std::cout << "\n";
    for (auto &el : tt) {
        std::cout << el.first << " ";
        for (auto &ell : el.second) {
            std::cout << ell.second << "    ";
        }
        std::cout << "\n";
    }
    setReaction(str);
}

ChemicalReaction::ChemicalReaction (const ChemicalReaction &reaction) : atoms(reaction.atoms), substances(reaction.substances) {
    input = reaction.input;
    output = reaction.output;
    A = reaction.A;
    n = reaction.n;
    E = reaction.E;
}

ChemicalReaction::ChemicalReaction (ChemicalReaction &&reaction) : atoms(reaction.atoms), substances(reaction.substances) {
    input = std::move(reaction.input);
    output = std::move(reaction.output);
    A = reaction.A;
    n = reaction.n;
    E = reaction.E;
}

ChemicalReaction::~ChemicalReaction () {}

//H + OH => H2O
void ChemicalReaction::setReaction (const std::string &str) {
    std::string tmp = str;
    //std::remove_if(tmp.begin(), tmp.end(), isspace);
    tmp.erase(std::remove_if(tmp.begin(), tmp.end(), isspace), tmp.end());
    std::cout << tmp << "\n";
    uint64_t div = tmp.find('=');
    std::string inputStr = tmp.substr(0, div), outputStr = tmp.substr(div + 2, tmp.size() - div);
    std::cout << inputStr << "||" << outputStr << "\n";
    uint64_t i = 0;
    //std::string &str = {inputStr, outputStr};
    while (i < inputStr.size()) {
        uint64_t coeff = readCoeff(inputStr, i);
        std::string sub = readSubstance(inputStr, i);
        auto it = std::find(substances.cbegin(), substances.cend(), sub);
        uint64_t idx = std::distance(substances.cbegin(), it);
        if (idx == substances.size()) {
            throw std::logic_error("operation \"setReaction\": substance \"" + sub + "\" not presented in list");
        }
        input.push_back({idx, coeff});
        ++i;
    }
    i = 0;
    while (i < outputStr.size()) {
        uint64_t coeff = readCoeff(outputStr, i);
        std::string sub = readSubstance(outputStr, i);
        auto it = std::find(substances.cbegin(), substances.cend(), sub);
        uint64_t idx = std::distance(substances.cbegin(), it);
        if (idx == substances.size()) {
            throw std::logic_error("operation \"setReaction\": substance \"" + sub + "\" not presented in list");
        }
        output.push_back({idx, coeff});
        ++i;
    }
    std::cout << "==COEFFS INPUT==\n";
    for (auto el : input) {
        std::cout << el.first << "==" << el.second << "\n";
    }
    std::cout << "==COEFFS OUTPUT==\n";
    for (auto el : output) {
        std::cout << el.first << "==" << el.second << "\n";
    }
    std::cout << "==CHECK FOR CORRECTNESS\n";
    if (checkForCorrect()) {
        std::cout << "fine\n";
    } else {
        std::cout << "not fine\n";
    }
}

void ChemicalReaction::setInput (const std::vector<std::pair<uint64_t, uint64_t>> &in) {
    input = in;
}

void ChemicalReaction::setOutput (const std::vector<std::pair<uint64_t, uint64_t>> &out) {
    output = out;
}

void ChemicalReaction::setParameters (double A, double n, double E) {
    this->A = A;
    this->n = n;
    this->E = E;
}

const std::vector<std::pair<uint64_t, uint64_t>> &ChemicalReaction::getInput () const {
    return input;
}

const std::vector<std::pair<uint64_t, uint64_t>> &ChemicalReaction::getOutput () const {
    return output;
}

std::tuple<double, double, double> ChemicalReaction::getParameters () const {
    return std::make_tuple(A, n, E);
}

ChemicalReaction &ChemicalReaction::operator= (const ChemicalReaction &reaction) {
    input = reaction.input;
    output = reaction.output;
    A = reaction.A;
    n = reaction.n;
    E = reaction.E;
    return *this;
}

ChemicalReaction &ChemicalReaction::operator= (ChemicalReaction &&reaction) {
    input = std::move(reaction.input);
    output = std::move(reaction.output);
    A = reaction.A;
    n = reaction.n;
    E = reaction.E;
    return *this;
}


ChemicalSystem::ChemicalSystem () {}

ChemicalSystem::~ChemicalSystem () {}

void ChemicalSystem::setAtomList (const std::vector<std::string> &list) {
    atoms = list;
}

void ChemicalSystem::setSubstanceList (const std::vector<std::string> &list) {
    substances = list;
}

void ChemicalSystem::addReaction (const std::string &reaction) {
    reactions.push_back(ChemicalReaction(reaction, atoms, substances));
}

void ChemicalSystem::addReaction (const ChemicalReaction &reaction) {
    reactions.push_back(reaction);
}

void ChemicalSystem::setPressure (double rho) {
    this->rho = rho;
}

double ChemicalSystem::getPressure () const {
    return rho;
}

void ChemicalSystem::setTemperature (double T) {
    this->T = T;
}

double ChemicalSystem::getTemperature () const {
    return T;
}

void ChemicalSystem::setReactionParameters (uint64_t i, double A, double n, double E) {
    reactions[i].setParameters(A, n, E);
}

std::tuple<double, double, double> ChemicalSystem::getReactionParameters (uint64_t i) const {
    return reactions[i].getParameters();
}

uint64_t ChemicalSystem::getCount () const {
    return reactions.size();
}

ChemicalReaction ChemicalSystem::operator[] (uint64_t i) const {
    return reactions[i];
}

ChemicalReaction &ChemicalSystem::operator[] (uint64_t i) {
    return reactions[i];
}

double K (double A, double n, double T, double E) {
    return A * std::pow(T, n) * std::exp(-E / T);
}

std::tuple<std::vector<double>, std::vector<double>> readInfo (const std::string &str) {
    std::vector<double> phi, entalpil;
    return std::make_tuple(phi, entalpil);
}

std::vector<std::function<double(const std::vector<double> &)>> rightPartGen (const ChemicalSystem &system) {
    std::vector<std::function<double(const std::vector<double> &)>> rightPart;
    double A, n, E;
    double rho, T;
    rho = system.getPressure();
    T = system.getTemperature();
    double R = 8.31;
    std::vector<double> phi, entalpil;
    std::tie(phi, entalpil) = readInfo("");
    double phi1, phi2, phi3, phi4, phi5, phi6, phi7;
    auto Phi0 = [=] (double t) -> double {
        return phi1 + phi2 * std::log(t) + phi3 * std::pow(t, -2) + phi4 * std::pow(t, -1) + phi5 * t + phi6 * std::pow(t, 2) + phi7 * std::pow(t, 3);
    };
    double entalpil1, entalpil2;
    auto G0 = [=] (double t) -> double {
        return entalpil1 - entalpil2 - t * Phi0(t);
    };
    for (uint64_t i = 0; i < system.getCount(); ++i) {
        std::tie(A, n, E) = system[i].getParameters();
        double constK = K(A, n, T, E);
        // auto F = [=] (const std::vector<double> &args) -> double {
        //     double right = constK, left = constK;
        //     const auto &in = system[i].getInput();
        //     const auto &out = system[i].getOutput();
        //     for (uint64_t j = 0; j < in.size(); ++j) {
        //         for (uint64_t k = 0; k < in[j].second; ++k) {
        //             right *= args[to_underlying(in[j].first)] * rho;
        //         }
        //     }
        //     for (uint64_t j = 0; i < out.size(); ++j) {
        //         for (uint64_t k = 0; k < out[j].second; ++k) {
        //             left *= args[to_underlying(out[j].first)] * rho;
        //         }
        //     }
        //     return right - left;
        // };
        auto Fin = [=] (const std::vector<double> &args) -> double {
            double ans = constK;
            const auto &in = system[i].getInput();
            for (uint64_t j = 0; j < in.size(); ++j) {
                //for (uint64_t k = 0; k < in[j].second; ++k) {
                    ans *= std::pow(args[in[j].first] * rho, in[j].second);
                //}
            }
            return ans;
        };
        auto Fout = [=] (const std::vector<double> &args) -> double {
            double ans = constK;
            const auto &out = system[i].getInput();
            for (uint64_t j = 0; i < out.size(); ++j) {
                //for (uint64_t k = 0; k < out[j].second; ++k) {
                    ans *= std::pow(args[out[j].first] * rho, out[j].second);
                //}
            }
            return ans;
        };
        //ans.push_back(F);
        rightPart.push_back(Fin);
        rightPart.push_back(Fout);
    }
    return rightPart;
}