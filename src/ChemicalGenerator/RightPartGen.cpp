#include "RightPartGen.hpp"

double PhiFunction::operator() (double t) const {
    for (uint64_t i = 0; i < T.size(); ++i) {
        if (t >= T[i].first && t <= T[i].second) {
            return phi[i][0] + phi[i][1]*std::log(t) + phi[i][2]*std::pow(t, -2) + phi[i][3]*std::pow(t, -1) + phi[i][4]*t + phi[i][5]*std::pow(t, 2) + phi[i][6]*std::pow(t, 3);
        }
    }
    return 0;
}

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

bool ChemicalReaction::checkForCorrect () const {
    // std::vector<uint64_t> in(atoms.size(), 0), out(atoms.size(), 0);
    // for (uint64_t i = 0; i < atoms.size(); ++i) {
    //     for (uint64_t j = 0; j < input.size(); ++j) {
    //         in[i] += input[j].second * table[atoms[i]][substances[input[j].first]];
    //     }
    //     for (uint64_t j = 0; j < output.size(); ++j) {
    //         out[i] += output[j].second * table[atoms[i]][substances[output[j].first]];
    //     }
    // }
    return input == output;
}

//std::vector<std::pair<Atom, double>> atomMas;

//Matrix<uint64_t> tableOfSubstance(to_underlying(Atom::COUNT_OF_ATOMS), to_underlying(Substance::COUNT_OF_SUBSTANCE));

ChemicalReaction::ChemicalReaction (const std::string &str, const std::vector<std::string> &atoms, const std::vector<std::string> &substances) : atoms(atoms), substances(substances) {
    //initTable();
    // std::cout << "  ";
    // //std::pair<std::string, std::unordered_map<std::string, uint64_t>> &elt = tt;
    // for (auto &elt : tt) {
    //     for (auto &ell : elt.second) {
    //         std::cout << ell.first << "    ";
    //     }
    //     break;
    // }
    
    // std::cout << "\n";
    // for (auto &el : tt) {
    //     std::cout << el.first << " ";
    //     for (auto &ell : el.second) {
    //         std::cout << ell.second << "    ";
    //     }
    //     std::cout << "\n";
    // }
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
    input.resize(substances.size(), 0);
    output.resize(substances.size(), 0);
    std::string tmp = str;
    //std::remove_if(tmp.begin(), tmp.end(), isspace);
    tmp.erase(std::remove_if(tmp.begin(), tmp.end(), isspace), tmp.end());
    std::cout << tmp << "\n";
    //std::find(tmp.begin(), tmp.end(), "<==>");
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
        input[idx] = coeff;
        //input.push_back({idx, coeff});
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
        output[idx] = coeff;
        //output.push_back({idx, coeff});
        ++i;
    }
    std::cout << "==COEFFS INPUT==\n";
    for (uint64_t i = 0; i < substances.size(); ++i) {
        std::cout << substances[i] << ": " << input[i] << "\n";
    }
    std::cout << "==COEFFS OUTPUT==\n";
    for (uint64_t i = 0; i < substances.size(); ++i) {
        std::cout << substances[i] << ": " << output[i] << "\n";
    }
    std::cout << "==CHECK FOR CORRECTNESS\n";
    if (checkForCorrect()) {
        std::cout << "fine\n";
    } else {
        std::cout << "not fine\n";
    }
}

void ChemicalReaction::setInput (const std::vector<uint64_t> &in) {
    input = in;
}

void ChemicalReaction::setOutput (const std::vector<uint64_t> &out) {
    output = out;
}

void ChemicalReaction::setParameters (double A, double n, double E) {
    this->A = A;
    this->n = n;
    this->E = E;
}

const std::vector<uint64_t> &ChemicalReaction::getInput () const {
    return input;
}

const std::vector<uint64_t> &ChemicalReaction::getOutput () const {
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

std::ostream &operator<< (std::ostream &out, const ChemicalReaction &reaction) {
    uint64_t i;
    for (i = 0; i < reaction.input.size(); ++i) {
        if (reaction.input[i] != 0) {
            if (reaction.input[i] > 1) {
                out << reaction.input[i];
            }
            out << reaction.substances[i];
            break;
        }
    }
    for (++i; i < reaction.input.size(); ++i) {
        if (reaction.input[i] > 0) {
            //if (i < reaction.input.size() - 1) {
            out << " + ";
            //}
            if (reaction.input[i] > 1) {
                out << reaction.input[i];
            }
            out << reaction.substances[i];
        }
    }
    out << " <==> ";
    for (i = 0; i < reaction.output.size(); ++i) {
        if (reaction.output[i] != 0) {
            if (reaction.output[i] > 1) {
                out << reaction.output[i];
            }
            out << reaction.substances[i];
            break;
        }
    }
    for (++i; i < reaction.output.size(); ++i) {
        if (reaction.output[i] > 0) {
            //if (i < reaction.input.size() - 1) {
            out << " + ";
            //}
            if (reaction.output[i] > 1) {
                out << reaction.output[i];
            }
            out << reaction.substances[i];
        }
    }
    // for (uint64_t i = 0; i < reaction.output.size(); ++i) {
    //     if (reaction.output[i] > 1) {
    //         out << reaction.output[i];
    //     }
    //     if (reaction.output[i] > 0) {
    //         out << reaction.substances[i];
    //         if (i < reaction.output.size() - 1) {
    //             out << " + ";
    //         }
    //     }
    // }
    return out;
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

void ChemicalSystem::setPressure (double P) {
    this->P = P;
}

double ChemicalSystem::getPressure () const {
    return P;
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

void ChemicalSystem::initFromFile (const std::string &filename) {
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::logic_error("ChemicalSystem::initFromFile: cant open file with name \"" + filename + "\"");
    }
    uint64_t countOfSubstances, countOfAtoms, count;
    std::string str;
    char tmpChar;
    double mass;
    file >> countOfSubstances >> countOfAtoms;
    atoms.resize(countOfAtoms);
    substances.resize(countOfSubstances);
    H0.resize(countOfSubstances);
    phi.resize(countOfSubstances);
    //phi = std::vector<double>(countOfSubstances, std::vector<double>(7));
    for (uint64_t i = 0; i < countOfAtoms; ++i) {
        file >> atoms[i];
    }
    for (uint64_t i = 0; i < countOfAtoms; ++i) {
        file >> mass;
        atom_mass[atoms[i]] = mass;
    }
    for (uint64_t i = 0; i < countOfSubstances; ++i) {
        file >> count >> substances[i] >> mass;
        substance_mass[substances[i]] = mass;
        for (uint64_t j = 0; j < countOfAtoms; ++j) {
            file >> count;
            table[atoms[j]][substances[i]] = count;
        }
        std::pair<double, double> pairH0;
        //std::vector<std::vector<double>> tmpPhi(count, std::vector<double>(7));
        //std::vector<std::pair<double, double>> temp(count);
        std::pair<double, double> temp;
        file >> pairH0.first >> pairH0.second >> count;
        H0[i] = pairH0;
        phi[i].phi.resize(count);
        phi[i].T.resize(count);
        for (uint64_t j = 0; j < count; ++j) {
            phi[i].phi[j].resize(7);
            for (uint64_t k = 0; k < phi[i].phi[0].size(); ++k) {
                file >> phi[i].phi[j][k];
            }
            file >> phi[i].T[j].first >> phi[i].T[j].second;
        }
    }

    file.close();
}

void ChemicalSystem::printInfo (std::ostream &out) const {
    out << "Chemical info\n\n";

    out << "P   = " << P << "Pa\n";
    out << "T   = " << T << "K\n";
    //out << "rho = " << rho << "kg/m3\n\n";
    
    out << "Atoms:\n";
    for (auto const &el : atom_mass) {
        out << el.first << " " << el.second << "\n";
    }
    
    out << "\nSubstances:\n";
    for (auto const &el : substance_mass) {
        out << el.first << " " << el.second << "\n";
    }

    out << "\nInfo about substances:\n";
    for (uint64_t i = 0; i < substances.size(); ++i) {
        out << i + 1 << ") " << substances[i] << "\n";
        out << "delta(H0(T0)) = " << H0[i].first << "\n";
        out << "H0(T0) - H0(0) = " << H0[i].second << "\n";
        out << "Temperature diapasones: " << phi[i].T.size() << "\n";
        for (uint64_t j = 0; j < phi[i].phi.size(); ++j) {
            out << "\tphi:";
            for (uint64_t k = 0; k < phi[i].phi[j].size(); ++k) {
                out << " " << phi[i].phi[j][k];
            }
            out << "\n\tT: " << phi[i].T[j].first << "K ~ " << phi[i].T[j].second << "K\n";
        }
    }

    //вывод констант скоростей при заданной температуре
    out << "\nReactions:\n";
    for (uint64_t i = 0; i < reactions.size(); ++i) {
        double A, n, E;
        std::tie(A, n, E) = reactions[i].getParameters();
        out << i + 1 << ") " << reactions[i] << "\n";
        out << "\tA = " << A << "\n";
        out << "\tn = " << n << "\n";
        out << "\tE = " << E << "\n";
    }
}

void ChemicalSystem::getODUSystem () const {}

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

std::vector<std::function<double(const std::vector<double> &)>> ChemicalSystem::rightPartGen () {
    std::vector<std::function<double(const std::vector<double> &)>> rightPart;
    double A, n, E;
    //double rho, T;
    ChemicalSystem &system = *this;
    //rho = system.getDensity();
    //P = 
    //T = system.getTemperature();
    double R = 8.31;
    // std::vector<double> phi, entalpil;
    // std::tie(phi, entalpil) = readInfo("");
    // double phi1, phi2, phi3, phi4, phi5, phi6, phi7;
    // auto Phi0 = [=] (double t) -> double {
    //     return phi1 + phi2 * std::log(t) + phi3 * std::pow(t, -2) + phi4 * std::pow(t, -1) + phi5 * t + phi6 * std::pow(t, 2) + phi7 * std::pow(t, 3);
    // };
    // double entalpil1, entalpil2;
    // auto G0 = [=] (double t) -> double {
    //     return entalpil1 - entalpil2 - t * Phi0(t);
    // };
    std::vector<std::function<double (double)>> GFunc;
    for (uint64_t i = 0; i < system.substances.size(); ++i) {
        auto func = [=] (double t) -> double {
            return system.H0[i].first + system.H0[i].second + system.phi[i](t);
        };
        GFunc.push_back(func);
    }
    std::cout << "creating reactions speeds\n";
    std::vector<std::function<double (const std::vector<double> &)>> Wright(system.getCount()), Wleft(system.getCount());
    for (uint64_t i = 0; i < system.getCount(); ++i) {
        std::cout << "Reaction " << i + 1 << "\n";
        std::tie(A, n, E) = system[i].getParameters();
        double Kright = K(A, n, T, E);
        double Kleft = 0;
        const auto &in = system[i].getInput();
        const auto &out = system[i].getOutput();
        for (uint64_t j = 0; j < in.size(); ++j) {
            double v = in[j] - out[j];
            v *= GFunc[j](T) / (R * T) + std::log(R * T / P);
            Kleft += v;
        }
        Kleft = std::exp(Kleft);
        std::cout << "K right: " << Kright << "\n";
        std::cout << "K left: " << Kleft << "\n";
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
        // auto W = [=] (const std::vector<double> &args) -> double {
        //     double right = Kright;
        //     double left = Kleft;
        //     const auto &in = system[i].getInput();
        //     const auto &out = system[i].getOutput();
        //     for (uint64_t j = 0; j < in.size(); ++j) {
        //         //for (uint64_t k = 0; k < in[j].second; ++k) {
        //             right *= std::pow(args[in[j].first] * rho, in[j].second);
        //             left *= std::pow(args[out[j].first] * rho, out[j].second);
        //         //}
        //     }
        //     return right - left;
        // };
        auto Fin = [=] (const std::vector<double> &args) -> double {
            double ans = Kright;
            double rho = 0;
            for (uint64_t j = 0; j < in.size(); ++j) {
                rho += args[j + 1];
            }
            rho = P / (R * T * rho);
            //const auto &in = system[i].getInput();
            for (uint64_t j = 0; j < in.size(); ++j) {
                //for (uint64_t k = 0; k < in[j].second; ++k) {
                    ans *= std::pow(args[j + 1] * rho, in[j]);
                //}
            }
            return ans;
        };
        auto Fout = [=] (const std::vector<double> &args) -> double {
            double ans = Kleft;
            double rho = 0;
            for (uint64_t j = 0; j < in.size(); ++j) {
                rho += args[j + 1];
            }
            rho = P / (R * T * rho);
            //const auto &out = system[i].getOutput();
            for (uint64_t j = 0; i < out.size(); ++j) {
                //for (uint64_t k = 0; k < out[j].second; ++k) {
                    ans *= std::pow(args[j + 1] * rho, out[j]);
                //}
            }
            return ans;
        };
        Wright[i] = std::move(Fin);
        Wleft[i] = std::move(Fout);
        //ans.push_back(F);
        //rightPart.push_back(Fin);
        //rightPart.push_back(Fout);
    }
    std::cout << "creating substances speeds\n";
    for (uint64_t i = 0; i < system.substances.size(); ++i) {
        auto Wsub = [=] (const std::vector<double> &args) -> double {
            double ans = 0;
            double rho = 0;
            for (uint64_t j = 0; j < args.size(); ++j) {
                rho += args[j + 1];
            }
            for (uint64_t j = 0; j < Wright.size(); ++j) {
                ans += system[j].getInput()[i] * Wright[j](args);
            }
            for (uint64_t j = 0; j < Wleft.size(); ++j) {
                ans -= system[j].getOutput()[i] * Wleft[j](args);
            }
            return ans / rho;
        };
        rightPart.push_back(Wsub);
    }
    std::cout << "Chem size: " << rightPart.size() << "\n";
    std::cout << "Chemic func exmpl: " << rightPart[0]({0, 0.5, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0}) << "\n";

    std::cout << "done\n";
    return rightPart;
}