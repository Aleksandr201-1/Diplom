#include "ChemicalTask.hpp"

float128_t PhiFunction::operator() (float128_t t) const {
    for (uint64_t i = 0; i < T.size(); ++i) {
        if (t >= T[i].first && t < T[i].second) {
            t = t / 10000;
            return phi[i][0] + phi[i][1]*std::log(t) + phi[i][2]*std::pow(t, -2) + phi[i][3]*std::pow(t, -1) + phi[i][4]*t + phi[i][5]*std::pow(t, 2) + phi[i][6]*std::pow(t, 3);
        }
    }
    return 0;
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

bool ChemicalReaction::checkForCorrect () const {
    std::vector<uint64_t> in(atoms.size(), 0), out(atoms.size(), 0);
    for (uint64_t i = 0; i < atoms.size(); ++i) {
        for (uint64_t j = 0; j < input.size(); ++j) {
            in[i] += input[j] * table.at(atoms[i]).at(substances[j]);;
        }
        for (uint64_t j = 0; j < output.size(); ++j) {
            out[i] += output[j] * table.at(atoms[i]).at(substances[j]);;
        }
    }
    return input == output;
}

ChemicalReaction::ChemicalReaction (const std::string &str,
                                    const std::vector<std::string> &atoms,
                                    const std::vector<std::string> &substances,
                                    const std::vector<std::string> &additives,
                                    const std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &table) : atoms(atoms), substances(substances), additives(additives), table(table) {
    setReaction(str);
}

ChemicalReaction::ChemicalReaction (const ChemicalReaction &reaction) : atoms(reaction.atoms), substances(reaction.substances), additives(reaction.additives), table(reaction.table) {
    input = reaction.input;
    output = reaction.output;
    input_add = reaction.input_add;
    output_add = reaction.output_add;
    A = reaction.A;
    n = reaction.n;
    E = reaction.E;
}

ChemicalReaction::ChemicalReaction (ChemicalReaction &&reaction) : atoms(std::move(reaction.atoms)), substances(std::move(reaction.substances)), additives(std::move(additives)), table(std::move(reaction.table)) {
    input = std::move(reaction.input);
    output = std::move(reaction.output);
    input_add = std::move(reaction.input_add);
    output_add = std::move(reaction.output_add);
    A = reaction.A;
    n = reaction.n;
    E = reaction.E;
}

ChemicalReaction::~ChemicalReaction () {}

void ChemicalReaction::setReaction (const std::string &str) {
    input.resize(substances.size(), 0);
    output.resize(substances.size(), 0);
    std::string tmp = str;
    //std::remove_if(tmp.begin(), tmp.end(), isspace);
    tmp.erase(std::remove_if(tmp.begin(), tmp.end(), isspace), tmp.end());
    std::cout << tmp << "\n";
    //std::find(tmp.begin(), tmp.end(), "<==>");
    uint64_t div = tmp.find('<');
    std::string inputStr = tmp.substr(0, div);
    div = tmp.find('>');
    std::string outputStr = tmp.substr(div + 1, tmp.size() - div);
    std::cout << inputStr << "||" << outputStr << "\n";
    uint64_t i = 0;
    //std::string &str = {inputStr, outputStr};
    while (i < inputStr.size()) {
        uint64_t coeff = readCoeff(inputStr, i);
        std::string sub = readSubstance(inputStr, i);
        auto it = std::find(substances.cbegin(), substances.cend(), sub);
        uint64_t idx = std::distance(substances.cbegin(), it);
        if (idx == substances.size()) {
            it = std::find(additives.cbegin(), additives.cend(), sub);
            idx = std::distance(additives.cbegin(), it);
            if (idx == additives.size()) {
                throw std::logic_error("operation \"setReaction\": substance \"" + sub + "\" not presented in list");
            } else {
                input_add.push_back(idx);
            }
        } else {
            input[idx] += coeff;
        }
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
            it = std::find(additives.cbegin(), additives.cend(), sub);
            idx = std::distance(additives.cbegin(), it);
            if (idx == additives.size()) {
                throw std::logic_error("operation \"setReaction\": substance \"" + sub + "\" not presented in list");
            } else {
                output_add.push_back(idx);
            }
        } else {
            output[idx] += coeff;
        }
        //output.push_back({idx, coeff});
        ++i;
    }
    //std::cout << "==COEFFS INPUT==\n";
    for (uint64_t i = 0; i < substances.size(); ++i) {
        std::cout << substances[i] << ": " << input[i] << "\n";
    }
    //std::cout << "==COEFFS OUTPUT==\n";
    for (uint64_t i = 0; i < substances.size(); ++i) {
        std::cout << substances[i] << ": " << output[i] << "\n";
    }
    //std::cout << "==CHECK FOR CORRECTNESS\n";
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

void ChemicalReaction::setParameters (float128_t A, float128_t n, float128_t E) {
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

const std::vector<uint64_t> &ChemicalReaction::getInputAdditive () const {
    return input_add;
}

const std::vector<uint64_t> &ChemicalReaction::getOutputAdditive () const {
    return output_add;
}

std::tuple<float128_t, float128_t, float128_t> ChemicalReaction::getParameters () const {
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
        //std::cout << reaction.input_add.size() << "\n";
        //if (reaction.input_add.size() > 0) {
            //std::cout << reaction.input_add[0] << "\n";
            //std::cout << reaction.additives.size() << "\n";
            //std::cout << reaction.additives.size() << "\n";
            //std::cout << reaction.additives[reaction.input_add[0]] << "\n";
        //}
    for (uint64_t j = 0; j < reaction.input_add.size(); ++j) {
        out << " + ";
        out << reaction.additives[reaction.input_add[j]];
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
    for (uint64_t j = 0; j < reaction.output_add.size(); ++j) {
        out << " + ";
        out << reaction.additives[reaction.output_add[j]];
    }
    return out;
}

const float128_t ChemicalSystem::R = 8.3144;
const float128_t ChemicalSystem::P0 = 101325;

std::string ChemicalSystem::readAtom (const std::string &sub, uint64_t &i) const {
    std::string ans;
    ans += sub[i];
    ++i;
    if (i < sub.size() && sub[i] >= 'a' && sub[i] <= 'z') {
        ans += sub[i];
        ++i;
    }
    return ans;
}

uint64_t ChemicalSystem::readCoeff (const std::string &sub, uint64_t &i) const {
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

std::vector<uint64_t> ChemicalSystem::getSubstanceContent (const std::string &sub) const {
    std::vector<uint64_t> ans(atoms.size(), 0);
    uint64_t i = 0;
    while (i < sub.size()) {
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

void ChemicalSystem::initTable () {
    for (uint64_t i = 0; i < substances.size(); ++i) {
        auto coeffs = getSubstanceContent(substances[i]);
        for (uint64_t j = 0; j < atoms.size(); ++j) {
            table[atoms[j]][substances[i]] = coeffs[j];
        }
    }
}


ChemicalSystem::ChemicalSystem () {}

ChemicalSystem::ChemicalSystem (const ChemicalSystem &system) {
    table = system.table;
    atoms = system.atoms;
    atom_mass = system.atom_mass;
    substances = system.substances;
    substance_mass = system.substance_mass;
    additives = system.additives;
    additive_configs = system.additive_configs;
    H0 = system.H0;
    GFunc = system.GFunc;
    concentrations = system.concentrations;
    phi = system.phi;
    reactions = system.reactions;
    T = system.T;
    P = system.P;
}

ChemicalSystem::ChemicalSystem (ChemicalSystem &&system) {
    table = std::move(system.table);
    atoms = std::move(system.atoms);
    atom_mass = std::move(system.atom_mass);
    substances = std::move(system.substances);
    substance_mass = std::move(system.substance_mass);
    additives = std::move(system.additives);
    additive_configs = std::move(system.additive_configs);
    H0 = std::move(system.H0);
    GFunc = std::move(system.GFunc);
    concentrations = std::move(system.concentrations);
    phi = std::move(system.phi);
    reactions = std::move(system.reactions);
    T = std::move(system.T);
    P = std::move(system.P);
}

ChemicalSystem::ChemicalSystem (const std::string &filename) {
    initFromFile(filename);
}

ChemicalSystem::~ChemicalSystem () {}

void ChemicalSystem::setAtomList (const std::vector<std::string> &list) {
    atoms = list;
}

std::vector<std::string> ChemicalSystem::getAtomList () const {
    return atoms;
}

void ChemicalSystem::setSubstanceList (const std::vector<std::string> &list) {
    substances = list;
}

std::vector<std::string> ChemicalSystem::getSubstanceList () const {
    return substances;
}

void ChemicalSystem::addReaction (const std::string &reaction, float128_t A, float128_t n, float128_t E) {
    ChemicalReaction react(reaction, atoms, substances, additives, table);
    react.setParameters(A, n, E);
    reactions.push_back(react);
}

void ChemicalSystem::addReaction (const ChemicalReaction &reaction) {
    reactions.push_back(reaction);
}

void ChemicalSystem::addAdditive (const std::string &name, const std::vector<float128_t> &additive) {
    if (additive.size() != substances.size()) {
        throw std::runtime_error("ChemicalSystem::addAdditive: additive.size() != substances.size()");
    }
    additives.push_back(name);
    std::cout << "added: " << additives[0] << "\n";
    additive_configs[name] = additive;
}

void ChemicalSystem::setPressure (float128_t P) {
    this->P = P;
}

float128_t ChemicalSystem::getPressure () const {
    return P;
}

void ChemicalSystem::setTemperature (float128_t T) {
    this->T = T;
}

float128_t ChemicalSystem::getTemperature () const {
    return T;
}

std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &ChemicalSystem::getTable () {
    return table;
}

void ChemicalSystem::setReactionParameters (uint64_t i, float128_t A, float128_t n, float128_t E) {
    reactions[i].setParameters(A, n, E);
}

std::tuple<float128_t, float128_t, float128_t> ChemicalSystem::getReactionParameters (uint64_t i) const {
    return reactions[i].getParameters();
}

void ChemicalSystem::setConcentrations (const std::vector<float128_t> &conc) {
    concentrations = conc;
}

std::vector<float128_t> ChemicalSystem::getY0 () const {
    auto gamma = concentrations;
    float128_t rho = 0;
    for (uint64_t j = 0; j < gamma.size(); ++j) {
        rho += gamma[j] * substance_mass.find(substances[j])->second;
    }
    for (uint64_t j = 0; j < gamma.size(); ++j) {
        gamma[j] = gamma[j] / rho;
    }
    return gamma;
}

const std::vector<std::function<float128_t (float128_t)>> &ChemicalSystem::getGFunc () const {
    return GFunc;
}

void ChemicalSystem::initFromFile (const std::string &filename) {
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::logic_error("ChemicalSystem::initFromFile: cant open file with name \"" + filename + "\"");
    }
    uint64_t countOfSubstances, countOfAtoms, count;
    std::string str;
    char tmpChar;
    float128_t mass;
    file >> countOfSubstances >> countOfAtoms;
    atoms.resize(countOfAtoms);
    substances.resize(countOfSubstances);
    H0.resize(countOfSubstances);
    phi.resize(countOfSubstances);
    for (uint64_t i = 0; i < countOfAtoms; ++i) {
        file >> atoms[i];
    }
    for (uint64_t i = 0; i < countOfAtoms; ++i) {
        file >> mass;
        atom_mass[atoms[i]] = mass / 1000;
    }
    for (uint64_t i = 0; i < countOfSubstances; ++i) {
        file >> count >> substances[i] >> mass;
        substance_mass[substances[i]] = mass / 1000;
        for (uint64_t j = 0; j < countOfAtoms; ++j) {
            file >> count;
            table[atoms[j]][substances[i]] = count;
        }
        std::pair<float128_t, float128_t> pairH0;
        file >> pairH0.first >> pairH0.second >> count;
        pairH0.first  *= 1000;
        pairH0.second *= 1000;
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
    for (uint64_t i = 0; i < substances.size(); ++i) {
        auto func = [=] (float128_t t) -> float128_t {
            return H0[i].first - H0[i].second - t * phi[i](t);
        };
        GFunc.push_back(func);
    }

    file.close();
}

void ChemicalSystem::printInfo (std::ostream &out) const {
    const ChemicalSystem &system = *this;
    out << "Chemical info\n\n";

    out << "P = " << P << "Pa\n"
           "T = " << T << "K\n"
           "R = " << R << "\n\n";
    
    out << "Atoms:\n";
    for (const auto &el : atom_mass) {
        out << std::setw(5) << el.first << "  " << el.second << "\n";
    }
    
    out << "\nSubstances:\n";
    for (const auto &el : substance_mass) {
        out << std::setw(5) << el.first << "  " << el.second << "\n";
    }

    // std::vector<std::function<float128_t (float128_t)>> GFunc;
    // for (uint64_t i = 0; i < system.substances.size(); ++i) {
    //     auto func = [=] (float128_t t) -> float128_t {
    //         return system.H0[i].first - system.H0[i].second - t * system.phi[i](t);
    //     };
    //     GFunc.push_back(func);
    // }

    out << "\nInfo about substances:\n";
    for (uint64_t i = 0; i < substances.size(); ++i) {
        out << i + 1 << ") " << substances[i] << "\n";
        out << "delta(H0(T0))  = " << H0[i].first << "\n";
        out << "H0(T0) - H0(0) = " << H0[i].second << "\n";
        out << "Temperature diapasones: " << phi[i].T.size() << "\n";
        for (uint64_t j = 0; j < phi[i].phi.size(); ++j) {
            out << "\tphi:";
            for (uint64_t k = 0; k < phi[i].phi[j].size(); ++k) {
                out << " " << phi[i].phi[j][k];
            }
            out << "\n\tT: " << phi[i].T[j].first << "K ~ " << phi[i].T[j].second << "K\n";
        }
        out << "\tF0(" << T << ") = " << phi[i](T) << "\n";
        out << "\tG(" << T << ")  = " << GFunc[i](T) << "\n";
    }

    out << "\nAdditives:\n";
    for (const auto &el : additive_configs) {
        out << el.first << ":";
        for (auto num : el.second) {
            out << " " << num;
        }
        out << "\n";
    }

    //вывод констант скоростей при заданной температуре
    out << "\nReactions:\n";
    auto K = [] (float128_t A, float128_t n, float128_t T, float128_t E) {
        return A * std::pow(T, n) * std::exp(-E / T);
    };
    for (uint64_t i = 0; i < reactions.size(); ++i) {
        float128_t A, n, E;
        std::tie(A, n, E) = reactions[i].getParameters();
        float128_t Kright = K(A, n, T, E);
        float128_t Kleft = 0;
        auto input = reactions[i].getInput();
        auto output = reactions[i].getOutput();
        for (uint64_t j = 0; j < input.size(); ++j) {
            Kleft += (int64_t)(input[j] - output[j]) * (GFunc[j](T) / (R * T) + std::log(R * T / P));
        }
        std::string tmp;
        Kleft = std::exp(Kleft) * Kright;
        out << i + 1 << ") " << reactions[i] << "\n";
        out << std::left;
        tmp = "\tA = " + std::to_string(A);
        out << std::setw(20) << tmp;
        out << "\tKright(" << T << ") = " << Kright << "\n";
        tmp = "\tn = " + std::to_string(n);
        out << std::setw(20) << tmp;
        out << "\tKleft(" << T << ")  = " << Kleft << "\n";
        out << "\tE = " << E << "\n";
    }
}

uint64_t ChemicalSystem::getCount () const {
    return reactions.size();
}

ChemicalSystem &ChemicalSystem::operator= (const ChemicalSystem &system) {
    this->table = system.table;
    this->atoms = system.atoms;
    this->atom_mass = system.atom_mass;
    this->substances = system.substances;
    this->substance_mass = system.substance_mass;
    this->additives = system.additives;
    this->additive_configs = system.additive_configs;
    this->H0 = system.H0;
    this->GFunc = system.GFunc;
    this->concentrations = system.concentrations;
    this->phi = system.phi;
    this->reactions = system.reactions;
    this->T = system.T;
    this->P = system.P;
    return *this;
}

ChemicalSystem &ChemicalSystem::operator= (ChemicalSystem &&system) {
    this->table = std::move(system.table);
    this->atoms = std::move(system.atoms);
    this->atom_mass = std::move(system.atom_mass);
    this->substances = std::move(system.substances);
    this->substance_mass = std::move(system.substance_mass);
    this->additives = std::move(system.additives);
    this->additive_configs = std::move(system.additive_configs);
    this->H0 = std::move(system.H0);
    this->GFunc = std::move(system.GFunc);
    this->concentrations = std::move(system.concentrations);
    this->phi = std::move(system.phi);
    this->reactions = std::move(system.reactions);
    this->T = std::move(system.T);
    this->P = std::move(system.P);
    return *this;
}

ChemicalReaction ChemicalSystem::operator[] (uint64_t i) const {
    return reactions[i];
}

ChemicalReaction &ChemicalSystem::operator[] (uint64_t i) {
    return reactions[i];
}

std::tuple<std::vector<float128_t>, std::vector<float128_t>> readInfo (const std::string &str) {
    std::vector<float128_t> phi, entalpil;
    return std::make_tuple(phi, entalpil);
}

void ChemicalSystem::rightPartGen () {
    //std::vector<std::function<float128_t(const std::vector<float128_t> &)>> rightPart;
    float128_t A, n, E;
    const ChemicalSystem &system = *this;
    auto K = [] (float128_t A, float128_t n, float128_t T, float128_t E) {
        return A * std::pow(T, n) * std::exp(-E / T);
    };
    // std::vector<std::function<float128_t (float128_t)>> GFunc;
    // for (uint64_t i = 0; i < system.substances.size(); ++i) {
    //     auto func = [=] (float128_t t) -> float128_t {
    //         return system.H0[i].first - system.H0[i].second - t * system.phi[i](t);
    //     };
    //     GFunc.push_back(func);
    // }
    // std::cout << "GFunc:\n";
    // for (uint64_t i = 0; i < GFunc.size(); ++i) {
    //     std::cout << "func " << i + 1 << ": " << GFunc[i](T) << "\n";
    // }
    std::cout << "creating reactions speeds\n";
    std::vector<std::function<float128_t (const std::vector<float128_t> &)>> Wright(system.getCount()), Wleft(system.getCount());
    for (uint64_t i = 0; i < system.getCount(); ++i) {
        std::cout << "\nReaction " << i + 1 << "\n";
        std::tie(A, n, E) = system[i].getParameters();
        float128_t Kright = K(A, n, T, E);
        float128_t Kleft = 0;
        auto in = system[i].getInput();
        auto out = system[i].getOutput();
        //const float128_t P0 = 101'325;
        for (uint64_t j = 0; j < in.size(); ++j) {
            Kleft += (int64_t)(in[j] - out[j]) * (GFunc[j](T) / (R * T) + std::log(R * T / P0));
            //std::cout << "curr Kleft: " << Kleft << "\n";
            //std::cout << "for " << substances[j] << ": " << (int64_t)(in[j] - out[j]) << "\n";
        }
        //std::cout << "SUM: " << (GFunc[0](T) - 2*GFunc[1](T) + GFunc[3](T)) / (R * T) << "\n";
        Kleft = std::exp(-Kleft) * Kright;
        //std::cout << "K right: " << Kright << "\n";
        //std::cout << "K left: " << Kleft << "\n";
        //exit(0);
        //args = t, gamma_1, gamma_2, ..., gamma_n, rho, T
        auto Fin = [=] (const std::vector<float128_t> &args) -> float128_t {
            float128_t Rho = args[args.size() - 2];
            float128_t Temp = args[args.size() - 1];
            float128_t ans = K(A, n, Temp, E);
            float128_t tmp;
            //std::cout << "Kin: " << ans << "\n";
            auto input = system[i].getInput();
            auto input_add = system[i].getInputAdditive();
            for (uint64_t j = 0; j < in.size(); ++j) {
                tmp = std::pow(args[j + 1] * Rho, input[j]);
                // if (tmp != 0.0) {
                //     ans *= tmp;
                // }
                ans *= tmp;
            }
            for (uint64_t j = 0; j < input_add.size(); ++j) {
                tmp = 0;
                auto coeff = additive_configs.find(additives[input_add[j]])->second;
                for (uint64_t k = 0; k < coeff.size(); ++k) {
                    //std::cout << k << ": " << coeff[k] << " * " << args[k + 1] << "\n"; 
                    tmp += coeff[k] * args[k + 1];
                }
                //std::cout << "in M: " << tmp * Rho << "\n";
                //std::cout << "sum = " << tmp << "\n";
                ans *= tmp * Rho;
                //exit(0);
            }
            //exit(0);
            return ans;
        };

        auto Fout = [=] (const std::vector<float128_t> &args) -> float128_t {
            float128_t Rho = args[args.size() - 2];
            float128_t Temp = args[args.size() - 1];
            float128_t ans = 0;
            auto input = system[i].getInput();
            auto output = system[i].getOutput();
            for (uint64_t j = 0; j < in.size(); ++j) {
                ans += (int64_t)(input[j] - output[j]) * (GFunc[j](Temp) / (R * Temp) + std::log(R * Temp / P0));
            }
            ans = std::exp(-ans) * K(A, n, T, E);
            //std::cout << "Kout: " << ans << "\n";
            float128_t tmp;
            auto output_add = system[i].getOutputAdditive();
            for (uint64_t j = 0; j < output.size(); ++j) {
                //for (uint64_t k = 0; k < out[j].second; ++k) {
                    tmp = std::pow(args[j + 1] * Rho, output[j]);
                    // if (tmp != 0.0) {
                    //     ans *= tmp;
                    // }
                    ans *= tmp;
                    //ans *= std::pow(args[j + 1] * rho, out[j]);
                //}
            }
            for (uint64_t j = 0; j < output_add.size(); ++j) {
                tmp = 0;
                auto coeff = additive_configs.find(additives[output_add[j]])->second;
                for (uint64_t k = 0; k < coeff.size(); ++k) {
                    tmp += coeff[k] * args[k + 1];
                }
                //std::cout << "out M: " << tmp * Rho << "\n";
                ans *= tmp * Rho;
            }
            return ans;
        };
        Wright[i] = std::move(Fin);
        Wleft[i] = std::move(Fout);
    }
    std::cout << "Wright size: " << Wright.size() << "\n";
    std::cout << "Wleft size: " << Wleft.size() << "\n";
    std::cout << "creating substances speeds\n";
    for (uint64_t i = 0; i < system.substances.size(); ++i) {
        auto Wsub = [=] (const std::vector<float128_t> &args) -> float128_t {
            float128_t ans = 0;
            float128_t Rho = args[args.size() - 2];
            // for (uint64_t j = 1; j < args.size(); ++j) {
            //     rho += args[j];
            // }
            // rho = P / (R * T * rho);
            //std::cout << "www\n";
            //std::cout << "REACTION " << i + 1 << "\n";
            for (uint64_t j = 0; j < Wright.size(); ++j) {
                auto in = system[j].getInput();
                auto out = system[j].getOutput();
                float128_t right = Wright[j](args), left = Wleft[j](args);
                ans += out[i] * right;
                ans -= in[i]  * right;
                ans += in[i]  * left;
                ans -= out[i] * left;
                auto input_add = system[j].getInputAdditive();
                auto output_add = system[j].getOutputAdditive();
                for (uint64_t k = 0; k < input_add.size(); ++k) {
                    //std::cout << "got M in reaction " << j + 1 << "to work\n";
                    ans += additive_configs[additives[output_add[k]]][i] * right;
                    ans -= additive_configs[additives[input_add[k]]][i]  * right;
                    ans += additive_configs[additives[input_add[k]]][i]  * left;
                    ans -= additive_configs[additives[output_add[k]]][i] * left;
                }
            }
            //37 000
            //
            //std::cout << "W" << i << ": " << ans / Rho << "\n";
            //exit(0);
            return ans / Rho;
        };
        ode_system.push_back(Wsub);
    }
    //std::cout << "Chem size: " << rightPart.size() << "\n";
    //std::cout << "Chemic func exmpl: " << rightPart[0]({0, 0.5, 0, 0, 0.5, 0, 0}) << "\n";

    //std::cout << "done\n";
    //return rightPart;
}