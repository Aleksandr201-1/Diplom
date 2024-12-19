#include "ChemicalTask.hpp"

const std::map<ReactionType, std::string> reaction_types = {
    {ReactionType::ISOTERM_CONST_RHO, "IsotermConstRho"},
    {ReactionType::ADIABAT_CONST_RHO, "AdiabatConstRho"},
    {ReactionType::ISOTERM_CONST_P, "IsotermConstP"},
    {ReactionType::ADIABAT_CONST_P, "AdiabatConstP"}
};

std::string reactionTypeToString (ReactionType method) {
    return enumToString(method, reaction_types);
}

ReactionType stringToReactionType (const std::string &str) {
    return stringToEnum(str, reaction_types);
}

double PhiFunction::operator() (double t) const {
    for (uint64_t i = 0; i < T.size(); ++i) {
        if (t >= T[i].first && t < T[i].second) {
            t = t / 10000;
            return phi[i][0] + phi[i][1]*std::log(t) + phi[i][2]*std::pow(t, -2) + phi[i][3]*std::pow(t, -1) + phi[i][4]*t + phi[i][5]*std::pow(t, 2) + phi[i][6]*std::pow(t, 3);
        }
    }
    return 0;
}

double PhiFunction::der1 (double t) const {
    for (uint64_t i = 0; i < T.size(); ++i) {
        if (t >= T[i].first && t < T[i].second) {
            t = t / 10000;
            return (phi[i][1] / t - 2*phi[i][2]*std::pow(t, -3) - phi[i][3]*std::pow(t, -2) + phi[i][4] + 2*phi[i][5]*t + 3*phi[i][6]*std::pow(t, 2)) / 10000;
        }
    }
    return 0;
}

double PhiFunction::der2 (double t) const {
    for (uint64_t i = 0; i < T.size(); ++i) {
        if (t >= T[i].first && t < T[i].second) {
            t = t / 10000;
            return (-phi[i][1] / (t*t) + 6*phi[i][2]*std::pow(t, -4) + 2*phi[i][3]*std::pow(t, -3) + 2*phi[i][5] + 6*phi[i][6]*t) / 100000000;
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

const std::vector<uint64_t> &ChemicalReaction::getInputAdditive () const {
    return input_add;
}

const std::vector<uint64_t> &ChemicalReaction::getOutputAdditive () const {
    return output_add;
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

const double ChemicalSystem::R = 8.3144;
const double ChemicalSystem::P_ATM = 101325;

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


ChemicalSystem::ChemicalSystem () {
    T = P = 0;
}

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

ChemicalSystem::ChemicalSystem (const std::string &filename) : ChemicalSystem() {
    initFromFile(filename);
}

ChemicalSystem::~ChemicalSystem () {}

void ChemicalSystem::setAtomList (const std::vector<std::string> &list) {
    atoms = list;
}

std::vector<std::string> ChemicalSystem::getAtomList () const {
    return atoms;
}

const std::unordered_map<std::string, double> &ChemicalSystem::getAtomMasses () const {
    return atom_mass;
}

void ChemicalSystem::setSubstanceList (const std::vector<std::string> &list) {
    substances = list;
}

std::vector<std::string> ChemicalSystem::getSubstanceList () const {
    return substances;
}

const std::unordered_map<std::string, double> &ChemicalSystem::getSubstanceMasses () const {
    return substance_mass;
}

void ChemicalSystem::addReaction (const std::string &reaction, double A, double n, double E) {
    ChemicalReaction react(reaction, atoms, substances, additives, table);
    react.setParameters(A, n, E);
    reactions.push_back(react);
}

void ChemicalSystem::addReaction (const ChemicalReaction &reaction) {
    reactions.push_back(reaction);
}

void ChemicalSystem::changeReaction (uint64_t i, const std::string &reaction, double A, double n, double E) {
    ChemicalReaction react(reaction, atoms, substances, additives, table);
    react.setParameters(A, n, E);
    reactions[i] = react;
}

void ChemicalSystem::changeReaction (uint64_t i, const ChemicalReaction &reaction) {
    reactions[i] = reaction;
}

void ChemicalSystem::deleteReaction (uint64_t i) {
    reactions.erase(reactions.begin() + i);
}

void ChemicalSystem::clearReactions () {
    reactions.clear();
}

void ChemicalSystem::addAdditive (const std::string &name, const std::vector<double> &additive) {
    if (additive.size() != substances.size()) {
        throw std::runtime_error("ChemicalSystem::addAdditive: additive.size() != substances.size()");
    }
    additives.push_back(name);
    std::cout << "added: " << additives[0] << "\n";
    additive_configs[name] = additive;
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

std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &ChemicalSystem::getTable () {
    return table;
}

void ChemicalSystem::setReactionParameters (uint64_t i, double A, double n, double E) {
    reactions[i].setParameters(A, n, E);
}

std::tuple<double, double, double> ChemicalSystem::getReactionParameters (uint64_t i) const {
    return reactions[i].getParameters();
}

void ChemicalSystem::setConcentrations (const std::vector<double> &conc, ConcentrationMode mode) {
    if (mode == ConcentrationMode::MOLAR_MASS) {
        concentrations = conc;
    } else if (mode == ConcentrationMode::PERCENT) {
        concentrations = conc;
        double rho = 0;
        for (uint64_t j = 0; j < concentrations.size(); ++j) {
            rho += concentrations[j] * substance_mass.find(substances[j])->second;
        }
        for (uint64_t j = 0; j < concentrations.size(); ++j) {
            concentrations[j] = concentrations[j] / rho;
        }
    }
}
// std::vector<double> ChemicalSystem::getConcentrations () const {
//     return concentrations;
// }

double ChemicalSystem::getRho () const {
    auto &gamma = concentrations;
    double rho = 0;
    for (uint64_t j = 0; j < gamma.size(); ++j) {
        rho += gamma[j];// * substance_mass.find(substances[j])->second;
    }
    return P / (R * T * rho);
}
 
std::vector<double> ChemicalSystem::getY0 () const {
    // auto gamma = concentrations;
    // double rho = 0;
    // for (uint64_t j = 0; j < gamma.size(); ++j) {
    //     rho += gamma[j] * substance_mass.find(substances[j])->second;
    // }
    // for (uint64_t j = 0; j < gamma.size(); ++j) {
    //     gamma[j] = gamma[j] / rho;
    // }
    // return gamma;
    return concentrations;
}

const std::vector<std::function<double (double)>> &ChemicalSystem::getGFunc () const {
    return GFunc;
}

const std::vector<PhiFunction> &ChemicalSystem::getPhiFunc () const {
    return phi;
}

const std::vector<std::function<double (double)>> &ChemicalSystem::getEnthalpy () const {
    return enthalpy;
}

void ChemicalSystem::initFromFile (const std::string &filename) {
    std::cout << "init from file " << filename << "\n";
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
        std::pair<double, double> pairH0;
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
        auto gfunc = [=] (double t) -> double {
            return H0[i].first - H0[i].second - t * this->phi[i](t);
        };
        GFunc.push_back(gfunc);
        auto efunc = [=] (double t) -> double {
            return gfunc(t) - t * (-t * phi[i].der1(t) - phi[i](t)); //experimental
            //return gfunc(t) - t * derivative(gfunc, t, 0.01, DiffConfig::POINTS2_ORDER1_WAY3); //normal
        };
        enthalpy.push_back(efunc);
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

    // std::vector<std::function<double (double)>> GFunc;
    // for (uint64_t i = 0; i < system.substances.size(); ++i) {
    //     auto func = [=] (double t) -> double {
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
    auto K = [] (double A, double n, double T, double E) {
        return A * std::pow(T, n) * std::exp(-E / T);
    };
    for (uint64_t i = 0; i < reactions.size(); ++i) {
        double A, n, E;
        std::tie(A, n, E) = reactions[i].getParameters();
        double Kright = K(A, n, T, E);
        double Kleft = 0;
        auto input = reactions[i].getInput();
        auto output = reactions[i].getOutput();
        for (uint64_t j = 0; j < input.size(); ++j) {
            Kleft += (int64_t)(input[j] - output[j]) * (GFunc[j](T) / (R * T) + std::log(R * T / P));
        }
        std::string tmp;
        Kleft = std::exp(-Kleft) * Kright;
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

ChemicalReaction ChemicalSystem::at (uint64_t i) const {
    return reactions[i];
}

ChemicalReaction &ChemicalSystem::at (uint64_t i) {
    return reactions[i];
}

ChemicalReaction ChemicalSystem::operator[] (uint64_t i) const {
    return reactions[i];
}

ChemicalReaction &ChemicalSystem::operator[] (uint64_t i) {
    return reactions[i];
}

std::tuple<std::vector<double>, std::vector<double>> readInfo (const std::string &str) {
    std::vector<double> phi, entalpil;
    return std::make_tuple(phi, entalpil);
}

void ChemicalSystem::rightPartGen () {
    ode_system.clear();
    Y.clear();
    double A, n, E;
    const ChemicalSystem &system = *this;
    auto K = [] (double A, double n, double T, double E) {
        return A * std::pow(T, n) * std::exp(-E / T);
    };
    Wright.clear();
    Wleft.clear();
    Wright.resize(system.getCount());
    Wleft.resize(system.getCount());
    for (uint64_t i = 0; i < system.getCount(); ++i) {
        std::tie(A, n, E) = system[i].getParameters();
        double Kright = K(A, n, T, E);
        double Kleft = 0;
        auto &in = system[i].getInput();
        auto &out = system[i].getOutput();
        for (uint64_t j = 0; j < in.size(); ++j) {
            Kleft += (int64_t)(in[j] - out[j]) * (GFunc[j](T) / (R * T) + std::log(R * T / P_ATM));
        }
        Kleft = std::exp(-Kleft) * Kright;
        //args = t, gamma_1, gamma_2, ..., gamma_n, rho, T
        auto Fin = [=] (const std::vector<double> &args) -> double {
            double Rho = args[args.size() - 2];
            double Temp = args[args.size() - 1];
            double ans = K(A, n, Temp, E);
            double tmp = 0;
            auto &input = this->at(i).getInput();
            auto &input_add = this->at(i).getInputAdditive();
            for (uint64_t j = 0; j < in.size(); ++j) {
                tmp = std::pow(args[j + 1] * Rho, input[j]);
                ans *= tmp;
            }
            for (uint64_t j = 0; j < input_add.size(); ++j) {
                tmp = 0;
                auto &coeff = this->additive_configs.find(this->additives[input_add[j]])->second;
                for (uint64_t k = 0; k < coeff.size(); ++k) {
                    tmp += coeff[k] * args[k + 1];
                }
                ans *= tmp * Rho;
            }
            return ans;
        };

        auto Fout = [=] (const std::vector<double> &args) -> double {
            double Rho = args[args.size() - 2];
            double Temp = args[args.size() - 1];
            double ans = 0;
            auto &input = this->at(i).getInput();
            auto &output = this->at(i).getOutput();
            for (uint64_t j = 0; j < in.size(); ++j) {
                ans += (int64_t)(input[j] - output[j]) * (GFunc[j](Temp) / (R * Temp) + std::log(R * Temp / P_ATM));
            }
            ans = std::exp(-ans) * K(A, n, Temp, E);
            double tmp = 0;
            auto &output_add = this->at(i).getOutputAdditive();
            for (uint64_t j = 0; j < output.size(); ++j) {
                    tmp = std::pow(args[j + 1] * Rho, output[j]);
                    ans *= tmp;
            }
            for (uint64_t j = 0; j < output_add.size(); ++j) {
                tmp = 0;
                auto &coeff = this->additive_configs.find(this->additives[output_add[j]])->second;
                for (uint64_t k = 0; k < coeff.size(); ++k) {
                    tmp += coeff[k] * args[k + 1];
                }
                ans *= tmp * Rho;
            }
            return ans;
        };
        Wright[i] = std::move(Fin);
        Wleft[i] = std::move(Fout);
    }
    for (uint64_t i = 0; i < system.substances.size(); ++i) {
        auto Wsub = [=] (const std::vector<double> &args) -> double {
            double ans = 0;
            double Rho = args[args.size() - 2];
            for (uint64_t j = 0; j < this->Wright.size(); ++j) {
                auto &in = this->at(j).getInput();
                auto &out = this->at(j).getOutput();
                double right = this->Wright[j](args), left = this->Wleft[j](args);
                ans += out[i] * right;
                ans -= in[i]  * right;
                ans += in[i]  * left;
                ans -= out[i] * left;
            }
            return ans / Rho;
        };
        ode_system.push_back(std::move(Wsub));
    }
}