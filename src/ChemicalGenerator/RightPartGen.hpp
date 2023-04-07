#ifndef RIGHT_PART_GEN_HPP
#define RIGHT_PART_GEN_HPP

#include <vector>
#include <cmath>
#include <tuple>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <General/General.hpp>
#include <General/Matrix.hpp>
// #include "../General/General.hpp"
// #include "../General/Matrix.hpp"

template <typename E>
constexpr auto to_underlying(E e) noexcept {
    return static_cast<std::underlying_type_t<E>>(e);
}

// enum class Atom : uint64_t {
//     H,
//     O,
//     COUNT_OF_ATOMS
// };

// enum class Substance : uint64_t {
//     H2O,
//     OH,
//     H2,
//     O2,
//     H,
//     O,
//     COUNT_OF_SUBSTANCE
// };

//Matrix<uint64_t> tableOfSubstance(to_underlying(Atom::COUNT_OF_ATOMS), to_underlying(Substance::COUNT_OF_SUBSTANCE));

struct PhiFunction {
    std::vector<std::vector<double>> phi;
    std::vector<std::pair<double, double>> T;

    double operator() (double t) const;
};

class ChemicalSystem;

class ChemicalReaction {
    private:
        std::string readAtom (const std::string &sub, uint64_t &i) const;
        std::string readSubstance (const std::string &sub, uint64_t &i) const;
        uint64_t readCoeff (const std::string &sub, uint64_t &i) const;
        std::vector<uint64_t> getSubstanceContent (const std::string &sub) const;
        bool checkForCorrect () const;
    public:
        //ChemicalReaction ();
        ChemicalReaction (const std::string &str, const std::vector<std::string> &atoms, const std::vector<std::string> &substances);
        ChemicalReaction (const ChemicalReaction &reaction);
        ChemicalReaction (ChemicalReaction &&reaction);
        ~ChemicalReaction ();

        //static void setAtomList (const std::vector<std::string> &list);
        //static void setSubstanceList (const std::vector<std::string> &list);

        void setReaction (const std::string &str);

        void setInput (const std::vector<uint64_t> &in);
        void setOutput (const std::vector<uint64_t> &out);
        void setParameters (double A, double n, double E);

        const std::vector<uint64_t> &getInput () const;
        const std::vector<uint64_t> &getOutput () const;

        std::tuple<double, double, double> getParameters () const;

        ChemicalReaction &operator= (const ChemicalReaction &reaction);
        ChemicalReaction &operator= (ChemicalReaction &&reaction);

        friend std::ostream &operator<< (std::ostream &out, const ChemicalReaction &reaction);

        friend ChemicalSystem;
    private:
        //static std::vector<std::string> atoms;
        //static std::vector<std::string> substances;
        const std::vector<std::string> &atoms;
        const std::vector<std::string> &substances;
        //std::vector<std::pair<uint64_t, uint64_t>> input, output;
        std::vector<uint64_t> input, output;
        double A, n, E;
};

class ChemicalSystem {
    // private:
    //     std::string readAtom (const std::string &sub, uint64_t &i) const;
    //     std::string readSubstance (const std::string &sub, uint64_t &i) const;
    //     uint64_t readCoeff (const std::string &sub, uint64_t &i) const;
    //     std::vector<uint64_t> getSubstanceContent (const std::string &sub) const;
    //     void initTable ();
    //     bool checkForCorrect () const;
    public:
        ChemicalSystem ();
        ~ChemicalSystem ();

        void setAtomList (const std::vector<std::string> &list);
        void setSubstanceList (const std::vector<std::string> &list);

        void addReaction (const std::string &reaction);
        void addReaction (const ChemicalReaction &reaction);
        void setPressure (double P);
        double getPressure () const;
        void setTemperature (double T);
        double getTemperature () const;
        void setReactionParameters (uint64_t i, double A, double n, double E);
        std::tuple<double, double, double> getReactionParameters (uint64_t i) const;

        void initFromFile (const std::string &filename);

        void printInfo (std::ostream &out) const;

        void getODUSystem () const;

        uint64_t getCount () const;

        std::vector<std::function<double(const std::vector<double> &)>> rightPartGen ();

        ChemicalReaction operator[] (uint64_t i) const;
        ChemicalReaction &operator[] (uint64_t i);

    private:
        std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> table;
        std::vector<std::string> atoms;
        std::unordered_map<std::string, double> atom_mass;
        //std::vector<double> atom_mass;
        std::vector<std::string> substances;
        std::unordered_map<std::string, double> substance_mass;
        //std::vector<double> substance_mass;
        std::vector<std::pair<double, double>> H0;
        //std::vector<std::vector<double>> phi;
        std::vector<PhiFunction> phi;
        std::vector<ChemicalReaction> reactions;
        double T, P;
};

double K (double A, double n, double T, double E);

std::vector<std::function<double(const std::vector<double> &)>> rightPartGen (const ChemicalSystem &system);

std::vector<std::function<double(const std::vector<double> &)>> gammaGen (const ChemicalSystem &system);

#endif