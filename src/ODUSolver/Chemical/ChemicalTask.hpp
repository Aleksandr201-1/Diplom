#ifndef CHEMICAL_TASK_HPP
#define CHEMICAL_TASK_HPP

#include <vector>
#include <cmath>
#include <tuple>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <Math/Matrix.hpp>
#include <NumericMethods/Differentiation.hpp>
#include <ODUSolver/Task.hpp>

enum class ReactionType {
    ISOTERM_CONST_RHO,
    ADIABAT_CONST_RHO,
    ISOTERM_CONST_P,
    ADIABAT_CONST_P,
    ERROR
};

std::string reactionTypeToString (ReactionType method);

ReactionType stringToReactionType (const std::string &str);

enum class ConcentrationMode {
    PERCENT,
    MOLAR_MASS
};

struct PhiFunction {
    std::vector<std::vector<double>> phi;
    std::vector<std::pair<double, double>> T;

    double operator() (double t) const;

    double der1 (double t) const;

    double der2 (double t) const;
};

class ChemicalSystem;

class ChemicalReaction {
    private:
        //std::string readAtom (const std::string &sub, uint64_t &i) const;
        std::string readSubstance (const std::string &sub, uint64_t &i) const;
        uint64_t readCoeff (const std::string &sub, uint64_t &i) const;
        //std::vector<uint64_t> getSubstanceContent (const std::string &sub) const;
        bool checkForCorrect () const;
    public:
        ChemicalReaction (const std::string &str, const std::vector<std::string> &atoms,
                          const std::vector<std::string> &substances,
                          const std::vector<std::string> &additives,
                          const std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &table);
        ChemicalReaction (const ChemicalReaction &reaction);
        ChemicalReaction (ChemicalReaction &&reaction);
        ~ChemicalReaction ();

        void setReaction (const std::string &str);

        void setInput (const std::vector<uint64_t> &in);
        void setOutput (const std::vector<uint64_t> &out);
        void setParameters (double A, double n, double E);

        const std::vector<uint64_t> &getInput () const;
        const std::vector<uint64_t> &getOutput () const;
        const std::vector<uint64_t> &getInputAdditive () const;
        const std::vector<uint64_t> &getOutputAdditive () const;

        std::tuple<double, double, double> getParameters () const;

        ChemicalReaction &operator= (const ChemicalReaction &reaction);
        ChemicalReaction &operator= (ChemicalReaction &&reaction);

        friend std::ostream &operator<< (std::ostream &out, const ChemicalReaction &reaction);

        friend ChemicalSystem;
    private:
        const std::vector<std::string> &atoms;
        const std::vector<std::string> &substances;
        const std::vector<std::string> &additives;
        const std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &table;
        //const std::unordered_map<std::string, std::vector<double>> &additives;
        std::vector<uint64_t> input, output;
        std::vector<uint64_t> input_add, output_add;
        //std::vector<std::iterator<double>> ff;
        double A, n, E;
};

class ChemicalSystem : public Task {
    private:
        std::string readAtom (const std::string &sub, uint64_t &i) const;
        //std::string readSubstance (const std::string &sub, uint64_t &i) const;
        uint64_t readCoeff (const std::string &sub, uint64_t &i) const;
        std::vector<uint64_t> getSubstanceContent (const std::string &sub) const;
        void initTable ();
    //     bool checkForCorrect () const;
    public:
        ChemicalSystem ();
        ChemicalSystem (const ChemicalSystem &system);
        ChemicalSystem (ChemicalSystem &&system);
        ChemicalSystem (const std::string &filename);
        ~ChemicalSystem ();

        void setAtomList (const std::vector<std::string> &list);
        std::vector<std::string> getAtomList () const;
        const std::unordered_map<std::string, double> &getAtomMasses () const;
        void setSubstanceList (const std::vector<std::string> &list);
        std::vector<std::string> getSubstanceList () const;
        const std::unordered_map<std::string, double> &getSubstanceMasses () const;

        void addReaction (const std::string &reaction, double A, double n, double E);
        void addReaction (const ChemicalReaction &reaction);
        void changeReaction (uint64_t i, const std::string &reaction, double A, double n, double E);
        void changeReaction (uint64_t i, const ChemicalReaction &reaction);
        void deleteReaction (uint64_t i);
        void clearReactions ();

        void addAdditive (const std::string &name, const std::vector<double> &additive);
        void setPressure (double P);
        double getPressure () const;
        void setTemperature (double T);
        double getTemperature () const;
        std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &getTable ();
        void setReactionParameters (uint64_t i, double A, double n, double E);
        std::tuple<double, double, double> getReactionParameters (uint64_t i) const;
        void setConcentrations (const std::vector<double> &conc, ConcentrationMode mode);
        //std::vector<double> getConcentrations () const;
        double getRho () const;
        std::vector<double> getY0 () const;
        const std::vector<std::function<double (double)>> &getGFunc () const;
        const std::vector<PhiFunction> &getPhiFunc () const;
        const std::vector<std::function<double (double)>> &getEnthalpy () const;


        void initFromFile (const std::string &filename);

        void printInfo (std::ostream &out) const;

        //ChemicalTask getChemicalTask () const;

        uint64_t getCount () const;

        //std::vector<std::function<double(const std::vector<double> &)>> rightPartGen () const;
        void rightPartGen ();

        ChemicalReaction at (uint64_t i) const;
        ChemicalReaction &at (uint64_t i);

        ChemicalSystem &operator= (const ChemicalSystem &reaction);
        ChemicalSystem &operator= (ChemicalSystem &&reaction);

        ChemicalReaction operator[] (uint64_t i) const;
        ChemicalReaction &operator[] (uint64_t i);

        static const double R, P_ATM;

    private:
        std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> table;
        std::vector<std::string> atoms;
        std::unordered_map<std::string, double> atom_mass;
        std::vector<std::string> substances;
        std::unordered_map<std::string, double> substance_mass;
        std::vector<std::string> additives;
        std::unordered_map<std::string, std::vector<double>> additive_configs;
        std::vector<std::pair<double, double>> H0;
        std::vector<std::function<double (double)>> GFunc, enthalpy;
        std::vector<std::function<double (const std::vector<double> &)>> Wright, Wleft;
        std::vector<double> concentrations;
        std::vector<PhiFunction> phi;
        std::vector<ChemicalReaction> reactions;
        double T, P;
};

#endif