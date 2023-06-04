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
#include <General/General.hpp>
#include <General/Matrix.hpp>
#include <ODUSolver/Task.hpp>

enum class ReactionType {
    ISOTERM_CONST_RHO,
    ADIABAT_CONST_RHO,
    ISOTERM_CONST_P,
    ADIABAT_CONST_P
};

struct PhiFunction {
    std::vector<std::vector<float128_t>> phi;
    std::vector<std::pair<float128_t, float128_t>> T;

    float128_t operator() (float128_t t) const;

    float128_t der1 (float128_t t) const;

    float128_t der2 (float128_t t) const;
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
        void setParameters (float128_t A, float128_t n, float128_t E);

        const std::vector<uint64_t> &getInput () const;
        const std::vector<uint64_t> &getOutput () const;
        const std::vector<uint64_t> &getInputAdditive () const;
        const std::vector<uint64_t> &getOutputAdditive () const;

        std::tuple<float128_t, float128_t, float128_t> getParameters () const;

        ChemicalReaction &operator= (const ChemicalReaction &reaction);
        ChemicalReaction &operator= (ChemicalReaction &&reaction);

        friend std::ostream &operator<< (std::ostream &out, const ChemicalReaction &reaction);

        friend ChemicalSystem;
    private:
        const std::vector<std::string> &atoms;
        const std::vector<std::string> &substances;
        const std::vector<std::string> &additives;
        const std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &table;
        //const std::unordered_map<std::string, std::vector<float128_t>> &additives;
        std::vector<uint64_t> input, output;
        std::vector<uint64_t> input_add, output_add;
        //std::vector<std::iterator<float128_t>> ff;
        float128_t A, n, E;
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
        const std::unordered_map<std::string, float128_t> &getAtomMasses () const;
        void setSubstanceList (const std::vector<std::string> &list);
        std::vector<std::string> getSubstanceList () const;
        const std::unordered_map<std::string, float128_t> &getSubstanceMasses () const;

        void addReaction (const std::string &reaction, float128_t A, float128_t n, float128_t E);
        void addReaction (const ChemicalReaction &reaction);
        void changeReaction (uint64_t i, const std::string &reaction, float128_t A, float128_t n, float128_t E);
        void changeReaction (uint64_t i, const ChemicalReaction &reaction);
        void deleteReaction (uint64_t i);
        void clearReactions ();

        void addAdditive (const std::string &name, const std::vector<float128_t> &additive);
        void setPressure (float128_t P);
        float128_t getPressure () const;
        void setTemperature (float128_t T);
        float128_t getTemperature () const;
        std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> &getTable ();
        void setReactionParameters (uint64_t i, float128_t A, float128_t n, float128_t E);
        std::tuple<float128_t, float128_t, float128_t> getReactionParameters (uint64_t i) const;
        void setConcentrations (const std::vector<float128_t> &conc);
        std::vector<float128_t> getY0 () const;
        const std::vector<std::function<float128_t (float128_t)>> &getGFunc () const;
        const std::vector<PhiFunction> &getPhiFunc () const;


        void initFromFile (const std::string &filename);

        void printInfo (std::ostream &out) const;

        //ChemicalTask getChemicalTask () const;

        uint64_t getCount () const;

        //std::vector<std::function<float128_t(const std::vector<float128_t> &)>> rightPartGen () const;
        void rightPartGen ();

        ChemicalSystem &operator= (const ChemicalSystem &reaction);
        ChemicalSystem &operator= (ChemicalSystem &&reaction);

        ChemicalReaction operator[] (uint64_t i) const;
        ChemicalReaction &operator[] (uint64_t i);

        static const float128_t R, P0;

    private:
        std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> table;
        std::vector<std::string> atoms;
        std::unordered_map<std::string, float128_t> atom_mass;
        std::vector<std::string> substances;
        std::unordered_map<std::string, float128_t> substance_mass;
        std::vector<std::string> additives;
        std::unordered_map<std::string, std::vector<float128_t>> additive_configs;
        std::vector<std::pair<float128_t, float128_t>> H0;
        std::vector<std::function<float128_t (float128_t)>> GFunc;
        std::vector<float128_t> concentrations;
        std::vector<PhiFunction> phi;
        std::vector<ChemicalReaction> reactions;
        float128_t T, P;
};

#endif