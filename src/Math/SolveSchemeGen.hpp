#ifndef SOLVE_SCHEME_GEN_HPP
#define SOLVE_SCHEME_GEN_HPP

#include <Math/FuncMaker.hpp>

class SchemeGen {
    private:
        struct Coordinate {
            std::string name, step;
        };
        enum class NodeType {
            OPERATION,
            //DIFFERENTIATION,
            VALUE,
            VARIABLE
        };
        struct SGNode {
            public:
                SGNode (NodeType type);
                virtual ~SGNode ();
                friend SchemeGen;
            protected:
                NodeType type;
                uint64_t priority;
                std::vector<std::unique_ptr<SGNode>> leafs;
        };
    public:
        SchemeGen ();
        ~SchemeGen ();

        void addEquation (const std::string &equation);
        std::string translate () const;
        void setDerivativeRule ();
        void addCoord (const std::string &name, const std::string &step);

    private:
        std::vector<std::string> m_params;
        std::vector<Coordinate> m_coords;
        SGNode m_root;
};
/*
d(1/(rho*U*pow(y,nu)), x) = d(V/U, xi)
*/

#endif