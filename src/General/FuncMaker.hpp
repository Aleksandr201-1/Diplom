#ifndef FUNCTIONAL_TREE_HPP
#define FUNCTIONAL_TREE_HPP

#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>

// TODO:
// - исправить унарный оператор
// - добавить числовые константы
// - улучшить нахождение коэффициентов
// - улучшить печать функции
// - добавить конвертацию в массив
// - изменить формат добавления новых функций и их аналогов

struct OperationStruct {
    std::vector<std::string> op_str;
    std::function<double (double, double)> func;
    //Operation op;
    uint64_t priority;

    OperationStruct (const std::vector<std::string> &op_str, const std::function<double (double, double)> &func, uint64_t priority);
    ~OperationStruct ();
};

struct ConstantValue {
    std::string val_name;
    double val;

    ConstantValue (const std::string &val_name, double val);
    ~ConstantValue ();
};

//PRIORITY
//val, ()       0
//^             1
//sin, cos, ... 2
//*, /, %       3
//+, -          4
enum class NodeType {
    OPERATION,
    VALUE,
    VARIABLE
};

enum class Operation {
    PLUS,   // +
    MINUS,  // -
    MUL,    // *
    DIV,    // /
    MOD,    // %
    POW,    // ^, **
    SQRT,   // sqrt
    SIN,    // sin
    COS,    // cos
    TAN,    // tg
    CTG,    // ctg
    SINH,   // sinh
    COSH,   // cosh
    TANH,   // tanh
    CTH,    // cth
    ASIN,   // arcsin
    ACOS,   // arccos
    ATAN,   // arctg
    ACOT,   // arcctg
    LOG,    // log_10
    LN,     // log_e, ln
    EXP,    // exp
    ABS,    // abs, ||
    NOT_AN_OPERATION
};

class FunctionalTree;

class FunctionalTreeNode {
    public:
        FunctionalTreeNode (NodeType type);
        virtual ~FunctionalTreeNode ();
        friend FunctionalTree;
    protected:
        NodeType type;
        uint64_t priority;
        std::unique_ptr<FunctionalTreeNode> left, right;
};

class OperationNode : public FunctionalTreeNode {
    public:
        OperationNode (Operation op);
        ~OperationNode ();
        friend FunctionalTree; 
    private:
        Operation op;
};

class ValueNode : public FunctionalTreeNode {
    public:
        ValueNode (double val);
        ~ValueNode ();
        friend FunctionalTree;
    private:
        double val;
};

class VariableNode : public FunctionalTreeNode {
    public:
        VariableNode (uint64_t idx);
        ~VariableNode ();
        friend FunctionalTree;
    private:
        uint64_t idx;
};

//класс для представления строки в функцию
class FunctionalTree {
    private:
        using NodePtr = std::unique_ptr<FunctionalTreeNode>;
    private:
        void inputCheck (const std::vector<std::string> &vars) const; 
        std::string readOperation (const std::string &func, uint64_t &i) const;
        std::string readWord (const std::string &func, uint64_t &i) const;
        double readNumber (const std::string &func, uint64_t &i) const;
        std::string readInbrace (const std::string &func, uint64_t &i) const;
        Operation getOperation (const std::string &str) const;
        double getConstant (const std::string &str) const;
        uint64_t getPriority (Operation op) const;
        double useOperation (Operation op, double x, double y) const;
        double getVal (const NodePtr &node, const std::vector<double> &X) const;
        void addToTree (NodePtr &root, NodePtr &toAdd);
        NodePtr buildTree (const std::string &func);
        void printTree (const NodePtr &node, std::ostream &out) const;
        void printFunc (const NodePtr &node, std::ostream &out) const;
        NodePtr copyTree (const NodePtr &node) const;
        FunctionalTree (const NodePtr &node);
        FunctionalTree (NodePtr &&tree);
    public:
        FunctionalTree ();
        FunctionalTree (const std::string &func);
        FunctionalTree (const std::string &func, const std::string &var);
        FunctionalTree (const std::string &func, const std::vector<std::string> &vars);
        FunctionalTree (const FunctionalTree &tree);
        FunctionalTree (FunctionalTree &&tree);
        ~FunctionalTree ();
        void reset (const std::string &func, const std::vector<std::string> &vars);
        double func (double x) const;
        double func (const std::vector<double> &X) const;
        double calculate () const;
        std::vector<std::string> getVarList () const;
        FunctionalTree getCoeff (uint64_t idx) const;
        FunctionalTree getCoeff (const std::string &param) const;
        FunctionalTree getDiv () const;
        void printTree () const;
        void printFunc () const;
        //void simplify ();
        FunctionalTree &operator= (const FunctionalTree &tree);
        FunctionalTree &operator= (FunctionalTree &&tree);
        double operator() (double x) const;
        double operator() (const std::vector<double> &X) const;

        //вывод
        friend std::ostream &operator<< (std::ostream &output, const FunctionalTree &tree);

        //чтение и запись из файла
        friend std::ifstream &operator>> (std::ifstream &file, FunctionalTree &tree);
        friend std::ofstream &operator<< (std::ofstream &file, const FunctionalTree &tree);
    private:
        //static const std::vector<std::string> operations;
        static const std::vector<OperationStruct> operations;
        static const std::vector<ConstantValue> const_val;
        static const uint64_t VARIABLE_LIMIT = 10;
        std::vector<std::string> vars;
        NodePtr root;
};

#endif