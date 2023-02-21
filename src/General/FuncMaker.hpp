#ifndef FUNCTIONAL_TREE_HPP
#define FUNCTIONAL_TREE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>

// TODO:
// V исправить унарный оператор
// V добавить числовые константы
// - улучшить нахождение коэффициентов
// V улучшить печать функции
// - добавить конвертацию в массив
// V изменить формат добавления новых функций и их аналогов
// - добавить поддержку функций от 2х и более переменных
// - добавить вложенные функции

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
    SIGN,   // sign
    NOT_AN_OPERATION
};

enum class Style {
    DEFAULT,
    GNUPLOT,
    LATEX
};

struct OperationStruct {
    std::vector<std::string> op_str;
    std::function<double (double, double)> func;
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
        //проверка имён переменных на корректность
        void inputCheck (const std::vector<std::string> &vars) const; 
        //чтение операции
        std::string readOperation (const std::string &func, uint64_t &i) const;
        //чтение слова
        std::string readWord (const std::string &func, uint64_t &i) const;
        //чтение числа
        double readNumber (const std::string &func, uint64_t &i) const;
        //чтение выражения в скобках
        std::string readInbrace (const std::string &func, uint64_t &i) const;
        //конвертация строки в операцию
        Operation getOperation (const std::string &str) const;
        //конвертация строки в числовую константу
        double getConstant (const std::string &str) const;
        //получение приоритета операции (0 - выполяться сразу, 4 - выполняться в конце)
        uint64_t getPriority (Operation op) const;
        // использование операции на 2х переменных
        double useOperation (Operation op, double x, double y) const;
        //вычисление выражение в поддереве node
        double getVal (const NodePtr &node, const std::vector<double> &X) const;
        //добавление узла в дерево
        void addToTree (NodePtr &tree, NodePtr &node);
        //построение дерева по строке
        NodePtr buildTree (const std::string &func);
        //вспомогательные функции для вывода дерева
        void printTree (const NodePtr &node, std::ostream &out) const;
        void printFunc (const NodePtr &node, std::ostream &out) const;
        //void printNode (const NodePtr &node) const;
        void toStringDefault (const NodePtr &node, std::string &str) const;
        //void toStringGNUPlot (const NodePtr &node, std::string &str) const;
        void toStringLatex (const NodePtr &node, std::string &str) const;
        //копирование дерева
        NodePtr copyTree (const NodePtr &node) const;
        //конструкторы копирования поддерева
        FunctionalTree (const NodePtr &node);
        FunctionalTree (NodePtr &&tree);
    public:
        //конструкторы
        FunctionalTree ();
        FunctionalTree (const std::string &func);
        FunctionalTree (const std::string &func, const std::string &var);
        FunctionalTree (const std::string &func, const std::vector<std::string> &vars);
        FunctionalTree (const FunctionalTree &tree);
        FunctionalTree (FunctionalTree &&tree);
        ~FunctionalTree ();
        //сброс дерева
        void reset (const std::string &func, const std::vector<std::string> &vars);
        //вызов функции
        double func (double x) const;
        double func (const std::vector<double> &X) const;
        //посчитать выражение без переменных
        double calculate () const;
        //список переменных
        std::vector<std::string> getVarList () const;
        //коэффициент при переменной
        FunctionalTree getCoeff (uint64_t idx) const;
        FunctionalTree getCoeff (const std::string &param) const;
        //FunctionalTree getDiv () const;
        //печать содержимого
        void printTree () const;
        void printFunc () const;
        std::string toString (Style style) const;
        //void simplify ();
        FunctionalTree &operator= (const FunctionalTree &tree);
        FunctionalTree &operator= (FunctionalTree &&tree);

        //оператор функции
        double operator() (double x) const;
        double operator() (const std::vector<double> &X) const;

        //вывод
        friend std::ostream &operator<< (std::ostream &output, const FunctionalTree &tree);

        //чтение и запись из файла
        friend std::ifstream &operator>> (std::ifstream &file, FunctionalTree &tree);
        friend std::ofstream &operator<< (std::ofstream &file, const FunctionalTree &tree);
    private:
        static const std::vector<OperationStruct> operations;
        static const std::vector<ConstantValue> const_val;
        static const uint64_t VARIABLE_LIMIT = 10;
        std::vector<std::string> vars;
        NodePtr root;
};

#endif