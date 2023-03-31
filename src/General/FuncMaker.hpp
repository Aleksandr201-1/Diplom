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
#include <map>

// TODO:
// V исправить унарный оператор
// V добавить числовые константы
// - улучшить нахождение коэффициентов
// V улучшить печать функции
// - добавить конвертацию в массив
// V изменить формат добавления новых функций и их аналогов
// - добавить поддержку функций от 2х и более переменных
// - добавить вложенные функции
// - добавить вложенные переменные

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

// enum class Operation {
//     PLUS,   // +
//     MINUS,  // -
//     MUL,    // *
//     DIV,    // /
//     MOD,    // %
//     POW,    // ^, **
//     SQRT,   // sqrt
//     SIN,    // sin
//     COS,    // cos
//     TAN,    // tg
//     CTG,    // ctg
//     SINH,   // sinh
//     COSH,   // cosh
//     TANH,   // tanh
//     CTH,    // cth
//     ASIN,   // arcsin
//     ACOS,   // arccos
//     ATAN,   // arctg
//     ACOT,   // arcctg
//     LOG,    // log_10
//     LN,     // log_e, ln
//     EXP,    // exp
//     ABS,    // abs, ||
//     SIGN,   // sign
//     NOT_AN_OPERATION
// };

struct OperationStruct {
    std::vector<std::string> op_str;
    std::function<double (double, double)> func;
    uint64_t priority;

    OperationStruct (const std::vector<std::string> &op_str, const std::function<double (double, double)> &func, uint64_t priority);
    ~OperationStruct ();
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
        OperationNode (uint64_t idx);
        ~OperationNode ();
        friend FunctionalTree; 
    private:
        uint64_t idx;
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
        uint64_t getOperation (const std::string &str) const;
        //конвертация строки в числовую константу
        double getConstant (const std::string &str) const;
        //конвертация строки в индекс перпеменной
        uint64_t getVariable (const std::string &str) const;
        //получение приоритета операции (0 - выполяться сразу, 4 - выполняться в конце)
        uint64_t getPriority (uint64_t op) const;
        // использование операции на 2х переменных
        double useOperation (uint64_t op, double x, double y) const;
        //вычисление выражение в поддереве node
        double calcNode (const NodePtr &node, const std::vector<double> &X) const;
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
        enum class Style {
            DEFAULT,
            GNUPLOT,
            LATEX
        };
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
        //добавить значение
        void setValue (const std::string &name, double val);
        //получение значения по имени (если значения не существует, то возвращает nan)
        double getValue (const std::string &name) const;
        //вызов функции для 1 аргумента
        double func (double x) const;
        //вызов функции для произвольного числа аргументов
        double func (const std::vector<double> &X) const;
        //посчитать выражение без переменных (если в функции есть переменные, то их значение берётся равным нулю)
        double calculate () const;
        //список переменных
        std::vector<std::string> getVariableList () const;
        //список значений
        std::vector<std::string> getValueList () const;
        //коэффициент при переменной по индексу
        FunctionalTree getCoeff (uint64_t idx) const;
        //коэффициент при переменной по названию
        FunctionalTree getCoeff (const std::string &param) const;
        //FunctionalTree getDiv () const;
        //печать содержимого в виде дерева
        void printTree () const;
        //печать содержимого в виде функции
        void printFunc () const;
        std::string toString (Style style) const;
        //void simplify ();
        //оператор копирования
        FunctionalTree &operator= (const FunctionalTree &tree);
        //оператор перемещения
        FunctionalTree &operator= (FunctionalTree &&tree);

        //оператор функции от 1 переменной
        double operator() (double x) const;
        //оператор функции для произвольного числа переменных
        double operator() (const std::vector<double> &X) const;

        //вывод
        friend std::ostream &operator<< (std::ostream &output, const FunctionalTree &tree);

        //чтение и запись из файла
        friend std::ifstream &operator>> (std::ifstream &file, FunctionalTree &tree);
        friend std::ofstream &operator<< (std::ofstream &file, const FunctionalTree &tree);
    private:
        static const std::vector<OperationStruct> operations; //список операций
        static const std::map<std::string, double> const_vals; //список константных значений (pi, e и т.д.)
        //static const std::vector<ConstantValue> const_val;
        //static const uint64_t VARIABLE_LIMIT = 10;
        std::map<std::string, double> vals; //список пользовательских значений
        std::vector<std::string> vars; //список переменных
        //std::vector<int> aa;
        NodePtr root;
};

#endif