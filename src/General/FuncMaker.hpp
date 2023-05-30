/**
    * @file FunctionalTree.hpp
    * 
    * Определяет класс `FunctionalTree` для хранения и вычисления аналитических функций.
**/

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
#include "General.hpp"

// TODO:
// V исправить унарный оператор
// V добавить числовые константы
// - улучшить нахождение коэффициентов
// V улучшить печать функции
// - добавить конвертацию в массив
// V изменить формат добавления новых функций и их аналогов
// - добавить поддержку функций от 2х и более переменных
// - добавить вложенные функции
// V добавить вложенные переменные

//PRIORITY
//val, ()       0
//^             1
//sin, cos, ... 2
//*, /, %       3
//+, -          4


/**
    * @enum NodeType
    * 
    * @brief Перечисление типов узлов дерева функции.
    *
    * @details Определяет три типа узлов в дереве функции: операцию, значение и переменную.
**/
enum class NodeType {
    OPERATION, /**Операция.*/
    VALUE, /**Значение.*/
    VARIABLE /**Переменная.*/
};

/**
    * @struct OperationStruct
    * 
    * @brief Структура для хранения информации об операции.
    *
    * @details Структура включает в себя список строковых представлений операции (например, "+"), функцию для выполнения операции
    * и приоритет операции.
 */
struct OperationStruct {
    std::vector<std::string> op_str; /**Список строковых представлений операции. */
    std::function<float128_t (float128_t, float128_t)> func; /**Функция для выполнения операции. */
    uint64_t priority; /**Приоритет операции. */

    /**
        * Создает новый объект структуры `OperationStruct`.
        * 
        * @param op_str Список строковых представлений операции.
        * @param func Функция для выполнения операции.
        * @param priority Приоритет операции.
     */
    OperationStruct (const std::vector<std::string> &op_str, const std::function<float128_t (float128_t, float128_t)> &func, uint64_t priority);
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
        ValueNode (float128_t val);
        ~ValueNode ();
        friend FunctionalTree;
    private:
        float128_t val;
};

class VariableNode : public FunctionalTreeNode {
    public:
        VariableNode (uint64_t idx);
        ~VariableNode ();
        friend FunctionalTree;
    private:
        uint64_t idx;
};

// class FunctionNode : public FunctionalTreeNode {
//     public:
//         FunctionNode (uint64_t idx);
//         ~FunctionNode ();
//         friend FunctionalTree;
//     private:
//         uint64_t idx;
// };

/**
    * @class FunctionalTree
    * 
    * @brief Класс для представления математических функций в виде бинарных деревьев.
    *
    * @details
    * Этот класс может использоваться для представления математических функций в виде бинарных деревьев.
    * Дерево строится из строкового представления функции, которое может содержать переменные,
    * константы и основные арифметические операции, такие как сложение, вычитание, умножение,
    * деление и возведение в степень. После построения дерева его можно вычислить для определенных значений переменных
    * или констант, а также вывести в различных форматах.
    * 
    * Использование:
    * @code{.cpp}
    * // список переменных
    * std::vector<std::string> args = {"x", "y"};
    * // инициализация объекта FunctionalTree
    * FunctionalTree func("x^2 + 2*x*y + y^2", args);
    * 
    * std::cout << "func(3, 4) = " << func({3.0, 4.0}) << '\n';
    * @endcode 
**/
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
        float128_t readNumber (const std::string &func, uint64_t &i) const;
        //чтение выражения в скобках
        std::string readInbrace (const std::string &func, uint64_t &i) const;
        //конвертация строки в операцию
        uint64_t getOperation (const std::string &str) const;
        //конвертация строки в числовую константу
        float128_t getConstant (const std::string &str) const;
        //конвертация строки в индекс перпеменной
        uint64_t getVariable (const std::string &str) const;
        //получение приоритета операции (0 - выполяться сразу, 4 - выполняться в конце)
        uint64_t getPriority (uint64_t op) const;
        // использование операции на 2х переменных
        float128_t useOperation (uint64_t op, float128_t x, float128_t y) const;
        //вычисление выражение в поддереве node
        float128_t calcNode (const NodePtr &node, const std::vector<float128_t> &X) const;
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
        void setValue (const std::string &name, float128_t val);
        //получение значения по имени (если значения не существует, то возвращает nan)
        float128_t getValue (const std::string &name) const;
        //вызов функции для 1 аргумента
        float128_t func (float128_t x) const;
        //вызов функции для произвольного числа аргументов
        float128_t func (const std::vector<float128_t> &X) const;
        //посчитать выражение без переменных (если в функции есть переменные, то их значение берётся равным нулю)
        float128_t calculate () const;
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
        float128_t operator() (float128_t x) const;
        //оператор функции для произвольного числа переменных
        float128_t operator() (const std::vector<float128_t> &X) const;

        //вывод
        friend std::ostream &operator<< (std::ostream &output, const FunctionalTree &tree);

        //чтение и запись из файла
        friend std::ifstream &operator>> (std::ifstream &file, FunctionalTree &tree);
        friend std::ofstream &operator<< (std::ofstream &file, const FunctionalTree &tree);
    private:
        static const std::vector<OperationStruct> operations; //список операций
        static const std::map<std::string, float128_t> const_vals; //список константных значений (pi, e и т.д.)
        //static const std::vector<ConstantValue> const_val;
        //static const uint64_t VARIABLE_LIMIT = 10;
        std::map<std::string, float128_t> vals; //список пользовательских значений
        std::vector<std::string> vars; //список переменных
        //std::vector<int> aa;
        NodePtr root;
};

#endif