#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <cmath>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <functional>
#include <iostream>
#include "FuncMaker.hpp"
#include "Matrix.hpp"

//структура, содержащая информацию по краевой задаче
struct Task {
    std::vector<FunctionalTree> trees; //массив функций, содержащий уравнение под индексом 0,
                                       //условие на левой границе под индексом 1, условие на правой границе под индексом 2
                                       //условие в начальный момент времени 
    std::vector<std::vector<double>> ux; //ux содержит пары коэффициентов ux и u для левого и правого условия
    std::vector<double> coeff, X; //coeff содержит коэффициенты уравнения при uxx, ux и u соответственно
                                  //X содержит левый и правый конец для оси x
    //uint64_t leftKind, rightKind;
};

const uint64_t ITERATION_CAP = 20;
const uint64_t PRECISION = 3;

void printVector(const std::vector<double> &vec);

//поиск машинного эпсилона
double findEpsillon ();

//сравнение 2-х double на равенство с заданной точностью
bool isEqual(double x, double y);

double derivative (const std::function<double(double)> &f, double x, uint64_t degree = 1);

//конвертация числа с плавающей точкой в строку до заданного знака
std::string toString (double val, uint64_t precision);

//функция приведения строки к форме, пригодной для парсера (удаление скобок по типу u(0) => u)
double stringFix (std::string &str);

//функция опработки ввода для получения информации по задаче
Task getTaskInfo(const std::vector<std::string> &system);

std::string readLine ();

#endif