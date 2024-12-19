/**
    * @file ButcherTable.hpp
    * 
    * Определяет функцию и перечисление для работы с таблицами Бутчера.
**/

#ifndef BUTCHER_TABLE_HPP
#define BUTCHER_TABLE_HPP

#include <algorithm>
#include <map>
#include <Math/Matrix.hpp>
#include <General/Enum.hpp>

/**
    * @enum SolveMethod
    * 
    * @brief Перечисление методов решения ОДУ.
    *
    * @details
    * Определяет три типа узлов в дереве функции: операцию, значение и переменную.
**/
enum class SolveMethod {
    RUNGE_KUTTA,
    FALBERG,
    CHESKINO,
    MERSON,
    RADO,
    GAUSS,
    LOBATTO,
    L_STABLE_DIAGONAL,
    DORMAN_PRINCE,
    ERROR
};

/**
    * Функция, которая преобразует тип SolveMethod в строку.
    * 
    * @param method метод решения ОДУ типа SolveMethod
    * 
    * @return строковое представление метода
**/
std::string solveMethodToString (SolveMethod method);

/**
    * Функция, которая преобразует строку в тип SolveMethod.
    * 
    * @param str строковое представление метода решения ОДУ
    * 
    * @return соответствующий метод типа SolveMethod
**/
SolveMethod stringToSolveMethod (const std::string &str);

/**
    * Функция, которая создаёт таблицу Бутчера для заданного метода, порядка и способа
    * и возвращает её в виде объекта типа Matrix<double>.
    * 
    * @param method метод решения ОДУ типа SolveMethod
    * @param order порядок метода
    * @param way способ
    * 
    * @return таблица Бутчера типа Matrix<double>
**/
Matrix<double> createButcherTable (SolveMethod method, uint64_t order, uint64_t way);

#endif