#ifndef A_HPP
#define A_HPP

#include "General.hpp"
#include "RUN.hpp"
#include "SI.hpp"
#include "LU.hpp"

//методы: явный, неявный и Кранка-Николаса
enum class Method {
    EXPLICIT,
    NOT_EXPLICIT,
    KRANK_NICOLAS
};

enum class SIDE {
    LEFT,
    RIGHT
};

//аппроксимация
enum class ApproxLevel {
    POINT2_ORDER1, //двухточечная первого порядка
    POINT2_ORDER2, //двухточечная второго порядка
    POINT3_ORDER2, //трёхточечная второго порядка
    NONE
};

/**
    * \brief Функция решения краевой задачи. Принимает следующие аргументы
    *
    *   
    * \param task Структура, содержащая задачу.
    * \param timeLimit Верхний предел времени.
    * \param xh Размер шага для X.
    * \param th Размер шага для T.
    * \param method Метод решения краевой задачи.
    * \param approx Метод аппроксимации производной.
    * 
    * \return Таблица функции U(x, t)
    *
**/
std::vector<std::vector<double>> SolveIBVP (const Task &task, double timeLimit, double xh, double th, Method method, ApproxLevel approx);

#endif