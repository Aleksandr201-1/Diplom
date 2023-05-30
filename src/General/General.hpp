/**
    * @file General.hpp
    *
    * @brief Общие функции и типы данных для вычислительных задач.
**/
#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <cmath>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <functional>
#include <iostream>

using float128_t = long double;

const uint64_t ITERATION_CAP = 200;
const uint64_t PRECISION = 2;

/**
    * @brief Выводит вектор со значениями в стандартный поток вывода.
    *
    * @param vec Вектор значений для вывода.
**/
void printVector(const std::vector<float128_t> &vec);

/**
    * @brief Проверяет, равны ли два числа с учетом погрешности.
    * 
    * @param x Первое число для сравнения.
    * @param y Второе число для сравнения.
    * 
    * @return true, если числа равны.
    * @return false, если числа не равны.
**/
bool isEqual(float128_t x, float128_t y);

/**
    * @brief Преобразует значение числа с плавающей точкой в строку с указанной точностью.
    * 
    * @param val Значение числа с плавающей точкой.
    * @param precision Точность, до которой нужно округлить число.
    * 
    * @return std::string Строка, содержащая преобразованное значение числа.
**/
std::string toString (float128_t val, uint64_t precision);

#endif