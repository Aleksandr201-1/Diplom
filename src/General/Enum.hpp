#ifndef ENUM_HPP
#define ENUM_HPP

#include <string>
#include <map>
#include <algorithm>

/**
 * Преобразует элемент перечисления в строковое представление.
 *
 * @param enumEl Элемент перечисления, который нужно преобразовать в строку.
 * @param map Хеш таблица, которая связывает каждое значение перечисления с его соответствующим строковым представлением.
 * 
 * @return Строковое представление заданного элемента перечисления или "Error", если элемент не найден в карте.
**/
template <class T>
std::string enumToString (T enumEl, const std::map<T, std::string> &map) {
    auto it = map.find(enumEl);
    if (it != map.end()) {
        return it->second;
    }
    return "Error";
}

/**
 * Преобразует строковое представление в элемент перечисления.
 *
 * @param str Строковое представление элемента перечисления.
 * @param map Хеш таблица, которая связывает каждое значение перечисления с его соответствующим строковым представлением.
 * 
 * @return Элемент перечисления, соответствующий заданной строке, или значение типа T по умолчанию, если строка не найдена в карте.
**/
template <class T>
T stringToEnum (const std::string &str, const std::map<T, std::string> &map) {
    // лямбда-функция для проверки, имеет ли запись в карте заданное строковое значение
    auto check = [str] (const auto &pair) -> bool {
        return pair.second == str;
    };
    // находим первую запись в карте, значение которой соответствует заданной строке
    auto result = std::find_if(map.begin(), map.end(), check);
    if (result != map.end()) {
        return result->first;
    }
    // если соответствующая запись не найдена, возвращаем значение типа T по умолчанию
    return T::ERROR;
}

#endif