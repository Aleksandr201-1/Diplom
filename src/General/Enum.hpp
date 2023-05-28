#ifndef ENUM_HPP
#define ENUM_HPP

#include <string>
#include <map>

// template <class T>
// std::string enumToString (T enumEl, const std::map<T, std::string> &map) {
//     auto check = [enumEl] (const auto &pair) -> bool {
//         return pair.second == enumEl;
//     };
//     auto result = std::find_if(map.begin(), map.end(), check);
//     if (result != map.end()) {
//         return result->first;
//     }
//     return "NotAMethod";
// }

// template <class T>
// T stringToEnum (const std::string &str, const std::map<T, std::string> &map) {
//     auto it = map.find(str);
//     if (it != map.end()) {
//         return it->second;
//     }
//     return T::ERROR;
// }

template <class T>
std::string enumToString (T enumEl, const std::map<T, std::string> &map) {
    auto it = map.find(enumEl);
    if (it != map.end()) {
        return it->second;
    }
    return "Error";
}

template <class T>
T stringToEnum (const std::string &str, const std::map<T, std::string> &map) {
    auto check = [str] (const auto &pair) -> bool {
        return pair.second == str;
    };
    auto result = std::find_if(map.begin(), map.end(), check);
    if (result != map.end()) {
        return result->first;
    }
    return T::ERROR;
}

#endif