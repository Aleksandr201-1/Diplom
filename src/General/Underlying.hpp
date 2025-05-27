#ifndef UNDERLYING_HPP
#define UNDERLYING_HPP

#include <type_traits>

template <typename T>
constexpr auto to_underlying(T e) noexcept {
    return static_cast<std::underlying_type_t<T>>(e);
}

#endif