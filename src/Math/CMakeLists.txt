add_library(${MATH_LIB} "")

target_sources(${MATH_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/BasicMath.cpp
        ${CMAKE_CURRENT_LIST_DIR}/FuncMaker.cpp
        ${CMAKE_CURRENT_LIST_DIR}/LongDouble.cpp
        ${CMAKE_CURRENT_LIST_DIR}/LongInt.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/BasicMath.hpp
        ${CMAKE_CURRENT_LIST_DIR}/FuncMaker.hpp
        ${CMAKE_CURRENT_LIST_DIR}/LongDouble.hpp
        ${CMAKE_CURRENT_LIST_DIR}/LongInt.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Matrix.hpp
)

target_include_directories(${MATH_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})