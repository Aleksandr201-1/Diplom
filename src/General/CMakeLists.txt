add_library(${GENERAL_LIB} "")

target_sources(${GENERAL_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/FuncMaker.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ButcherTable.cpp
        ${CMAKE_CURRENT_LIST_DIR}/LSM.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Differentiation.cpp
        ${CMAKE_CURRENT_LIST_DIR}/General.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/FuncMaker.hpp
        ${CMAKE_CURRENT_LIST_DIR}/ButcherTable.hpp
        ${CMAKE_CURRENT_LIST_DIR}/LSM.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Differentiation.hpp
        ${CMAKE_CURRENT_LIST_DIR}/General.hpp
)

target_include_directories(${GENERAL_LIB}
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})