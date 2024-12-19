add_library(${CRANK_NICOLSON_SOLVER_LIB} "")

target_sources(${CRANK_NICOLSON_SOLVER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/CrankNicolson.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/CrankNicolson.hpp
)

target_link_libraries(${CRANK_NICOLSON_SOLVER_LIB} 
    PRIVATE
        ${ODE_GENERAL_LIB}
        ${GENERAL_LIB}
        ${MATH_LIB}
        ${NUMERIC_METHODS_LIB}
        ${CHEMICAL_SOLVER_LIB}
)

target_include_directories(${CRANK_NICOLSON_SOLVER_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})