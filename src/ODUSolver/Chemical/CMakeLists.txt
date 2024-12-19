add_library(${CHEMICAL_SOLVER_LIB} "")

target_sources(${CHEMICAL_SOLVER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ChemicalTask.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ChemicalSolver.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/ChemicalTask.hpp
        ${CMAKE_CURRENT_LIST_DIR}/ChemicalSolver.hpp
)

target_link_libraries(${CHEMICAL_SOLVER_LIB} 
    PRIVATE
        ${ODE_GENERAL_LIB}
        ${NUMERIC_METHODS_LIB}
        ${BUTCHER_LIB}
)

target_include_directories(${CHEMICAL_SOLVER_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})