add_library(${KOSHI_SOLVER_LIB} "")

target_sources(${KOSHI_SOLVER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/KoshiTask.cpp
        ${CMAKE_CURRENT_LIST_DIR}/KoshiSolver.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/KoshiTask.hpp
        ${CMAKE_CURRENT_LIST_DIR}/KoshiSolver.hpp
)

target_link_libraries(${KOSHI_SOLVER_LIB} 
    PRIVATE
        ${ODE_GENERAL_LIB}
)

target_include_directories(${KOSHI_SOLVER_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})