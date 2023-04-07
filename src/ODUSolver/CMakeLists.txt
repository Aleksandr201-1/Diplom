add_library(${ODU_SOLVER_LIB} "")

target_sources(${ODU_SOLVER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/RungeKutta.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/RungeKutta.hpp
)

target_link_libraries(${ODU_SOLVER_LIB} 
    PRIVATE
        ${GENERAL_LIB}
)

target_include_directories(${ODU_SOLVER_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})