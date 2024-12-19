add_library(${BUTCHER_LIB} "")

target_sources(${BUTCHER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ButcherTable.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/ButcherTable.hpp
)

target_include_directories(${BUTCHER_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

target_link_libraries(${BUTCHER_LIB} 
    PRIVATE
        ${MATH_LIB}
        ${GENERAL_LIB}
)

add_definitions(${COMPILE_FLAGS})