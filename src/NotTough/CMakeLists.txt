add_library(NotToughLIB "")

target_sources(NotToughLIB
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/RungeKutta.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/RungeKutta.hpp
)

target_link_libraries(NotToughLIB 
    PRIVATE
        GeneralLIB
)

target_include_directories(NotToughLIB
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})