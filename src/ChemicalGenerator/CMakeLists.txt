add_library(ChemicalGenLIB "")

target_sources(ChemicalGenLIB
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/RightPartGen.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/RightPartGen.hpp
)

target_include_directories(ChemicalGenLIB
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})