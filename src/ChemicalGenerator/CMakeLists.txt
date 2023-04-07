add_library(${CHEMICAL_GENERATOR_LIB} "")

target_sources(${CHEMICAL_GENERATOR_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/RightPartGen.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/RightPartGen.hpp
)

target_link_libraries(${CHEMICAL_GENERATOR_LIB} 
    PRIVATE
        ${GENERAL_LIB}
)

target_include_directories(${CHEMICAL_GENERATOR_LIB}
    PRIVATE
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})