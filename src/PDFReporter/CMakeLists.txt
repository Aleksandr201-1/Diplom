add_library(${PDF_REPORTER_LIB} "")

target_sources(${PDF_REPORTER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ReportGenerator.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/ReportGenerator.hpp
)

target_link_libraries(${PDF_REPORTER_LIB} 
    PRIVATE
        ${KOSHI_SOLVER_LIB}
        ${CHEMICAL_SOLVER_LIB}
        ${BUTCHER_LIB}
        ${NUMERIC_METHODS_LIB}
)

target_include_directories(${PDF_REPORTER_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})