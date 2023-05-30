add_executable(ODUSolverConsole main.cpp)

target_link_libraries(ODUSolverConsole
    PRIVATE
        ${KOSHI_SOLVER_LIB}
        ${CHEMICAL_SOLVER_LIB}
        ${PDF_REPORTER_LIB}
        ${INPUT_LIB}
)

target_include_directories(ODUSolverConsole
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
        ${APPLICATION_PATH}
)

add_definitions(${COMPILE_FLAGS})