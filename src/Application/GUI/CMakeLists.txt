add_executable(ODUSolverGUI
    GUITest.cpp
)

target_link_libraries(ODUSolverGUI
    PRIVATE
        ${PDF_REPORTER_LIB}
        ${KOSHI_SOLVER_LIB}
        ${CHEMICAL_SOLVER_LIB}
        ${GUI_LIB}
)

target_include_directories(ODUSolverGUI
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
        ${APPLICATION_PATH}
)

add_definitions(${COMPILE_FLAGS})