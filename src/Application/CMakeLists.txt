set(APPLICATION_PATH ${CMAKE_CURRENT_LIST_DIR})

function(build_app TARGET_NAME SOURCES EXTRA_LIBS)
    add_executable(${TARGET_NAME} ${SOURCES})
    target_link_libraries(${TARGET_NAME}
        PRIVATE
            ${PDF_REPORTER_LIB}
            ${KOSHI_SOLVER_LIB}
            ${CRANK_NICOLSON_SOLVER_LIB}
            ${CHEMICAL_SOLVER_LIB}
            ${EXTRA_LIBS}
            ${NUMERIC_METHODS_LIB}
            ${BUTCHER_LIB}
    )
    target_include_directories(${TARGET_NAME} 
        PUBLIC
            ${MODULES_INCLUDE_PATH}
            ${CMAKE_CURRENT_LIST_DIR}
            ${APPLICATION_PATH}
    )
    # set_target_properties(${TARGET_NAME}
    #     PROPERTIES
    #         ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TEST_DIR}/${CMAKE_INSTALL_LIBDIR}
    #         LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TEST_DIR}/${CMAKE_INSTALL_LIBDIR}
    #         RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TEST_DIR}/${CMAKE_INSTALL_BINDIR}
    # )
    add_definitions(${COMPILE_FLAGS})
endfunction()

add_subdirectory(Console)
add_subdirectory(GUI)