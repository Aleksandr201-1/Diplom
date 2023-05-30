set(MODULE_TEST_PATH ${CMAKE_CURRENT_LIST_DIR}/ModuleTest)

#set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
#set(CMAKE_DISABLE_SOURCE_CHANGES ON)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

function(build_module_tester TARGET_NAME SOURCES LIBRARIES)
    add_executable(${TARGET_NAME} ${MODULE_TEST_PATH}/src/${SOURCES})
    target_link_libraries(${TARGET_NAME}
        PRIVATE
            ${LIBRARIES}
    )
    target_include_directories(${TARGET_NAME} 
        PUBLIC
            ${MODULES_INCLUDE_PATH}
    )
    # set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    #     ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    #     set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    #     ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    #     set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    #     ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # set_target_properties(${TARGET_NAME}
    #     PROPERTIES
    #         ARCHIVE_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/lib
    #         LIBRARY_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/lib
    #         RUNTIME_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/bin
    # )
    set_target_properties(${TARGET_NAME}
        PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/${CMAKE_INSTALL_LIBDIR}
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/${CMAKE_INSTALL_LIBDIR}
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/${CMAKE_INSTALL_BINDIR}
    )
endfunction()

build_module_tester(FuncMakerTest FuncMakerTest.cpp ${GENERAL_LIB})
build_module_tester(Calculator FuncMakerCalculator.cpp ${GENERAL_LIB})
build_module_tester(ChemicalGeneratorTest ChemicalGeneratorTest.cpp ${CHEMICAL_SOLVER_LIB})
build_module_tester(DifferentiationTest DifferentiationTest.cpp ${GENERAL_LIB})

add_definitions(${COMPILE_FLAGS})