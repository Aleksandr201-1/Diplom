set(MODULE_TEST_PATH ${CMAKE_BINARY_DIR}/test/ModuleTest)

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
    set_target_properties(${TARGET_NAME}
        PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/lib
            LIBRARY_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/lib
            RUNTIME_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/bin
    )
endfunction()

build_module_tester(FuncMakerTest FuncMakerTest.cpp ${GENERAL_LIB})
build_module_tester(Calculator FuncMakerCalculator.cpp ${GENERAL_LIB})
build_module_tester(ChemicalGeneratorTest ChemicalGeneratorTest.cpp ${CHEMICAL_SOLVER_LIB})
build_module_tester(DifferentiationTest DifferentiationTest.cpp ${GENERAL_LIB})

add_definitions(${COMPILE_FLAGS})