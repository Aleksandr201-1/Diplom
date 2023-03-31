set(TEST_SRC ModuleTest/src)

add_executable(FuncMakerTest ${TEST_SRC}/FuncMakerTest.cpp)

#add_subdirectory(General)
#add_subdirectory(PDF-reporter)
#add_subdirectory(Tough)
#add_subdirectory(NotTough)
#add_subdirectory(ChemicalGenerator)

target_link_libraries(FuncMakerTest
    PRIVATE
        GeneralLIB
)

set_target_properties(FuncMakerTest
    PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/ModuleTest/lib
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/ModuleTest/lib
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/ModuleTest/bin
)

add_definitions(${COMPILE_FLAGS})