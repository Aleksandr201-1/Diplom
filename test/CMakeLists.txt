set(MODULE_TEST_PATH ${CMAKE_BINARY_DIR}/test/ModuleTest)

function(build_module_tester TARGET_NAME SOURCES LIBRARIES)
    add_executable(${TARGET_NAME} ${MODULE_TEST_PATH}/src/${SOURCES})
    target_link_libraries(${TARGET_NAME}
        PRIVATE
            ${LIBRARIES}
            #${GENERAL_LIB}
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



	# add_library(${MODULE_NAME} OBJECT ${MODULE_NAME}.hpp ${MODULE_NAME}.cpp)
	# add_executable(${MODULE_NAME}Exmpl ${MODULE_NAME}Exmpl.cpp ${MODULE_NAME}.hpp ${MODULE_NAME}.cpp)
	# target_link_libraries(${MODULE_NAME}Exmpl ${SFML_LIBS})
	# if (NOT ${ARGN} EQUAL "")
	# 	#add_library(${MODULE_NAME} OBJECT ${MODULE_NAME}.hpp ${MODULE_NAME}.cpp)
	# 	target_link_libraries(${MODULE_NAME}Exmpl ${ARGN})
	# endif()
endfunction()

#FuncMakerTest
# add_executable(FuncMakerTest ${MODULE_TEST_PATH}/src/FuncMakerTest.cpp)
# target_link_libraries(FuncMakerTest
#     PRIVATE
#         ${GENERAL_LIB}
# )
# target_include_directories(FuncMakerTest 
#     PUBLIC 
#         ${MODULES_INCLUDE_PATH}
# )
build_module_tester(FuncMakerTest FuncMakerTest.cpp ${GENERAL_LIB})
build_module_tester(ChemicalGeneratorTest ChemicalGeneratorTest.cpp "${GENERAL_LIB};${CHEMICAL_GENERATOR_LIB}")

#ChemicalGeneratorTest
# add_executable(ChemicalGeneratorTest ${MODULE_TEST_PATH}/src/ChemicalGeneratorTest.cpp)
# target_link_libraries(ChemicalGeneratorTest
#     PRIVATE
#         ${GENERAL_LIB}
#         ${CHEMICAL_GENERATOR_LIB}
# )
# target_include_directories(ChemicalGeneratorTest
#     PUBLIC
#         ${MODULES_INCLUDE_PATH}
# )

#include(GNUInstallDirs)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
#   ${MODULE_TEST_PATH}/${CMAKE_INSTALL_LIBDIR})
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
#   ${MODULE_TEST_PATH}/${CMAKE_INSTALL_LIBDIR})
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
#   ${MODULE_TEST_PATH}/${CMAKE_INSTALL_BINDIR})

#set_target_properties(FuncMakerTest
#    PROPERTIES
#        ARCHIVE_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/lib
#        LIBRARY_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/lib
#        RUNTIME_OUTPUT_DIRECTORY ${MODULE_TEST_PATH}/bin
#)

add_definitions(${COMPILE_FLAGS})