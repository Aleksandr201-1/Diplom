cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

#platform choice
set(WINDOWS 1)
set(LINUX 2)
set(DEBUG 1)

set(PLATFORM ${LINUX})

if(${PLATFORM} EQUAL ${WINDOWS})
  message(STATUS "Working for Window ver.")
  set(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++)
else()
  message(STATUS "Working for Linux ver.")
  set(CMAKE_CXX_COMPILER g++)
endif()

# libraries names
set(MODULES_INCLUDE_PATH ${CMAKE_BINARY_DIR}/src)
set(GENERAL_LIB GeneralLIB)
set(ODE_GENERAL_LIB ODEGeneralLIB)
set(KOSHI_SOLVER_LIB KoshiSolverLIB)
set(CHEMICAL_SOLVER_LIB ChemicalSolverLIB)
set(PDF_REPORTER_LIB PDFReporterLIB)
set(INPUT_LIB InputLIB)
#set(CHEMICAL_GENERATOR_LIB ChemicalGeneratorLIB)

set(COMPILE_WARNINGS -Wextra -Wall -Wno-sign-compare -Wno-unused-result)

#debugging
if(${DEBUG} EQUAL 1)
    message(STATUS "Its debugging time!")
    #set(COMPILE_FLAGS -Wextra -Wall -Wno-sign-compare -Wno-unused-result -O3)
    set(COMPILE_FLAGS -O3 ${COMPILE_WARNINGS})
else()
    set(COMPILE_FLAGS -Werror ${COMPILE_WARNINGS} -O3)
endif()

project(Diplom LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# defines targets and sources
add_subdirectory(src)
add_subdirectory(test)

# contains an "external" library we will link to
add_subdirectory(external)

# enable testing and define tests
#enable_testing()
#add_subdirectory(tests)