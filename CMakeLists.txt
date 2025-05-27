cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

#platform choice
set(WINDOWS 1)
set(LINUX 2)
set(DEBUG 0)

#set(PLATFORM ${WINDOWS})

set(SFML_WINDOWS_LIBS sfml-graphics-2 sfml-window-2 sfml-audio-2 sfml-system-2 opengl32)
set(SFML_LINUX_LIBS sfml-graphics sfml-window sfml-audio sfml-system GL)
set(SFML_VERSION 2.5.1)

if(${PLATFORM} EQUAL ${WINDOWS})
    message(STATUS "Working for Windows ver.")
    set(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++)
    set(SFML_LIBS ${SFML_WINDOWS_LIBS})
    #set(SFML_PATH_TO_LIBS -L/${CMAKE_CURRENT_LIST_DIR}/external/Windows)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/../games/SFML-${SFML_VERSION}/include)
    link_directories(${CMAKE_CURRENT_LIST_DIR}/../games/SFML-${SFML_VERSION}/bin)
    #link_directories(${CMAKE_CURRENT_LIST_DIR}/external/Windows)
else()
    message(STATUS "Working for Linux ver.")
    set(CMAKE_CXX_COMPILER g++)
    set(SFML_LIBS ${SFML_LINUX_LIBS})
    set(SFML_PATH_TO_LIBS "")
endif()

# libraries names
set(MODULES_INCLUDE_PATH ${CMAKE_CURRENT_LIST_DIR}/src)
set(GENERAL_LIB GeneralLIB)
set(MATH_LIB MathLIB)
set(BUTCHER_LIB ButcherLIB)
set(NUMERIC_METHODS_LIB NMLIB)
set(ODE_GENERAL_LIB ODEGeneralLIB)
set(KOSHI_SOLVER_LIB KoshiSolverLIB)
set(CRANK_NICOLSON_SOLVER_LIB CrankNicolsonLIB)
set(CRANK_NICOLSON_SOLVER_LIB_2 CrankNicolsonLIB2)
set(CHEMICAL_SOLVER_LIB ChemicalSolverLIB)
set(PDF_REPORTER_LIB PDFReporterLIB)
set(INPUT_LIB InputLIB)
set(SF_EXTENSIONS_LIB SFExtensionsLIB)
set(GUI_LIB GUILIB)
#set(CHEMICAL_GENERATOR_LIB ChemicalGeneratorLIB)

set(COMPILE_WARNINGS -Wextra -Wall -Wno-sign-compare -Wno-unused-result)

#debugging
if(${DEBUG} EQUAL 1)
    message(STATUS "Its debugging time!")
    #set(COMPILE_FLAGS -Wextra -Wall -Wno-sign-compare -Wno-unused-result -O3)
    set(COMPILE_FLAGS ${COMPILE_WARNINGS})
else()
    #-Werror
    set(COMPILE_FLAGS ${COMPILE_WARNINGS} -O3)
endif()

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

project(Diplom LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
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
#add_subdirectory(test)

# contains an "external" library we will link to
add_subdirectory(external)

# enable testing and define tests
#enable_testing()
#add_subdirectory(tests)