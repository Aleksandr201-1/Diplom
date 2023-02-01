cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(WINDOWS 1)
set(LINUX 2)

set(PLATFORM ${LINUX})

if(${PLATFORM} EQUAL ${WINDOWS})
  message(STATUS "Working for Window ver.")
  set(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++)
else()
  message(STATUS "Working for Linux ver.")
  set(CMAKE_CXX_COMPILER g++)
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

# contains an "external" library we will link to
add_subdirectory(external)

# enable testing and define tests
#enable_testing()
#add_subdirectory(tests)