cmake_minimum_required(VERSION 3.8)

project(KrancNicolas VERSION 1.0 LANGUAGES CXX)

set(MAIN main)
set(MAIN_SOURCE General.hpp General.cpp FuncMaker.hpp FuncMaker.cpp Matrix.hpp RUN.hpp KrancNicolas.hpp KrancNicolas.cpp main.cpp)

add_executable(${MAIN} ${MAIN_SOURCE})

add_definitions(-std=c++11 -Wextra -Wall -Wno-sign-compare -Wno-unused-result -O3)