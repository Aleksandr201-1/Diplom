cmake_minimum_required(VERSION 3.8)

project(Diplom VERSION 0.0.1 LANGUAGES CXX)

#add_subdirectory()
set(MAIN main)

#set(COMPILE_KEYS -std=c++11 -Wextra -Wall -Werror -Wno-sign-compare -Wno-unused-result -O3)
set(COMPILE_KEYS -std=c++11 -Wextra -Wall -Wno-sign-compare -Wno-unused-result)

set(G General)
set(N NotTough)
set(T Tough)

target_include_directories(main PRIVATE ${G})
set(GENERAL_SOURCE ${G}/General.hpp ${G}/General.cpp ${G}/FuncMaker.hpp ${G}/FuncMaker.cpp)

target_include_directories(main PRIVATE ${N})
set(NOT_TOUGH_SOURCE ${N}/RungeKutta.hpp ${N}/RungeKutta.hpp)

#target_include_directories(main PRIVATE ${T})
#set(TOUGH_SOURCE ${T}/RungeKutta.hpp ${T}/RungeKutta.cpp)

set(MAIN_SOURCE ${GENERAL_SOURCE} ${NOT_TOUGH_SOURCE} NotTough.hpp main.cpp)

#set(BOOST_LIBS boost_iostreams boost_system boost_filesystem)

add_executable(${MAIN} ${MAIN_SOURCE})
#target_link_libraries(${MAIN} ${BOOST_LIBS})

add_definitions(${COMPILE_KEYS})