#cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

#project(ODUSolverGUI LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(QT_LIBS Qt5::Core Qt5::Widgets)

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# include(GNUInstallDirs)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
#   ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
#   ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
#   ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# QtCreator supports the following variables for Android, which are identical to qmake Android variables.
# Check http://doc.qt.io/qt-5/deployment-android.html for more information.
# They need to be set before the find_package(Qt5 ...) call.

#if(ANDROID)
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
#    if (ANDROID_ABI STREQUAL "armeabi-v7a")
#        set(ANDROID_EXTRA_LIBS
#            ${CMAKE_CURRENT_SOURCE_DIR}/path/to/libcrypto.so
#            ${CMAKE_CURRENT_SOURCE_DIR}/path/to/libssl.so)
#    endif()
#endif()

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
# find_package(Qt5 COMPONENTS PrintSupport REQUIRED)
# add_library(QCustomPlot STATIC QCustomPlot.cpp)
# target_link_libraries(QCustomPlot PRIVATE Qt5::PrintSupport)
# target_include_directories(QCustomPlot INTERFACE .)

if(ANDROID)
    add_library(ODUSolverGUI SHARED
        main.cpp
        #QCustomPlot.cpp
        #QCustomPlot.hpp
        MainWindow.cpp
        MainWindow.hpp
        MainWindow.ui
    )
else()
    add_executable(ODUSolverGUI
        main.cpp
        #QCustomPlot.cpp
        #QCustomPlot.hpp
        MainWindow.cpp
        MainWindow.hpp
        MainWindow.ui
    )
endif()

target_link_libraries(ODUSolverGUI
    PRIVATE
        ${QT_LIBS}
        ${KOSHI_SOLVER_LIB}
        ${CHEMICAL_SOLVER_LIB}
        #QCustomPlot
        #qwt
)

target_include_directories(ODUSolverGUI
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
        ${APPLICATION_PATH}
)

add_definitions(${COMPILE_FLAGS})