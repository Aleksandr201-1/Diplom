add_library(GeneralLIB "")

target_sources(GeneralLIB
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/FuncMaker.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ButcherTable.cpp
        ${CMAKE_CURRENT_LIST_DIR}/LSM.cpp
        ${CMAKE_CURRENT_LIST_DIR}/General.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/FuncMaker.hpp
        ${CMAKE_CURRENT_LIST_DIR}/ButcherTable.hpp
        ${CMAKE_CURRENT_LIST_DIR}/LSM.hpp
        ${CMAKE_CURRENT_LIST_DIR}/General.hpp
)

target_include_directories(GeneralLIB
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
)