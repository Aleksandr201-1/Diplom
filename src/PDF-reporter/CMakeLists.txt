add_library(ReportLIB "")

target_sources(ReportLIB
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ReportGenerator.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/ReportGenerator.hpp
)

target_include_directories(ReportLIB
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
)