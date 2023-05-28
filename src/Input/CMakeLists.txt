add_library(${INPUT_LIB} "")

target_sources(${INPUT_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/InputHandle.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/InputHandle.hpp
)

target_link_libraries(${INPUT_LIB} 
    PRIVATE
        ${PDF_REPORTER_LIB}
)

target_include_directories(${INPUT_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})