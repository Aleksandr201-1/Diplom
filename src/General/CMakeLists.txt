add_library(${GENERAL_LIB} "")

target_sources(${GENERAL_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/FloatToString.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/FloatToString.hpp
        ${CMAKE_CURRENT_LIST_DIR}/CrashReporter.hpp
        ${CMAKE_CURRENT_LIST_DIR}/ErrorCodes.hpp
        ${CMAKE_CURRENT_LIST_DIR}/GlobalVars.hpp
        ${CMAKE_CURRENT_LIST_DIR}/ResourceHolder.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Underlying.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Enum.hpp
)

target_link_libraries(${GENERAL_LIB} 
    PRIVATE
        ${C_LIBS}
)

target_include_directories(${GENERAL_LIB}
    PRIVATE
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})