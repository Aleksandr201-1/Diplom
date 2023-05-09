add_library(${ODE_GENERAL_LIB} "")

add_subdirectory(Koshi)
add_subdirectory(Chemical)

target_sources(${ODE_GENERAL_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Task.cpp
        ${CMAKE_CURRENT_LIST_DIR}/IterationAlgo.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ToughDet.cpp
        #${CMAKE_CURRENT_LIST_DIR}/Koshi/KoshiTask.cpp
        #${CMAKE_CURRENT_LIST_DIR}/Koshi/KoshiSolver.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Task.hpp
        ${CMAKE_CURRENT_LIST_DIR}/IterationAlgo.hpp
        ${CMAKE_CURRENT_LIST_DIR}/ToughDet.hpp
        #${CMAKE_CURRENT_LIST_DIR}/Koshi/KoshiTask.hpp
        #${CMAKE_CURRENT_LIST_DIR}/Koshi/KoshiSolver.hpp
)

target_link_libraries(${ODE_GENERAL_LIB} 
    PRIVATE
        ${GENERAL_LIB}
)

target_include_directories(${ODE_GENERAL_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})