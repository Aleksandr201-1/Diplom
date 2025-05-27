add_library(${NUMERIC_METHODS_LIB} "")

target_sources(${NUMERIC_METHODS_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Integral.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Interpolation.cpp
        ${CMAKE_CURRENT_LIST_DIR}/CubeSpline.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Differentiation.cpp
        ${CMAKE_CURRENT_LIST_DIR}/LSM.tpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Integral.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Interpolation.hpp
        ${CMAKE_CURRENT_LIST_DIR}/CubeSpline.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Differentiation.hpp
        ${CMAKE_CURRENT_LIST_DIR}/LSM.hpp
        ${CMAKE_CURRENT_LIST_DIR}/SI.hpp
        ${CMAKE_CURRENT_LIST_DIR}/RUN.hpp
        ${CMAKE_CURRENT_LIST_DIR}/QR.hpp
        ${CMAKE_CURRENT_LIST_DIR}/LU.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Gauss.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Alg1Solve.hpp
)

target_include_directories(${NUMERIC_METHODS_LIB}
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

target_link_libraries(${NUMERIC_METHODS_LIB} 
    PRIVATE
        ${MATH_LIB}
        ${GENERAL_LIB}
)

add_definitions(${COMPILE_FLAGS})