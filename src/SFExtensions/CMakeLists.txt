add_library(${SF_EXTENSIONS_LIB} "")

target_sources(${SF_EXTENSIONS_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/EllipseShape.cpp
        ${CMAKE_CURRENT_LIST_DIR}/SquircleShape.cpp
        ${CMAKE_CURRENT_LIST_DIR}/RhombShape.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Line.cpp
        ${CMAKE_CURRENT_LIST_DIR}/SpriteManager.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/EllipseShape.hpp
        ${CMAKE_CURRENT_LIST_DIR}/SquircleShape.hpp
        ${CMAKE_CURRENT_LIST_DIR}/RhombShape.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Line.hpp
        ${CMAKE_CURRENT_LIST_DIR}/SpriteManager.hpp
)

target_link_libraries(${SF_EXTENSIONS_LIB} 
    PRIVATE
        ${SFML_LIBS}
)

target_include_directories(${SF_EXTENSIONS_LIB}
    PRIVATE
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})