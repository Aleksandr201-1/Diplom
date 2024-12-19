function(add_module_lib TARGET_NAME)
    set(OPTIONS "")
    set(ONE_ARG "")
    set(MULTI_ARG PRIVATE PUBLIC DEPENDS)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "${OPTIONS}" "${ONE_ARG}" "${MULTI_ARG}")
    if (NOT "${ARG_UNPARSED_ARGUMENTS}" STREQUAL "")
        message(FATAL_ERROR "add_module_lib: extra arguments are ${ARG_UNPARSED_ARGUMENTS}")
    endif()
    #if (NOT ARG_PRIVATE)
    #    message(FATAL_ERROR "add_module_lib: PRIVATE sources not provided")
    #endif()
    add_library(${TARGET_NAME} "")
    if (ARG_PRIVATE AND ARG_PUBLIC)
        target_sources(${TARGET_NAME}
            PRIVATE
                ${ARG_PRIVATE}
            PUBLIC
                ${ARG_PUBLIC}
        )
    endif()
    target_link_libraries(${TARGET_NAME} 
        PRIVATE
            ${SFML_LIBS}
            ${GENERAL_LIB}
            ${ARG_DEPENDS}
    )
    target_include_directories(${TARGET_NAME}
        PRIVATE
            ${MODULES_INCLUDE_PATH}
            ${CMAKE_CURRENT_LIST_DIR}
    )
    add_definitions(${COMPILE_FLAGS})
endfunction()

#libraries
add_subdirectory(General)
add_subdirectory(Math)
add_subdirectory(Butcher)
add_subdirectory(NumericMethods)
add_subdirectory(ODUSolver)
add_subdirectory(PDFReporter)
add_subdirectory(Input)
add_subdirectory(SFExtensions)
add_subdirectory(GUI)

#executable
add_subdirectory(Application)