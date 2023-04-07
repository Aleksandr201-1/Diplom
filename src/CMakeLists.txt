#if(${PLATFORM} EQUAL ${WINDOWS})
#set(MAIN main)

add_executable(main main.cpp)

#set(BOOST_LIBS boost_iostreams boost_system boost_filesystem)
#set(BOOST_LIBS )

add_subdirectory(General)
add_subdirectory(PDFReporter)
add_subdirectory(ODUSolver)
add_subdirectory(ChemicalGenerator)

target_link_libraries(main
    PRIVATE
        ${GENERAL_LIB}
        ${PDF_REPORTER_LIB}
        ${ODU_SOLVER_LIB}
        ${CHEMICAL_GENERATOR_LIB}
        #${BOOST_LIBS}
)

target_include_directories(main
    PUBLIC
        ${MODULES_INCLUDE_PATH}
        ${CMAKE_CURRENT_LIST_DIR}
)

add_definitions(${COMPILE_FLAGS})