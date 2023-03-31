#if(${PLATFORM} EQUAL ${WINDOWS})
#set(MAIN main)

add_executable(main main.cpp)

#set(BOOST_LIBS boost_iostreams boost_system boost_filesystem)
#set(BOOST_LIBS )

add_subdirectory(General)
add_subdirectory(PDF-reporter)
add_subdirectory(Tough)
add_subdirectory(NotTough)
add_subdirectory(ChemicalGenerator)

target_link_libraries(main
    PRIVATE
        GeneralLIB
        ReportLIB
        #ToughLIB
        NotToughLIB
        ChemicalGenLIB
        #${BOOST_LIBS}
)

add_definitions(${COMPILE_FLAGS})