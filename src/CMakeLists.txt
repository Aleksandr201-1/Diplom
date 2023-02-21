#if(${PLATFORM} EQUAL ${WINDOWS})
#set(MAIN main)

add_executable(main main.cpp)

#set(COMPILE_FLAGS -Wextra -Wall -Werror -Wno-sign-compare -Wno-unused-result -O3)
set(COMPILE_FLAGS -Wextra -Wall -Wno-sign-compare -Wno-unused-result)



#set(BOOST_LIBS boost_iostreams boost_system boost_filesystem)
set(BOOST_LIBS )

add_subdirectory(General)
add_subdirectory(PDF-reporter)
add_subdirectory(Tough)
add_subdirectory(NotTough)

target_link_libraries(main
    PRIVATE
        GeneralLIB
        ReportLIB
        #ToughLIB
        NotToughLIB
        ${BOOST_LIBS}
)

add_definitions(${COMPILE_FLAGS})