build_app(ODUSolverConsole main.cpp ${INPUT_LIB})
build_app(crank crankTest.cpp "${INPUT_LIB};${GUI_LIB}")