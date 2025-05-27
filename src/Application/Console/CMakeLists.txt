#build_app(ODUSolverConsole main.cpp ${INPUT_LIB})
build_app(crank crankTest.cpp "${INPUT_LIB};${CRANK_NICOLSON_SOLVER_LIB_2};${GUI_LIB}")
build_app(experimental experimental.cpp "${INPUT_LIB};${CRANK_NICOLSON_SOLVER_LIB};${GUI_LIB}")