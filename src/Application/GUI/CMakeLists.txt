build_app(ODUSolverGUI GUITest.cpp ${GUI_LIB})
build_app(gui_chem gui_chem.cpp "${INPUT_LIB};${GUI_LIB}")