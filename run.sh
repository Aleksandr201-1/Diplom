#!/bin/bash
#METHOD=Gauss
METHOD=RungeKutta
#METHOD=Falberg
ORDER=4
WAY=1
ITERATION=SI
APPROX=0.01
TASK_ORDER=2
TEST_NUM=_tough3
#TEST_NUM=14

./bin/main -m $METHOD -o $ORDER -w $WAY -i $ITERATION -a $APPROX < test/Order$TASK_ORDER/test$TEST_NUM.txt; pdflatex --output-directory=report report/report.tex