#!/bin/bash
METHOD=Gauss
#METHOD=RungeKutta
#METHOD=Falberg
#METHOD=DormanPrince
#METHOD=LStableDiagonal
#METHOD=Lobatto
#METHOD=Rado
ORDER=6
WAY=1
#ITERATION=SI
ITERATION=Zeidel
#ITERATION=Newton
APPROX=0.001
TASK_ORDER=2
TEST_NUM=_tough3
#TEST_NUM=1

./bin/main -m $METHOD -o $ORDER -w $WAY -i $ITERATION -a $APPROX < test/ODUTest/Order$TASK_ORDER/test$TEST_NUM.txt
#./bin/main -m $METHOD -o $ORDER -w $WAY -i $ITERATION -a $APPROX  --system < test/ODUTest/Systems/test$TEST_NUM.txt
pdflatex --output-directory=report report/report.tex