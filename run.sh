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
APPROX=0.00001
TASK_ORDER=2
#TEST_NUM=_tough3
#TASK_TYPE=KoshiSystem
#TASK_TYPE=Chemical
TEST_NUM=1

#run -m Gauss -o 6 -w 1 -i Zeidel -a 0.001 -tt Chemical < test/ChemicTest/test1.txt
#--track-origins=yes
#valgrind ./bin/main -m $METHOD -o $ORDER -w $WAY -i $ITERATION -a $APPROX -tt Koshi < test/ODUTest/Order$TASK_ORDER/test$TEST_NUM.txt
#valgrind ./bin/main -m $METHOD -o $ORDER -w $WAY -i $ITERATION -a $APPROX -tt KoshiSystem < test/ODUTest/Systems/test$TEST_NUM.txt
./bin/main -m $METHOD -o $ORDER -w $WAY -i $ITERATION -a $APPROX -tt Chemical < test/ChemicTest/test$TEST_NUM.txt
pdflatex --output-directory=report report/report.tex