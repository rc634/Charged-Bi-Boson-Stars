The file QBSParams.hpp contains the physical (and some numerical) parameters you may wish to change. After chaning the code must be compiled before running. When the code successfuly runs the plotting python script can be used.

The default example has Qtotal/Ntotal = -3.03288411907781e-06 which corresponds to 99.9994349308342% of the charge being self cancelled out. Quite close to a neutral obeject.

Note: to run the .sh files the access permissions may need relaxing with

chmod 777 *.sh

 
COMPILE WITH

./compile


COMPILE AND RUN WITH

./compile_and_run

COMPILE, RUN AND PLOT WITH (ASSUMING PYTHON3 COMMAND WORKS)

./compile_run_plot


------ Alternatively see below --------

COMPILE THE CODE WITH

g++ -O3 -o executable main.cpp


RUN THE C++ CODE WITH

./executable


GRAPH PLOT SOLUTION WITH PYTHON3 (I THINK AVOID PYTHON2)

python3 plot.py

or maybe

python plot.py
