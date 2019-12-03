#!/bin/bash
rm -rf RUNDIR
echo "Making directory named test for the current run"
mkdir RUNDIR
cp TI_CALC RUNDIR/ -v
cd RUNDIR
echo "Run Starts"
echo "==============================================================================================" 
echo
echo
echo
mpirun  -np $1 ./TI_CALC $2
