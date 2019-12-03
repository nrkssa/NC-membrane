#!/bin/bash
rm -rf RUNDIR
echo "Making directory named test for the current run"
mkdir RUNDIR
cp $1 RUNDIR/ -v
cd RUNDIR
echo "Run Starts"
echo "==============================================================================================" 
echo
echo
echo
mpirun  -np $2 ./$1 $3
