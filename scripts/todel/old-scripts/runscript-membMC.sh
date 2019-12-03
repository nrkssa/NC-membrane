#!/bin/bash
rm -rf RUNDIR
echo "Making directory named test for the current run"
mkdir RUNDIR
cp MEMB_MC RUNDIR/ -v
cd RUNDIR
echo "Run Starts"
echo "==============================================================================================" 
echo
echo
echo
mpirun -np $1 ./MEMB_MC $2


