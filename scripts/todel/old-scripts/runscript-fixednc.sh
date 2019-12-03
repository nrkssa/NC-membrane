#!/bin/bash
rm -rf RUNDIR
echo "Making directory named test for the current run"
mkdir RUNDIR
cp NC_MEMB_fixedNC RUNDIR/ -v
cd RUNDIR
echo "Run Starts"
echo "==============================================================================================" 
echo
echo
echo
./NC_MEMB_fixedNC


