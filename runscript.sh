#!/bin/bash
rm -vf RUNDIR/*
mkdir RUNDIR
cp $2 RUNDIR/ -v
cd RUNDIR
echo "Run Starts"
echo "==============================================================================================" 
echo
echo
echo
mpirun  -np $1 ./$2 $3 $4 $5
