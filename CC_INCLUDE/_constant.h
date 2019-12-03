#ifndef __CONSTANT_H__
#define __CONSTANT_H__
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <random>
#include <chrono>
#include <omp.h>
using namespace std;


//#define MAX_ENTRIES 10
#define MAX_CHAR 512
#define MAX_FILE_ENTRY 100

//@ All length in C++ are defined in base units of nm
#define base_length 1E-9
#define base_length_sq 1E-18

#define BUFFER_SIZE_SMALL 1048576                                                  // 1MB buffer   
#define BUFFER_SIZE_LARGE 10485760                                                 // 10MB buffer 

#define PI 3.141592653589793238462643                                              // pi
#define TWOPI 2*3.141592653589793238462643                                         // 2pi

#define sqrttwo  sqrt(2.0)
#define sqrtthree sqrt(3.0)
#define oneosqrt2 (1.0/sqrttwo)
#define onemoneosqrt2 (1.0-oneosqrt2)

#define kb 1.3806503e-23     // Boltzmann constant in SI units
#define Na 6.0221415e23;     // Avagadro's Number	
#define R 8.31447;           // Universal gas constant
#define COMPARE_TOL 1.0                                       // do not set it too low. The floating point errors will show up in vesicle rotations.

#define max_curvature_binned 1.0
#define maxantigen_ver 10
#define maxantigen_tri 30

#define str_mean_curv "mean_curv"
#define  str_distance "distance"
#define  str_rdf "rdf"
#define  str_current_block "current_block"
#define  str_running_avg "running_average"
#define  advance_dir "advance"
#define  recede_dir "recede"
#define  xy_proj "xy"

#define  str_red "\033[1;31m"                                                   // strings to display screen output in red,blue,green and black
#define  str_blue "\033[1;34m"
#define  str_green "\033[2;32m"
#define  str_black "\033[0m"
#define  endmarker "\n\t\t\t =============================="


//typedef struct{
// char* token[MAX_ENTRIES];
// int n;
//} parser;

#endif
