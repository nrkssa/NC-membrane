#ifndef __DECLARATIONS_H__
#define __DECLARATIONS_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include "_constant.h"
#include <limits>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <iterator>
#include <typeinfo>
#include <array>
#include "MersenneTwister.h"
using namespace std;
typedef std::numeric_limits< double > dbl;

//@ structure to store normal and householder matrix
struct tmatrix{
  double normal[3][1], hhm[3][3];
};

struct runparam{
  int curr_proc, total_proc;
  double nmcss, nmcse;
  int win_start, win_end, step_start, step_end, sampling_cutoff;
  string system, mpi_mode;
  int samp_freq, osdp_freq, hist_freq, antab_freq, conf_freq, ant_write_freq;
  int DTMCsteps;
  bool DTMCflag;
};


#endif
