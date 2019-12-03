#ifndef __DECLARATIONS_H__
#define __DECLARATIONS_H__
#include <cmath>
#include <limits>
#include <array>
#include <omp.h>
#include "_constant.h"
using namespace std;
#include "MersenneTwister.h"
typedef std::numeric_limits< double > dbl;

//@ structure to store normal and householder matrix
struct tmatrix{
  double normal[3][1], hhm[3][3];
};

struct runparam{
  int curr_proc, total_proc, ensno;
  double nmcss, nmcse;
  int win_start, win_end, step_start, step_end, sampling_cutoff;
  string system, mpi_mode;
  int samp_freq, osdp_freq, antab_freq, conf_freq, ant_write_freq, ant_traj_freq;
  int DTMCsteps;
  bool DTMCflag;
};

struct ncstruct{
	double xc, yc, zc, anglex, angley, anglez, radius, soft_radius;
	double L, H, temperature, kBT;
	double *x_antibody, *y_antibody, *z_antibody;
	double *x_npb_antibody, *y_npb_antibody, *xf_antibody, *yf_antibody, *zf_antibody;
	double *xf_npb_antibody, *yf_npb_antibody;
	bool *bonded;                         
	int num_ab, num_abtype, index;
	double *ab_size, *ab_radius, *ab_type_conc, *soft_rad_ab;
	string ncshape,ab_arrangement,biasdir;
	int *bondedto, *ab_type, *tnum_ab_type;
	
	//@Counters for translation and rotation
	int trans_counter[4], rot_counter[4];
	double stepsize[4];
	int multivalency, nbonded, nbroken, bond_counters[3];
	double tstep_ub, tstep_b, rstep_ub, rstep_b, adjust_interval;                                               // adaptive step size
	
	//@ AB data
	int nabtype;
	double *ablength, *abradius;
	
	//@ Biasing modes
	char bias_mode;
	double kbias,biasref,binsize,lambda;
	double  rstart, hstart, rref, href;
	double meanr[3],biasparm[4];
};

#endif
