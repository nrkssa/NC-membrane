#ifndef __DATAINPUT_H__
#define __DATAINPUT_H__
#include "_declarations.h"
class _datainput {

//@ parameters for the antigens (recptors)  
int num_antigens, type_ant, *num_type_ant; 
double *conc_ant, *radius_ant, *length_ant, *kflex;
double min_antlen, max_antlen;
double kBT, binsize_rdf;
std::string pattern_ant, run_system, mpi_mode;
int k, num_moves, num_dtmc_moves, sampling_interval, N, iRES, nver, ntr, num_linkcells, hist_samp_freq;
double periodic_box_length, periodic_box_height;
int num_nanocarriers, num_eq_moves;
int win_start, win_end;
double periodic_box_length_height_ratio, memb_init_bond_length, membrane_init_z;
double temperature;
int window_conf_interval, window_antab_data_interval, osd_print_interval, antigen_memb_write_interval,antigen_traj_interval;
runparam rparam;

public:
  _datainput() {}
  void init();
  //@
  void read_antigen_parameters();
  void dump_antigen_parameters(ofstream &ofile);
  void read_antigen_parameters(ifstream &ifile);
  void print_message_ant();
  void read_system_geometry();
  void dump_system_geometry(ofstream &ofile);
  void read_system_geometry(ifstream &ifile);
  double get_max_ncsize();
  int get_types_of_antigens();
  int* get_number_of_antigen_types();
  double* get_antigen_radius();
  double get_antigen_radius(int i);
  double* get_antigen_length();
  double get_antigen_length(int i);
  double* get_antigen_flexure();
  string get_antigen_pattern();
  double get_antigen_maxlength();
  double get_antigen_minlength();
  //@
  void setup_runparameters();
  void dump_run_parameters(ofstream &ofile);
  void read_run_parameters(int ensno, int window, ifstream &ifile);
  runparam get_runparameters();
  void print_runparam_message();
  
  //@
  void check_periodicbox_dimensions();
  void resetperiodic_box_dimensions();
  double getperiodic_box_length();
  double getperiodic_box_height();
  void resetperiodic_box_dimensions(double boxlength);
  
  int getnum_members(char ch);
  void setnum_members(char ch, int nval);
  int getnum_moves();
  int getnum_eq_moves();
  double gettemperature();
  double get_kBT();
  double getmembrane_init_z();
  int getN();
  void setN(int N_new);
  int getsampling_interval();
  double get_binsize_rdf();
  int get_window_conf_interval();
  int get_window_antab_data_interval();
  int get_antigen_memb_write_interval();
  int get_antigen_traj_interval();
  int get_antigen_traj_samples();
  int get_antigen_conf_interval();
  int get_osd_print_interval();
  double get_window_step_size();
  double get_antigen_vesicle_cutoff_dist();
  string get_approach_dir();
  void initialize_fortran_variables();
};
#endif
