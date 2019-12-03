#ifndef __INTERACTION_PARAMETERS_H__
#define __INTERACTION_PARAMETERS_H__
class _intparam{
int nabtype, nanttype;
double **kbond, **delg, **delr, **delrsq, min_delr, max_delr;
double param[3], L;
public: _intparam(){}
public: ~_intparam(){}
public: void init();
public: void dump_parameters(ofstream &ofile);
public: void read_parameters(ifstream &ifile);
public: void set_array_sizes();
public: void reset_box_dimensions();
public: void read_interaction_parameters();
public: void set_minmax_reaction_length();
public: double get_kbond(int i, int j);
public: double get_kbond(int ncnum, int abnum, int antnum);
public: double get_delg(int i, int j);
public: double get_delg(int ncnum, int abnum, int antnum);
public: double get_delr(int i, int j);
public: double get_delrsq(int i, int j);
public: double get_delr(int ncnum, int abnum, int antnum);
public: double get_delrsq(int ncnum, int abnum, int antnum);
public: double* get_bellbond_parameters(int ncnum, int abnum, int antnum);
public: double get_min_reaction_length();
public: double get_max_reaction_length();
};
#endif
