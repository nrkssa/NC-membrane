#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__
class _distribution{ 
char ch;
double *hist,L,H;
double *rdf,nsample_rdf;
double zmin,zmax,window_size,hist_binsize;
int Nbin_rdf;
double binsize_rdf;                                                                //number of bins
double reaction_dist_sq;
int num_members;
double radius,soft_radius,rdf_avgnum_antigens,rdf_antigen_cutoff,avg_antigen_density,antigen_length;
double factor, ig;
int radial_counter;
struct vertex *verptr;
struct antigen *antptr;
ofstream outfile;
public: _distribution(){}
public: void init ();
public: void reset_box_dimensions();
public: void assign_binvalue(double **dist_name, int nbin,  int column, double minval, double maxval, double binsize);
public: void assign_binvalue(double **dist_name, int nbin,  int column, double value);
public: void assign_binvalue(double *dist_name, int nbin,  double value);
public: void transfer_binvalues(double **target, double **source, int nbin,  int column);
public: void compute_antigen_rdf_distribution(int nbin, double binsize, double *hist_array,double *nsample);
public: void write_antigen_rdf_distribution(int nbin, double binsize, double *hist_array,double nsample,std::ostringstream &ss);
public: double array_sum (double **hist_array, int nbin, int column);
public: void compute_rdf_distribution();
public: void write_rdf_distribution(int rank,int frame);
public: int  compute_multivalency(int ves_num);
public: void check_antigen_vesicle_dist(int ves_num,int stepno,double type_of_move, double type_of_trans,double ratio_m, double flip_move_ratio);
public: bool check_bondspring_dist(int mem_sel,char ch);
public: void initialize_zhistogram();
public: void initialize_antigen_distribution_histogram();
public: void calculate_zhistogram(double zvalue,double meanz);
public: void write_zhistogram(int rank,int frameno);
public: double** compute_histogram(double **array1,double **array2, int nrow, int ncol, double maxval, double minval, double binsize, int nbin);
public: double** compute_histogram(double **array1,int nrow, int ncol, double maxval, double minval, double binsize, int nbin);
public: double dist_receptor_ligand(int ncnum, int abnum, int recep_num);	
public: double dist_nc_receptor(int ncnum, int recep_num);
public: double  dist_nc_nc(int ncnum1, int ncnum2);
};
#endif
