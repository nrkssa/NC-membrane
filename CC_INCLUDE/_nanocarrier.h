#ifndef __NANOCARRIER_H__
#define __NANOCARRIER_H__

class _ANTIBODY{
int nabtype;
double *ablength, *abradius,*ab_type_conc;
public: _ANTIBODY(){}
public: void init(int index);
public: void init(int nab, double* abl, double* abr, double* abc);
public: void print_message();
public: int get_numtype();
public: double get_abradius(int i);
public: double get_ablength(int i);
public: void parse_abparameters(int index);	
};

class _NANOCARRIER {
struct vertex *verptr;
struct nanocarrier *ncpointer;
double xc, yc, zc, anglex, angley, anglez, radius, soft_radius;
double xcold, ycold, zcold, anglexold, angleyold, anglezold;
double L, H, temperature, kBT;
double *x_antibody, *y_antibody, *z_antibody,*xtold,*ytold, *ztold;
double *x_npb_antibody, *y_npb_antibody, *xf_antibody, *yf_antibody, *zf_antibody;
double *xbold,*ybold,*zbold;
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
_ANTIBODY abclass;

//@ Biasing modes
char bias_mode;
double kbias,biasref,binsize,lambda;
double  rstart, hstart, rref, href;
double meanr[3],biasparm[4];

public: _NANOCARRIER(){}
public: void init(int i);
public: void init(int ncindex,ifstream &ifile);

public: void dump_conf(ofstream &ofile);
public: void reset_box_dimensions();

public: void print_address();
public: void parse_nc_parameters();
public: void setnum_ab();
public: void update_fortran_pointers();
public: double compute_max_soft_radius();
private: void create_antibody();
private: void shift_to_COM_frame(int nab, double *xpos, double *ypos, double *zpos);
private: void set_base_tip_antibody(int i, double xpos, double ypos, double zpos);
public: void  compute_ves_AB_distance();
public: double apply_PBC(double rval);
public: double shift_PBC(double rval);
public: double compute_bond_length(int j, int k);

//@
public: char get_bias_mode();
public: int get_num_abtype();
public: double* get_shadowr();
public: double  get_shadowr_dist();
public: double get_shadowcurv();
public: int get_nshadowvert();
//@
public: void write_NC_configuration(string filename);
public: void shift_nc_position(double xpos, double ypos, double zpos);
public: void setL (double L1);
public: void setH (double H1);
public: void setradius(double radius1, double soft_radius1);
public: double getradius();
public: double getsoft_radius();
//@
public: void store_configuration();
public: void restore_configuration
	();
public: void setrc(double xc1, double yc1, double zc1);
public: void setxc(double xc1);
public: void setyc(double yc1);
public: void setzc(double zc1);
public: double getxc();
public: double getyc();
public: double getzc();
//@
public: void setanglex(double ax1);
public: void setangley(double ay1);
public: void setanglez(double az1);
public: void reset_euler_angles();
public: double getanglex();
public: double getangley();
public: double getanglez();
//@
public: int getnum_ab();
public: int* get_ab_types();
public: int get_ab_type(int k);
//@
public: void setxy_npb_antibody();
public: double getx_antibody(int k);
public: double gety_antibody(int k);
public: double getz_antibody(int k);
public: double getxf_antibody(int k);
public: double getyf_antibody(int k);
public: double getzf_antibody(int k);
public: double getx_npb_antibody(int k);
public: double gety_npb_antibody(int k);
public: double getxf_npb_antibody(int k);
public: double getyf_npb_antibody(int k);

//@
public: bool getbonded(int k);
public: void setbonded(int k, bool temp);
public: int getbondedto(int k);
public: void setbondedto(int k,  int temp);
public: bool does_bond_breaks(double tempx,double tempy,double tempz);
public: bool does_bond_breaks();
public: bool Does_Bond_Breaks();

//@
public: void translate (double tempx, double tempy, double tempz);
public: void rotate (double phi, double theta, double psi,bool transpose=false);
public: void rotate_mp (double phi, double theta, double psi,bool transpose=false);
public: void translate_OMP (double tempx, double tempy, double tempz);
public: void rotate_ab (double phi, double theta, double psi, int k);

//@
public: double get_rstart();
public: double get_hstart();
public: double get_kbias();
public: double get_biasref();
public: double get_lambda();
public: string get_biasdir();
public: double get_biasbinsize();
public: double* get_biasparameters();
public: void set_biasref(double refval);
public: void shift_biaswindow();

//@ Function for counters and step ab_size
public: void reset_counters(string cname);
public: void increment_counters(string cname, bool acceptance);
public: int* get_counter(string cname);
public: double* get_mc_stepsizes();
public: void reset_step_size(int mode);
public: double get_step_size(string cname);

//@ Function to handle multivalency calculations
public: void set_multivalency(int mval);
public: int get_multivalency();
public: void update_bond_counters(int dmval, int dnbon, int dnbr);
public: int* get_bond_counters();
public: void reset_bond_counters();
};


//@ This is a public class that interacts with the rest of the code

class _nanocarrier {
  
// Pointers to Fortran data
struct vertex *verptr;
struct antigen *antptr;
struct triangle *triptr;

//int gcount=0;
int num_vesicles, num_antigens, num_memb_vert, num_samp=100;
double L, H, temperature, beta, kBT;
int temp[2], *b;
double angle[2], dirNormal[3], *thetat, *phit, *Eo, *weight;
double *xpn, *ypn, *zpn;

_NANOCARRIER *v;

public: _nanocarrier(){}
public: void init();
public: void init(ifstream &ifile);
public: _nanocarrier(double L1, double H1, int num_vesicles1);
public: void set_array_sizes();
public: _NANOCARRIER getvesicle(int j);
public: void reset_box_dimensions();
public: void dump_conf(int i, ofstream &ofile);
public: void read_conf(int i, ifstream &ifile);

public: void print_address(int i);
public: int getnum_members();
public: double getL();
public: double getH();
public: void setvesicle (_NANOCARRIER v_temp, int i);

public: double getxc(int i);
public: double getyc(int i);
public: double getzc(int i);

public: void setxc(int i, double x);
public: void setyc(int i, double y);
public: void setzc(int i, double z);

public: double getanglex(int i);
public: double getangley(int i);
public: double getanglez(int i);

public: void setanglex(int i, double x);
public: void setangley(int i, double y);
public: void setanglez(int i, double z);

public: void setxy_npb_antibody(int ves_num);

public: void setradius(int i, double radius, double soft_radius1);
public: double getradius(int i);
public: double getsoft_radius(int i);

public: void setnum_ab(int i, int n);
public: int getnum_ab(int i);
public: int* get_ab_types(int i);
public: int get_ab_type(int i, int k);

public: void translate(int i, double tempx, double tempy, double tempz);
public: void rotate(int i, double theta, double psi, double phi, bool transpose=false);

public: void store_configuration(int i);
public: void restore_configuration(int i);
public: double getx_antibody(int i, int k);
public: double gety_antibody(int i, int k);
public: double getz_antibody(int i, int k);
public: double getxf_antibody(int i, int k);
public: double getyf_antibody(int i, int k);
public: double getzf_antibody(int i, int k);
public: double getx_npb_antibody(int i, int k);
public: double gety_npb_antibody(int i, int k);
public: double getxf_npb_antibody(int i, int k);
public: double getyf_npb_antibody(int i, int k);

public: bool getbonded(int i, int k);
public: bool getbonded(int i);
public: void setbonded(int i, int k, bool temp);
public: int getbondedto(int i, int k);
public: void setbondedto(int i, int k,  int temp);
public: bool does_bond_breaks(int i, double tempx, double tempy, double tempz);
public: bool does_bond_breaks(int i);
public: bool Does_Bond_Breaks(int i);
public: void updatebonds(int i);
public: void reset_antab_bonds(int i);

public: void  compute_ves_AB_distance(int );
public: double compute_bond_length(int i, int j, int k);
public: void reset_euler_angles(int i);

//@
public: void store_restore_shadow(int flag, int ncnum);
public: void compute_shadow(int ncnum);
public: double* get_shadowr(int i);
public: double get_shadowr_dist(int i);
public: double get_shadowcurv(int i);
public: int get_nshadowvert(int i);

public: int get_num_abtype(int i);
public: void shift_nc_position(int i, double xpos, double ypos, double zpos);
public: void write_NC_configuration(int i,string filename);

//@ Function for counters and step size maintenance
public: void reset_counters(int i, string cname, bool allflag = false);
public: void increment_counters(int i, string cname, bool acceptance);
public: int* get_counter(int i, string cname);
public: double get_step_size(int i, string cname);
public: double* get_mc_stepsizes(int i);

//@ Functions to handle multivalencies
public: int get_multivalency(int i);
public: void set_multivalency(int i, int mval);
public: void update_bond_counters(int i, int dmval, int dnbon, int dnbr);
public: int* get_bond_counters(int i);
public: void reset_bond_counters(int i);

//@ MC moves
public: void translation();
public: void rotation();
public: void rotation_mpparallel();
public: void bond_formation();
public: void bond_formation_mpaccelerate();
public: void bond_formation_mpparallel();
public: tmatrix construct_householder_matrix(int sel_antigen);
public: double calc_energy(int sel_ves, int sel_ab, int sel_ant, double dist, double flex_dist, int bound);
public: int select_from_Rosenbluth(double *w, double sumw);

//@
public: char get_bias_mode(int i);
public: double get_rstart(int i);
public: double get_hstart(int i);
public: double get_kbias(int i);
public: double get_biasref(int i);
public: double get_lambda(int i);
public: string get_biasdir(int i);
public: double get_biasbinsize(int i);
public: double* get_biasparameters(int i);
public: void set_biasref(int i, double refval);
public: bool is_calculation_TI();
public: void shift_biaswindow();
};
#endif
