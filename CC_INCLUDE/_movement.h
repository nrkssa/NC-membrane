#ifndef __MOVEMENT_H__
#define __MOVEMENT_H__
class _movement {
struct vertex *verptr;
struct antigen *antptr;
struct triangle *triptr;
struct nanocarrier *ncpointer;
double chem_cutoffu, chem_cutoffl, soft_radius, ant_size;
double max_ves_antigentip_dist;
double Dreaction,flexural_stiffness,spring_constant,eqm_bond_energy;
double ves_radius;
double tempx, tempy, tempz;
double rrx, rry, rrz;
double temperature,beta;
double energy_shear, energy_rot;
double dl,l_max, angle_max, delE;
int accept;
int accept_tcnob,accept_tcb,accept_tvnob,accept_tvb,accept_rot;
int num_tcnob,num_tcb,num_tvnob,num_tvb,num_rot,num_samp,M;
int *lscl, *head, *cell,*mlscl,*mhead,*mcell;
int sel_antigen, sel_ves, sel_ab, mem;
int num_link1cells, N_point, frame;
int *temp_bonded, temp_bond[2], temp;
char s7[10];
double shear_rate;
ofstream outfile7;
double memb_nc_cutoff_distance,max_AB_ANT_distance,max_bond_formation_distance;
double max_nc_memb_distance,max_nc_anttip_distance;
double *xpn,*ypn,*zpn, *Eo,*w, *thetat,*phit;
int number_of_move, *b,*bond_status;
double nbonded,nbreakage;
char bias_mode;
int current_multivalency, num_nc, num_antigens, nver, pbnum, ntr, nlink;

public: _movement(){}
public: void init();
public: void reset_bond_and_equilibrate(int nsteps, bool DTMCflag);
public: void translation(char ch);
public: void rotation(char ch);
public: int getmembrane_farthestz(int i,char ch);
public:	int getmembrane_farthestvertex(int i,char ch);
public: void check_antigen_antibody_dist(char *s, int ves);
public: void write_vtkfiles(int confno);
public: void bond_formation(int timestep);
public: double calc_energy(int mem_sel,double dist,double flex_dist, int bound);
public: int select(double *w,double sumw);
public: int get_accept();
public: double get_delE();
public: int getaccept_tcnob();
public: int getnum_tcnob();
public: int getaccept_tcb();
public: int getnum_tcb();
public: int getaccept_tvnob();
public: int getnum_tvnob();
public: int getaccept_tvb();
public: int getnum_tvb();
public: int getaccept_rot();
public: int getnum_rot();
public: int getN_point();
public: void reset_bond_counter();
public: int* get_bond_status();
};
#endif
