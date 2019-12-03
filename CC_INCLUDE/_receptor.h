#ifndef __RECEPTOR_H__
#define __RECEPTOR_H__

class _RECEPTOR {
struct vertex *verptr;
struct antigen *antptr;
struct triangle *triptr;
double xc, yc, zc;
double xt, yt, zt, beta;                                         // the antigen tip positions:Jin
float *trbx, *trby, *trbz, *trtx, *trty, *trtz, *tangt, *tangp;                                            // Store the tip and base position as a function of MC time
int  *tbond, *tblig, *tbnc;                                      // Store the base position and ligand pos as a function of time     
double L, radius, ant_length;                                    // length of cell
double memb_nc_cutoff_distance;
bool bonded;
int index ;
int bondedtoves, bondedtoab;
int type;
int ntsamples,ctsamples;

public: _RECEPTOR(){}
public: void init(int i);
public: void init(int i1, double* fv, int* i2, double* fv1, int* i3);
public: void reset_box_dimensions();
public: void dump_conf(ofstream &ofile);
public: void initialize_storage();
public: void initialize_storage_variables();
public: void store_trajectory(int window);
public: void dump_trajectory(int window);
public: void setradius(double radius1);
public: void setxc(double xc1);
public: void setyc(double yc1);
public: void setzc(double zc1);
public: double getradius();
public: double getlength();
public: int gettype();
public: double get_flexure();
public: double getxc();
public: double getyc();
public: double getzc();
public: double getxt();
public: double getyt();
public: double getzt();
public: void setxt(double xt1);
public: void setyt(double yt1);
public: void setzt(double zt1);
public: void set_theta(double theta);
public: void set_phi(double phi);
public: double get_theta();
public: double get_phi();
public: void setantigentip_position();
public: void setL(double L1);
public: bool getbonded();
public: void setbonded(bool temp);
public: int getbondedtoves();
public: void setbondedtoves(int temp);
public: int getbondedtoab();
public: void setbondedtoab(int temp);
public: void antigen_translation();
public: bool does_bond_breaks(double tempx,double tempy,double tempz);
public: bool does_bond_breaks();
public: void rotate(double phi, double theta, double psi);
};

class _receptor {
struct vertex *verptr;
struct antigen *antptr;
int num_vesicles, num_antigens, num_memb_vert;
int type_ant, *num_type_ant;
double *radius_ant, *length_ant, *kflex;
string ant_pattern;
double L, H;
double Dreaction, chem_cutoffu, chem_cutoffl, ant_length, soft_radius;
int temp[2];
double angle[2];
double dirNormal[3];
_RECEPTOR *c;
public: _receptor(){}
public: void init();
public: void init(ifstream &ifile);
public: _receptor(double L1, double H1, int num_antigens1);
public: void set_array_sizes();
public: void reset_box_dimensions();
public: void assign_antigen_types();
public: void store_trajectory(int window);
public: void dump_trajectory(int window);
public: void dump_conf(ofstream &ofile);
public: void read_conf(ifstream &ifile);
public: _RECEPTOR getantigen(int j);
public: int getnum_members();
public: double getL();
public: double getH();
public: void setantigen (_RECEPTOR c_temp, int i);
public: void setxc(int i, double x);
public: double getxc(int i);
public: void setyc(int i, double y);
public: double getyc(int i);
public: void setzc(int i, double z);
public: double getzc(int i);
public: double getantigenxt(int i);
public: double getantigenyt(int i);
public: double getantigenzt(int i);
public: double getxt(int i);
public: double getyt(int i);
public: double getzt(int i);
public: void setantigenxt(int i, double xt1);
public: void setantigenyt(int i, double yt1);
public: void setantigenzt(int i, double zt1);
public: void setantigenrt(int i, double xt1, double yt1, double zt1);
public: void setantigentip_position(int i);
public: double getradius(int i);
public: double getlength(int i);
public: int gettype(int i);
public: double get_flexure(int i);
public: void setradius(int i, double radius);
public: void antigen_translation(int i);
public: bool getbonded(int i);
public: void setbonded (int i, bool temp);
public: int* getbondedto(int i);
public: void setbondedto(int i,  int* temp);
public: void reset_antigen_bond(int i);
public: bool does_bond_breaks(int i, double tempx,double tempy,double tempz);
public: bool does_bond_breaks(int i);
public: void set_theta_phi(int i, double theta, double phi);
public: void set_theta(int i,double theta);
public: void set_phi(int i,double phi);
public: double* get_theta_phi(int i);
public: double get_theta(int i);
public: double get_phi(int i);
public: void updatebonds(int i);
public: void receptor_hopping();
public: void receptor_diffusion();
};
#endif
