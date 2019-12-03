#ifndef __POTENTIAL_H__
#define __POTENTIAL_H__
class _potential {
double L ,H, temperature, kBT;
struct vertex *verptr;
struct antigen *antptr;
struct triangle *triptr;
int num_nc, num_antigens, nver, pbnum, ntr, nlink;


public: _potential(){}
public: void init();
public: void reset_box_dimensions();
public: double dist_receptor_ligand(int ncnum, int abnum, int recep_num);
public: double dist_nc_receptor(int ncnum, int recep_num);
public: double calc (int mem_sel, char ch);
public: double calc (int mem_sel);                               //The reaction energy alone : Ram -30July12
public: double glycocalyx(int mem_sel);                          //compute energy due to glycocalyx int. with the immersed nano carrier
public: double reaction_energy(int ncnum);                       // get reaction energy for nanocarrier ncnum
public: double reaction(double dist, int ncnum, int abnum, int antnum, bool kBTconversion=false);
public: double change_reaction(double dist1, double dist2, int ncnum, int abnum, int antnum, bool kBTconversion=false);
public: double antigen_flexure_energy(int ant_num);
public: double total_energy();                                   //to compute the total energy of the system:Jin
public: double total_Renergy();                                  //to compute the total reaction energy of the system:Jin
public: double total_Benergy();                                  //to compute the total bending energy of the system:Jin
public: void doglyx_calc(int flag);                              //to set a flag to do glyx calc.
public: double compute_distance(double *pos1, double *pos2);     // compute the geometric distance between  pos1 and pos2
public: double compute_biasing_potential(int mem_sel, double valold, double valnew);
public: double compute_biasing_potential(int mem_sel, double *rold,double *rmeanold,double *rnew,double *rmeannew);
};
#endif
