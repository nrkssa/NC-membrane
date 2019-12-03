#ifndef __SELFAVOIDANCE_H__
#define __SELFAVOIDANCE_H__
class _selfavoidance {
double *radius, *soft_radius, *ant_radius;
int *cell, num_nc, num_antigens; 						// nearest linked cells array
double tempx, tempy, tempz; 							// temp variables to hold position
double L ,H;
double **antigen_vesicle_cutoff_dist;
double **antigen_antigen_cutoff,**vesicle_vesicle_cutoff,*vesicle_membrane_cutoff;
struct vertex *verptr;
struct antigen *antptr;
struct triangle *triptr;

public: _selfavoidance(){}
public: void init();
public: void reset_box_dimensions();
public: double sa_distance(char ch1, char ch2, int m1, int m2);
public: bool overlap (int mem_sel, char ch); 	 			                        // Checking overlap of one vesicle with all other vesicles
public: bool overlap (int mem_sel, char ch,char ch1);  		                                // Checking overlap of one vesicle with all other vesicles
public: bool overlap1 (int mem_sel, char ch); 	 			                        // Checking overlap of one vesicle with all other vesicles
public: bool overlap1 (int mem_sel, char ch,char ch1);  		                                // Checking overlap of one vesicle with all other vesicles
public: bool overlap (int mem_sel, char ch,double tempx,double tempy,double tempz);             // Checking overlap of one object with another when displaced by tempx,tempy and tempz
public: bool overlap (int mem_sel, char ch,char ch1,double tempx,double tempy,double tempz);    // Checking overlap of one object with another when displaced by tempx,tempy and tempz
public: bool overlap_vesicle_membrane(int mem_sel);			                        // Checking overlap of one vesicle with all other antigens
public: bool overlap_all (char ch,char ch1);  	                                                // Checking overlap of one vesicle with all other object of type ch1 even out of the linkcell
public: bool overlap_all1 (char ch,char ch1);  	                                                // Checking overlap of one vesicle with all other object of type ch1 even out of the linkcell
public: void overlap_all (int mem1, char ch,int mem2,char ch1);  	                        // Checking overlap of one vesicle with all other object of type ch1 even out of the linkcell
public: bool overlap_all (int mem_sel, char ch,char ch1);  		                        // Checking overlap of one vesicle with all other object of type ch1 even out of the linkcell
public: bool overlap_antigens(int mem1,int mem2);                 	                        // Checking overlap of one vesicle with all other object of type ch1 even out of the linkcell
};
#endif
