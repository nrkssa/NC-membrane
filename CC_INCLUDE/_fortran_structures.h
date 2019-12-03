// Replica of the fortran data structure defined in __module_datastruct.f
// Make sure every entry in the struct given below match variable by variable to those defined in Fortran type
// If some entries are missing or interchanged the code would start to generate wrong results

#ifndef __FORTRAN_STRUCTURES_H__
#define __FORTRAN_STRUCTURES_H__
#define maxnc 100

extern "C" int __module_datastruct_MOD_verno;
extern "C" int __module_datastruct_MOD_maxnc;
extern "C" int __module_datastruct_MOD_nver;
extern "C" int __module_datastruct_MOD_tlink;
extern "C" int __module_datastruct_MOD_pbnum;
extern "C" int __module_datastruct_MOD_ntr;
extern "C" int __module_datastruct_MOD_gsize;
extern "C" int __module_datastruct_MOD_num_antigens;
extern "C" int __module_datastruct_MOD_num_nanocarrier;
extern "C" int __module_datastruct_MOD_num_antibodies;
extern "C" double __module_datastruct_MOD_vesicle_radius;
extern "C" double __module_datastruct_MOD_vesicle_soft_radius;
extern "C" double __module_datastruct_MOD_antigen_radius;
extern "C" double __module_datastruct_MOD_eqm_bond_dist;
extern "C" double __module_datastruct_MOD_max_bond_stretch;
extern "C" double __module_datastruct_MOD_periodic_box_length;
extern "C" double __module_datastruct_MOD_periodic_box_height;
extern "C" double __module_datastruct_MOD_blen;
extern "C" int __module_datastruct_MOD_nlinkcells;
extern "C" double __module_datastruct_MOD_membrane_init_z;
extern "C" double __module_datastruct_MOD_depth;
extern "C" double __module_datastruct_MOD_period;
extern "C" char __module_datastruct_MOD_memb_geom[100];
extern "C" bool __module_datastruct_MOD_memb_restart_flag;
extern "C" int __module_datastruct_MOD_maxantigen;
extern "C" int __module_datastruct_MOD_curr_step;
extern "C" double __module_datastruct_MOD_kz;                         // kbias in units of k_BT
extern "C" double __module_datastruct_MOD_kh;                         // kbias in units of k_BT
extern "C" double __module_datastruct_MOD_rref;                       // current value of zref  
extern "C" double __module_datastruct_MOD_href;                       // current value of zref  
extern "C" bool __module_datastruct_MOD_debug_mode;
extern "C" int __module_datastruct_MOD_processor_number;
extern "C" double __module_datastruct_MOD_beta;


struct vertex {                                                      // C++ struct corresponding to TYPE vertex
 double vnor[3][1],vcoord[3][1],t1[3][1],t2[3][1];               
 int vneipt[10],vneitr[10],PBCver[10];             
 double L2G[3][3],HHM[3][3];                                         // Local to Global matrix
 int nonei,boundary,pbimno,pbmap[3],imver,nover,boxvert;    
 double mcur,cur1,cur2,czero,totarea;                     
 int antigen_flag,cellno,neigh;                                            // Antigen connected to a vertex
 int antigen_list[10],nantigen,shadownc[100],czero_flag;
 }; 
extern "C" struct vertex __module_datastruct_MOD_ver;

struct triangle {                                                  
double ar;                                                            // c++ strcut corresponding to TYPE triangle 
int pbflag,boxtriangle ;    
int nantigen;
int li[3],vert[3];    
double fnor[3][1];
double HHM[3][3];                                                     // householder matrix for the triangle based on its face normal
int antigen_list[30];                                                 // A triangle can accomodate 3*maxantigen (thrice the max. number of antigens on a vertex)
double vol;
};
extern "C" struct triangle __module_datastruct_MOD_tri;

struct antigen{                                                        // c++ struct for TYPE antigen  
int vertex;
double base_coord[3][1],tip_coord[3][1],theta,phi;
double disp1,disp2,length,radius,kflex;
int diffus_tri,ant_type;
};
extern "C" struct antigen __module_datastruct_MOD_antig;


struct nanocarrier{
int nshadow_vert,shadow_vertices[10000];     
double kbias, biasref, lambda, shadow_cutoff, shadow_cutoffsq;
double radius, soft_radius, multivalency;
double coord[3][1], meanr[3][1], meancurv, meandist;
char bias_mode;
int pbcx[10000], pbcy[10000];
};
extern "C" struct nanocarrier __module_datastruct_MOD_nc_f;
#endif

