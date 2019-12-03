#include "_declarations.h"
#include "_potential.h"
#include "_datainput.h"                                                            // definition for data reading from input.in
#include "_setup.h"                                                                // Lays down the datastructure for system variables
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_linklist.h"                                                             // Construct a linklist connecting space to nano carrier and memb.
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_histogram.h"
#include "_interaction_parameters.h"

void _potential :: init(){
  extern _datainput data;
  extern _histogram hist;
  extern _nanocarrier nc;
  extern _receptor re;
  //@
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();
  this->temperature = data.gettemperature();
  this->kBT = kb*temperature;
  //@
  this->num_nc = nc.getnum_members();
  this->num_antigens = re.getnum_members();
  
  //@
  this->verptr=&__module_datastruct_MOD_ver;                                                  // pointer to vertex type
  this->antptr=&__module_datastruct_MOD_antig;                                                // pointer to antigen type
  this->triptr=&__module_datastruct_MOD_tri;                                                  // pointer to triangle type
}

void _potential :: reset_box_dimensions(){
  extern _datainput data;
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();	
}

double _potential :: dist_receptor_ligand(int ncnum, int abnum, int recep_num){
	extern _nanocarrier nc;
	extern _receptor re;
	//@
	double rrx, rry, rrz;
    rrx = nc.getx_antibody(ncnum, abnum) - re.getxt(recep_num);               // distance between antibody and antigen tip
    rry = nc.gety_antibody(ncnum, abnum) - re.getyt(recep_num);
    rrz = nc.getz_antibody(ncnum, abnum) - re.getzt(recep_num);
	rrx = rrx - this->L*round(rrx/this->L);                                   // implements the PBC along x and y
	rry = rry - this->L*round(rry/this->L);
	return (rrx*rrx + rry*rry + rrz*rrz);                                      // real distance
}

double _potential :: dist_nc_receptor(int ncnum, int recep_num){
	extern _nanocarrier nc;
	extern _receptor re;
	//@
	double rrx, rry, rrz;
    rrx = nc.getxc(ncnum) - re.getxt(recep_num);               // distance between antibody and antigen tip
    rry = nc.getyc(ncnum) - re.getyt(recep_num);
  	rrz = nc.getzc(ncnum) - re.getzt(recep_num);
	rrx = rrx - this->L*round(rrx/this->L);                    // implements the PBC along x and y
	rry = rry - this->L*round(rry/this->L);
	return (rrx*rrx + rry*rry + rrz*rrz);                      // real distance
}

double _potential :: calc (int mem_sel, char ch){
  extern _nanocarrier nc;
  extern _receptor re;
  //@
  double dist,pot;
  int i, *temp, tem;

  pot = 0.0;
  if (ch == 'c' ){
   if(re.getbonded(mem_sel) == 1){                                            // To compute the energy due to bond stretching when antigen is translated
     temp = re.getbondedto(mem_sel);                                          // finds the antibody it is attached to
     dist = this->dist_receptor_ligand(temp[0],temp[1],mem_sel);
     pot = pot + this->reaction(dist,temp[0],temp[1],mem_sel);                // Contribution from the reaction energy
   }
   pot = pot + this->antigen_flexure_energy(mem_sel);                         // Contribution from the antigen flexure energy  
  }

  if (ch == 'v'){                                            // To compute the energy due to bond stretching when vesicle is rotated or translated
    for (i=0; i< nc.getnum_ab(mem_sel); i++){
     if (nc.getbonded(mem_sel, i)){
       tem = nc.getbondedto(mem_sel, i);                     // gives the antigen connected to ith antibody on mem_sel      
       dist = this->dist_receptor_ligand(mem_sel,i,tem);
       pot = pot + this->reaction(dist,mem_sel,i,tem)+ this->antigen_flexure_energy(tem);                          
   }
  }
 }
 return pot;
}

double _potential :: calc (int mem_sel){                                    // compute only the reaction energy for a vesicle
  extern _nanocarrier nc;
  double dist, pot;
  int i,tem;
  pot = 0.0;
  for (i=0; i<nc.getnum_ab(mem_sel); i++){
    if (nc.getbonded(mem_sel, i) == 1){
      tem = nc.getbondedto(mem_sel, i);                                                  // gives the antigen connected to ith antibody on mem_sel     
      dist = this->dist_receptor_ligand(mem_sel,i,tem);
      pot = pot + this->reaction(dist,mem_sel,i,tem);                           
    }    
  }
  return pot;
}

//@ reaction energy
double _potential :: total_energy(){
  double energy = 0.0;
  for (int i=0; i<num_nc; i++){ 
    energy += this->calc(i,'v');
  }
  return energy;
}

double _potential :: total_Renergy(){                                       // to compute the total reaction energy of the system:Jin
  double energy=0.0;
  for (int i=0; i<num_nc; i++){
    energy += this->calc(i);                                          // computes only the reaction energy part 
  }
  return energy;                                                      // reaction energy
}

double _potential :: total_Benergy(){                                       // to compute the total bending energy of the system:Jin
  extern _receptor re;
  double energy=0.0;
  for (int i=0; i<num_antigens; i++){
    if (re.getbonded(i))
      energy += this->antigen_flexure_energy(i);                  // bending energy
    }
  return energy;
}

double _potential :: reaction_energy(int ncnum){
    return this->calc(ncnum);
}

double _potential :: reaction(double dist, int ncnum, int abnum, int antnum, bool kBTconversion){
  extern _intparam intparam;
  extern _nanocarrier nc;
  double energy;
  double *iparam = intparam.get_bellbond_parameters(ncnum, abnum, antnum);              // NC/AB/Receptor dependent interaction parameters 
  energy = nc.get_lambda(ncnum)*(iparam[0] + 0.5*iparam[1]*dist*base_length_sq);        // if bonded a harmonic potential describes An-AB bond
  
  //@ Convert energy to units of kBT if requested
  if (kBTconversion)
    energy = energy/this->kBT;
  return energy;
}


double _potential :: change_reaction(double distold, double distnew, int ncnum, int abnum, int antnum, bool kBTconversion){
  extern _intparam intparam;
  extern _nanocarrier nc;
  double denergy;
  double *iparam = intparam.get_bellbond_parameters(ncnum, abnum, antnum);
  double lambda = nc.get_lambda(ncnum);
  
  if (distnew >= pow(iparam[2],2)){                                                                            
    // @bond breaks in new conformation
    // @This condition will not work since Fortran automatically rejects such a move
    denergy = -1.0*(lambda*(iparam[0] + 0.5*iparam[1]*distold*base_length_sq) + this->antigen_flexure_energy(antnum));
  }
  else{
    //@ Here the change in energy is only due change in the length of the receptor-ligand bonds
    denergy = lambda*(0.5*iparam[1]*(distnew-distold)*base_length_sq);           
  }
  
  //@ Convert energy to units of kBT if requested
  if (kBTconversion)
    denergy = denergy/this->kBT;
  
  return denergy;
  
}

double _potential :: antigen_flexure_energy(int ant_num){
  double dist_flex,kflex,length,theta;
  extern _receptor re;
  kflex = re.get_flexure(ant_num);
  length = re.getlength(ant_num);
  theta = re.get_theta(ant_num);
  dist_flex = length*(sin(theta));                         // compute flexed distance using value of theta stored
  return kflex*(dist_flex*dist_flex)*base_length_sq;                // flex. energy estimated from the deviation perp to the direction of normal
}

//@ biasing potential for bias_mode=H
double _potential :: compute_biasing_potential(int mem_sel, double valold, double valnew){
  extern _nanocarrier nc;
  double del_bias_ener = 0.0;
  double *bp = nc.get_biasparameters(mem_sel);
  del_bias_ener = 0.5*bp[0]*(pow((valnew-bp[1]),2)-pow((valold-bp[1]),2));
  return del_bias_ener;
}

//@ compute the geometric distance between any two position vectors pos1 and pos2
double _potential:: compute_distance(double *pos1, double *pos2){      
  double distance,dx,dy,dz;
  dx = pos1[0]-pos2[0];
  dx -= round(dx/this->L)*this->L;
  dy = pos1[1]-pos2[1];
  dy -= round(dy/this->L)*this->L;
  dz = pos1[2]-pos2[2];
  distance =sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
  return distance;
}

//biasing potential for bias_mode=Z
double _potential :: compute_biasing_potential(int mem_sel, double *rold, double *meanold, double *rnew, double *meannew){    
  extern _nanocarrier nc;
  double del_bias_ener = 0.0;
  double *bp = nc.get_biasparameters(mem_sel);
  del_bias_ener = 0.5*bp[0]*(pow((compute_distance(rnew,meannew)-bp[1]),2) - 
                             pow((compute_distance(rold,meanold)-bp[1]),2));
  return del_bias_ener;
}


