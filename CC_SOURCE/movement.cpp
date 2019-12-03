#include "_declarations.h"
#include "_movement.h"
#include "_datainput.h"                                            // definition for data reading from input.in
#include "_datawriter.h"                                           // definition for data reading from input.in
#include "_nanocarrier.h"                                          // Lays down the datastructure for system variables
#include "_receptor.h"
#include "_membrane.h"
#include "_linklist.h"                                            // Construct a linklist connecting space to nano carrier and memb.
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_interaction_parameters.h"

void _movement :: init(){
  extern _membrane me;
  extern _nanocarrier nc;
  extern _receptor re;
  
  //@
  this->nver = me.getnum_vertices();
  this->pbnum = me.getnum_verticestot();
  this->ntr = me.getnum_triangles();
  this->nlink = me.getnum_links();
  this->num_nc = nc.getnum_members();
  this->num_antigens = re.getnum_members();       

  //@
  this->verptr=&__module_datastruct_MOD_ver;
  this->antptr=&__module_datastruct_MOD_antig;
  this->triptr=&__module_datastruct_MOD_tri;
  this->ncpointer=&__module_datastruct_MOD_nc_f;
}

void _movement :: reset_bond_and_equilibrate(int nsteps, bool DTMCflag){          // freeze NC position and try to generate bonds 
  extern _receptor re;
  extern _nanocarrier nc;
  extern _membrane me;
  extern bool debug_mode;
 
  if (debug_mode){
  cout<<str_red<<"\n\t <RESET BONDS & EQUILIBRATE> "<<str_black;
  cout<<str_black<<"\t"<<nsteps<<" MC Steps "<<str_black;
  }
   
  //@ Break all antigen bonds
  for (int antnum=0 ; antnum < re.getnum_members() ; antnum++){    // over all antigens
    re.setantigentip_position(antnum);
    re.set_theta_phi(antnum,0.0,0.0);                              // reset theta and phi values for the antigen to zero
  }
 
  //@ Break all antibody bonds
  for (int ncnum=0; ncnum < nc.getnum_members(); ncnum++)
    nc.reset_antab_bonds(ncnum);
  
  // focus only on making bonds without translating the NC
  for (int i=0 ; i<nsteps ; i++){                                  // equilibrate the nanocarrier position and antigen positions for nsteps.
	nc.translation();
    nc.rotation();
    re.receptor_diffusion();
    if (((i%5000) == 0) && DTMCflag) me.membrane_montecarlo();
  }
  
  for (int i=0 ; i<nsteps ; i++){                                  // equilibrate the nanocarrier position and antigen positions for nsteps.
	  nc.rotation();
	  nc.bond_formation();
	  re.receptor_diffusion();
	  if (((i%5000) == 0) && DTMCflag) me.membrane_montecarlo();
  }
  
  if (debug_mode) cout<<str_red<<"\t <END OF EQUILIBRATION> \n"<<str_black<<endl;
}

void _movement :: check_antigen_antibody_dist(char *s, int ves){
  extern _datainput data;
  extern _nanocarrier nc;
  extern _receptor re;
  extern _intparam intparam;
  //@
  
  for (int ncnum=0; ncnum<num_nc; ncnum++){
    if (nc.getbonded(ncnum)){
      for(int abnum=0; abnum<nc.getnum_ab(ncnum); abnum++){
        if (nc.getbonded(ncnum,abnum)){
         int antnum = nc.getbondedto(ncnum,abnum);
		 double blen = 1.0;
         double blmax = intparam.get_delrsq(ncnum, abnum, antnum);
         if (blen > blmax){
          cout<<s<<" "<<abnum<<" "<<antnum<<" "<<blen<<" "<<blmax<<endl;
          cin.get();
         }
        }
      }
    }  
  }
}

void _movement :: write_vtkfiles(int confno){
  extern _datawriter w;
  w.VTK_membrane(confno);
  w.VTK_NanoCarrier(confno);
  w.VTK_Antigens(confno);
}
