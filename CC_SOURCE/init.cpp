#include "_declarations.h"
#include "_init.h"
#include "_datainput.h"                                           // definition for data reading from input.in
#include "_nanocarrier.h"                                         // Lays down the datastructure for system variables
#include "_potential.h"                                           // Computes various potential energies
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_selfavoidance.h"
#include "_linklist.h"
#include <algorithm>

//@ Initialization
void _init :: init(int mode){
extern _datainput data;
extern _nanocarrier  nc;
extern bool debug_mode;
if (debug_mode) cout << str_red<<" \n \t <Initializing Nanocarrier Position> "<<str_black<<endl<<endl;
this->num_members = nc.getnum_members();
this->L = data.getperiodic_box_length();
this->H = data.getperiodic_box_height();
if (mode ==1)
  this->random_distribution();
else if (mode==2)
  this->place_at_nearest_location();
if (debug_mode) {
	cout<<str_red<<"\n \t <End of Initialization>"<<str_black<<endl;                             
	cout<<endl<<endmarker<<endl<<endl<<endl;
}
}

void _init :: reset_box_dimensions(){
extern _datainput data;
this->L = data.getperiodic_box_length();
this->H = data.getperiodic_box_height();
}


//@ Randomly distribute the NC's at predefined offset distance above the membrane surface
void _init :: random_distribution(){
extern _linklist link1;
extern _selfavoidance sa;
extern _nanocarrier  nc;
extern MTRand mtrand1;
double xposnc, yposnc, zposnc, meanz, roffset;
int i=0;
do {
    roffset = nc.get_rstart(i);
    int ii = i+1;                                          // for transfer to Fortran
    bool accpflag = false;
    while (!accpflag){
      xposnc = mtrand1.rand(0.1*this->L,0.9*this->L);
      yposnc = mtrand1.rand(0.1*this->L,0.9*this->L);
      meanz = __module_datastruct_MOD_compute_closest_vertex(&ii, &xposnc,&yposnc);
      zposnc = min(this->H - nc.getsoft_radius(i), meanz + roffset);                               // choose a position between membrane and max box height
      do{
        nc.shift_nc_position(i,xposnc,yposnc,zposnc);
        link1.calc('v');
        if (sa.overlap(i,'v','v') || sa.overlap(i,'v','c') || sa.overlap(i,'v','m')){
         xposnc = mtrand1.rand(0.1*this->L,0.9*this->L);
         yposnc = mtrand1.rand(0.1*this->L,0.9*this->L);
         meanz = __module_datastruct_MOD_compute_closest_vertex(&ii, &xposnc,&yposnc);
         zposnc = min(this->H - nc.getsoft_radius(i), meanz + roffset);                               // choose a position between membrane and max box height
        }
        else{
         accpflag=true;
        }
      } while (!accpflag);
      
      i +=1;
    }
    cout<<"NC "<<i<<" -- zposnc "<<zposnc<<endl;
  } while(i<this->num_members);
}

//@ Randomly distribute the NC's at predefined offset distance above the membrane surface
void _init :: place_at_nearest_location(){
extern _linklist link1;
extern _selfavoidance sa;
extern _nanocarrier  nc;
extern _datainput data;
extern bool debug_mode;
extern MTRand mtrand1;
double xposnc, yposnc, zposnc, meanz, roffset;
int i=0;
do {
    roffset = nc.getsoft_radius(i) + data.get_antigen_minlength();
    int ii = i+1;                                                           // for transfer to Fortran
    bool accpflag = false;
	
    while (!accpflag){
      int ncounter = 0;  
      double dshift = 0.0;
      xposnc = mtrand1.rand(0.05*this->L,0.95*this->L);
      yposnc = mtrand1.rand(0.05*this->L,0.95*this->L);
      meanz = __module_datastruct_MOD_compute_closest_vertex(&ii, &xposnc,&yposnc);
      zposnc = min(this->H - nc.getsoft_radius(i), meanz + roffset);                               // choose a position between membrane and max box height
	  // xposnc = yposnc = 500.0;  zposnc = 1000.0;
      do{                                                                                          // Try adjusting the z direction for 10 steps 
        nc.shift_nc_position(i,xposnc,yposnc,zposnc + dshift);
        link1.calc('v');
        if (sa.overlap(i,'v','v') || sa.overlap(i,'v','c') || sa.overlap(i,'v','m')){
          dshift += nc.getradius(i)/100.0;
          ncounter++;
        }   
        else {
         accpflag=true;
         break;
        }
      } while ((ncounter<10) && !accpflag);                                            // If fails then try a new XY position  
    }
    i++;
	if (debug_mode) cout<<"\t \t --> NC "<<i<<"(x,y,z): "<<xposnc<<" "<<yposnc<<" "<<zposnc<<endl;
  } while(i<this->num_members);
}