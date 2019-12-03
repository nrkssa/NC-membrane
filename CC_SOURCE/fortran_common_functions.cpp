/* 5 October 12 - Make sure that the arguments to the functions called from Fortran are pointers since Fortran passes arguments as pointers by default */
#include "_declarations.h"
#include "_datainput.h"       	                                  // definition for data reading from input.in
#include "_linklist.h"
#include "_setup.h"   						  // Lays down the datastructure for system variables
#include "_fortran_structures.h"
#include "_potential.h"
#include "_datawriter.h"
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_membrane.h"

int getcell_error[1]={-111111};

extern "C" {
double antigenflexure_energy_(int  &sel_antigen){
extern _potential p;
return p.antigen_flexure_energy(sel_antigen-1);
}
}

extern "C"{                                                  //Compute the bonded state of the antigen
bool is_antigen_bonded_(int &antigenno){
extern _receptor re;
return re.getbonded(antigenno);
}
}

extern "C"{                                                  //Computes if an antigen-antibody bond is broken
bool does_bond_breaks_(int &antigenno){
extern _receptor re;
return re.does_bond_breaks(antigenno);
}
}

extern "C"{                                                  // Get the x coordinate of the object defined by ch
double get_xcoordinate_(int &objno,char &ch){
extern _setup a;
return a.getxc(objno,ch);
}
}

extern "C"{                                                  // Get the y coordinate of the object defined by ch
double get_ycoordinate_(int &objno,char &ch){
extern _setup a;
return a.getyc(objno,ch);
}
}

extern "C"{                                                 // Get the z coordinate of the object defined by ch
double get_zcoordinate_(int &objno,char &ch){
extern _setup a;
return a.getzc(objno,ch);
}
}


extern "C"{                                                                    //Compute the change in reaction energy due to antigen translation
double antigen_reaction_energy_change_(int &antigenno, double *drnew, double *drold){
int *temp;
double dxold,dyold,dzold,dxnew,dynew,dznew,distold,distnew;
extern _receptor re;
extern _nanocarrier nc;
extern _potential p;
extern _datainput data;
double L = data.getperiodic_box_length();

temp = re.getbondedto(antigenno);
dxold = drold[0] - nc.getx_antibody(temp[0],temp[1]);            // separation between antigen and antibody before translation
dxold -= round(dxold/L)*L;
dyold = drold[1] - nc.gety_antibody(temp[0],temp[1]);   
dyold -= round(dyold/L)*L;
dzold = drold[2] - nc.getz_antibody(temp[0],temp[1]);   
distold = (pow(dxold,2)+pow(dyold,2)+pow(dzold,2));
dxnew = drnew[0] - nc.getx_antibody(temp[0],temp[1]);            // separation between antigen and antibody after translation
dxnew -= round(dxnew/L)*L;
dynew = drnew[1] - nc.gety_antibody(temp[0],temp[1]);   
dynew -= round(dynew/L)*L;
dznew = drnew[2] - nc.getz_antibody(temp[0],temp[1]);   
distnew = (pow(dxnew,2)+pow(dynew,2)+pow(dznew,2));
return p.change_reaction(distold, distnew, temp[0],temp[1],antigenno,true);
}
}

extern "C"{                                                     // Compute the change in reaction energy due to antigen translation
double antigen_reaction_energy_(int &antigenno, double *drnew){
int *temp;
double dxnew,dynew,dznew,distnew;
extern _receptor re;
extern _nanocarrier nc;
extern _potential p;
extern _datainput data;
double L = data.getperiodic_box_length();
temp = re.getbondedto(antigenno);
dxnew = drnew[0] - nc.getx_antibody(temp[0],temp[1]);                                 // separation between antigen and antibody after translation
dxnew -= round(dxnew/L)*L;
dynew = drnew[1] - nc.gety_antibody(temp[0],temp[1]);   
dynew -= round(dynew/L)*L;
dznew = drnew[2] - nc.getz_antibody(temp[0],temp[1]);   
distnew = (pow(dxnew,2)+pow(dynew,2)+pow(dznew,2));
return p.reaction(distnew,temp[0],temp[1],antigenno,true);
}
}

extern "C"{                                                //pass the head linklist for an object
int* get_headlist_(char &ch){
extern _linklist link1;
return link1.gethead(ch);
}
}

extern "C"{                                                //  pass the head linklist for an object
void update_antigenbondedstate_(int &i){
extern _receptor re;
re.updatebonds(i);
}
}


extern "C"{                                                //pass the tail linklist for an object
int* get_lscllist_(char &ch){
extern _linklist link1;
return link1.getlscl(ch);
}
}

extern "C"{
int* getcells_(int &ringsize,int &i,char &ch){
extern _linklist link1;
if (ringsize == 1) return link1.get1ringcells(i-1,ch);
if (ringsize == 2) return link1.get2ringcells(i-1,ch);
cout << "Input error in "<<__func__<<endl;
cout << "Values are "<<ringsize<<" "<<i<<" "<<ch<<endl;
exit(0);
return getcell_error;                                      // This return statement is here just to complete the structure of the function and has no functional role.
}
}

extern "C"{
void print_celldata_(int &obj,char &ch){
extern _linklist link1;
int j,*head,*ring;
head=link1.gethead(ch);
ring=link1.get1ringcells(obj-1,ch);
for(j=0;j<27;j++){
  cout<<"cell number "<<j<<" "<<ring[j]<<" "<<head[ring[j]]<<endl;
}
}
}

extern "C"{
void update_linkcells_(int &vert,double &oldx, double &oldy,double &oldz,char &ch){
extern _linklist link1;
link1.calc(vert-1,oldx,oldy,oldz,ch);
}
}

extern "C"{
void print_cellnumber_(int &vert,char &ch){
extern _linklist link1;
cout << "AAA. link list for chosen member "<<vert-1<<" character "<<ch<<" is "<<link1.getcellnumber(vert-1,ch)<<endl<<endl;
}
}

extern "C"{
}

void write_vtkformat_(int &samp_no,char &ch){
extern _datawriter w;
if (ch=='c') w.VTK_Antigens(samp_no);
if (ch=='v') w.VTK_NanoCarrier(samp_no);
if (ch=='c') w.VTK_membrane(samp_no);
}

// Function required for parsing variable for antigen heterogeneity

extern "C"{
 int get_number_antigens_(){
  extern _datainput data;
  return data.getnum_members('c');
  }
}


extern "C"{
 int get_antigen_types_(){
  extern _datainput data;
  return data.get_types_of_antigens();
  }
}

extern "C"{
 int* get_num_type_antigens_(){
  extern _datainput data;
  return data.get_number_of_antigen_types();
  }
}

extern "C"{
 double* get_radius_antigens_(){
  extern _datainput data;
  return data.get_antigen_radius();
  }
}

extern "C"{
 double* get_length_antigens_(){
  extern _datainput data;
  return data.get_antigen_length();
  }
}

extern "C"{
 double* get_flexure_antigens_(){
  extern _datainput data;
  return data.get_antigen_flexure();
  }
}
