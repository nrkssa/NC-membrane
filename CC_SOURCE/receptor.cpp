#include "_declarations.h"
#include "_constant.h"
#include "_receptor.h"
#include "_nanocarrier.h"
#include "_membrane.h"
#include "_potential.h"
#include "_datainput.h"                       // definition for data reading from input.in
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_linklist.h"
#include "_selfavoidance.h"
#include "_datawriter.h"
#include "_debugger.h"
#include "_interaction_parameters.h"

void _RECEPTOR :: init(int i) {
  extern _datainput data;
  this->index= i;
  this->bonded = false;
  this->bondedtoves = -1;
  this->bondedtoab = -1;
  this->L = data.getperiodic_box_length();
  this->beta = 1.0/(kb*data.gettemperature());
  //@ Links to the fortran pointer for antigen, vertex and triangle
  this->antptr=&__module_datastruct_MOD_antig;
  this->verptr = &__module_datastruct_MOD_ver;
  this->triptr = &__module_datastruct_MOD_tri;
  this->initialize_storage();                        // Store time series of positions
}

void _RECEPTOR :: init(int i1, double* fv, int* i2, double* fv1, int* i3){
	extern _datainput data;
	this->antptr=&__module_datastruct_MOD_antig;
	this->verptr = &__module_datastruct_MOD_ver;
	
	this->index = i1;
	this->L = data.getperiodic_box_length();
	this->beta = 1.0/(kb*data.gettemperature());
	
	(this->antptr+this->index)->base_coord[0][0] = fv[0];
	(this->antptr+this->index)->base_coord[1][0] = fv[1];
	(this->antptr+this->index)->base_coord[2][0] = fv[2];
	(this->antptr+this->index)->tip_coord[0][0] = fv[3];
	(this->antptr+this->index)->tip_coord[1][0] = fv[4];
	(this->antptr+this->index)->tip_coord[2][0] = fv[5];
	(this->antptr+this->index)->theta = fv[6];
	(this->antptr+this->index)->phi = fv[7];
	(this->antptr+this->index)->disp1 = fv[8];
	(this->antptr+this->index)->disp2 = fv[9];
	(this->antptr+this->index)->vertex = i2[0];
	(this->antptr+this->index)->diffus_tri = i2[1];
	(this->antptr+this->index)->ant_type = i2[2];
	(this->antptr+this->index)->radius = fv1[0];
	(this->antptr+this->index)->length = fv1[1];
	(this->antptr+this->index)->kflex = fv1[2];
	if (i3[0]==0)
		this->bonded = false;
	else
		this->bonded = true;
	this->bondedtoves = i3[1];
	this->bondedtoab = i3[2];
	this->type = i3[3];
}

void _RECEPTOR :: reset_box_dimensions(){
  extern _datainput data;
  this->L = data.getperiodic_box_length();	
}

void _RECEPTOR :: dump_conf(ofstream &ofile){
	ofile<<this->index<<" "<<(this->antptr+this->index)->base_coord[0][0]<<" "<<
	       (this->antptr+this->index)->base_coord[1][0]<<" "<<
	       (this->antptr+this->index)->base_coord[2][0]<<" "<<
	       (this->antptr+this->index)->tip_coord[0][0]<<" "<<
	       (this->antptr+this->index)->tip_coord[1][0]<<" "<<
	       (this->antptr+this->index)->tip_coord[2][0]<<" "<<
	       (this->antptr+this->index)->theta<<" "<<
	       (this->antptr+this->index)->phi<<" "<<
	       (this->antptr+this->index)->disp1<<" "<<
	       (this->antptr+this->index)->disp2<<" "<<
	       (this->antptr+this->index)->vertex<<" "<<
	       (this->antptr+this->index)->diffus_tri<<" "<<
	       (this->antptr+this->index)->ant_type<<" "<<
	       (this->antptr+this->index)->radius<<" "<<
	       (this->antptr+this->index)->length<<" "<<
	       (this->antptr+this->index)->kflex<<" "<<
	       this->bonded<<" "<<
	       this->bondedtoves<<" "<<
	       this->bondedtoab<<" "<<this->type<<endl;
}

void _RECEPTOR :: initialize_storage(){
	extern _datainput data;
	this->ntsamples = 10000;
	this->ctsamples = -1;
	this->trbx = new(nothrow) float[this->ntsamples];
	this->trby = new(nothrow) float[this->ntsamples];
	this->trbz = new(nothrow) float[this->ntsamples];
	this->trtx = new(nothrow) float[this->ntsamples];
	this->trty = new(nothrow) float[this->ntsamples];
	this->trtz = new(nothrow) float[this->ntsamples];
	this->tangt = new(nothrow) float[this->ntsamples];
	this->tangp = new(nothrow) float[this->ntsamples];	
	this->tbond = new(nothrow) int[this->ntsamples];
	this->tblig = new(nothrow) int[this->ntsamples];
	this->tbnc = new(nothrow)  int[this->ntsamples];
	
	//@ often times, the values of ntsamples may be large for the memory of your machine
	if (!this->trbx || !this->trby || !this->trbz ||
	    !this->trtx || !this->trty || !this->trtz ||
	    !this->tangt || !this->tangp ||
	    !this->tbond || !this->tblig || !this->tbnc){
			cout<<str_red<<"Failure in memory allocation "<<__PRETTY_FUNCTION__<<endl;
			cout<<"exiting: change the values of ntsamples and proceed" <<endl;
			exit(1);
	    }
}


void _RECEPTOR :: initialize_storage_variables(){
	fill_n(this->trbx,this->ntsamples,0.0);
	fill_n(this->trby,this->ntsamples,0.0);
	fill_n(this->trbz,this->ntsamples,0.0);
	fill_n(this->trtx,this->ntsamples,0.0);
	fill_n(this->trty,this->ntsamples,0.0);
	fill_n(this->trtz,this->ntsamples,0.0);
	fill_n(this->tangt,this->ntsamples,0.0);
	fill_n(this->tangp,this->ntsamples,0.0);
	fill_n(this->tbond,this->ntsamples,-1);
	fill_n(this->tblig,this->ntsamples,-1);
	fill_n(this->tbnc,this->ntsamples,-1);
}

void _RECEPTOR :: store_trajectory(int window){
	this->ctsamples += 1; 
	this->trbx[ctsamples] =  (float) this->getxc();
	this->trby[ctsamples] =  (float) this->getyc();
	this->trbz[ctsamples] =  (float) this->getzc();
	this->trtx[ctsamples] =  (float) this->getxt();
	this->trty[ctsamples] =  (float) this->getyt();
	this->trtz[ctsamples] =  (float) this->getzt();
	this->tangt[ctsamples] = (float)  this->get_theta();
	this->tangp[ctsamples] = (float)  this->get_phi();
	if (this->getbonded()){
		this->tbond[ctsamples] = 1;
		this->tblig[ctsamples] = this->getbondedtoab();
		this->tbnc[ctsamples] = this->getbondedtoves();
	}
	else{
		this->tbond[ctsamples] = 0;
	}
	
	//@ we can only store a certain piece of the trajectory to the memory. 
	if (this->ctsamples == this->ntsamples-1){
		this->dump_trajectory(window);
		this->initialize_storage_variables();
		this->ctsamples=-1;
	}
}

void _RECEPTOR :: dump_trajectory(int window){
	stringstream fstream;
	fstream<<"../SYS_STATE/DUMP/receptor-"<<this->index<<"-ENS-"<<window<<".dat";
	ofstream dumpfile;
	dumpfile.open(fstream.str().c_str(), ios::out | ios::app | ios::binary);
	for (int i=0; i<this->ctsamples; i++){
		dumpfile<<trbx[i]<<" "<<trby[i]<<" "<<trbz[i]<<" "
				<<trtx[i]<<" "<<trty[i]<<" "<<trtz[i]<<" "
				<<tangt[i]<<" "<<tangt[i]<<" "
				<<tbond[i]<<" "<<tblig[i]<<" "<<tbnc[i]<<endl;
	}
	dumpfile.close();
	fstream.str("");
}


void _RECEPTOR :: setradius(double radius1) {
  (this->antptr+this->index)->radius = radius1;
}

void _RECEPTOR :: setxc(double xc1) {
  (this->antptr+index)->base_coord[0][0] = xc1;
}

void _RECEPTOR :: setyc(double yc1) {
  (this->antptr+this->index)->base_coord[1][0] = yc1;
}

void _RECEPTOR :: setzc(double zc1) {
  (this->antptr+this->index)->base_coord[2][0] = zc1;
}

double _RECEPTOR :: getradius() {
  return (this->antptr+this->index)->radius;
}

double _RECEPTOR :: getlength(){
  return (this->antptr+this->index)->length;
}

int _RECEPTOR :: gettype(){
  return (this->antptr+this->index)->ant_type;
}

double _RECEPTOR :: get_flexure(){
  return (this->antptr+this->index)->kflex;
}

double _RECEPTOR :: getxc() {
  return (this->antptr+this->index)->base_coord[0][0];
}

double _RECEPTOR :: getyc() {
  return (this->antptr+this->index)->base_coord[1][0];
}

double _RECEPTOR :: getzc() {
  return (this->antptr+this->index)->base_coord[2][0];
}

double _RECEPTOR :: getxt() {
  return (this->antptr+this->index)->tip_coord[0][0];
}

double _RECEPTOR :: getyt() {
  return (this->antptr+this->index)->tip_coord[1][0];
}

double _RECEPTOR :: getzt() {
  return (this->antptr+this->index)->tip_coord[2][0];
}

void _RECEPTOR :: setxt(double xt1) {          
  (this->antptr+this->index)->tip_coord[0][0] = xt1;
}
  
void _RECEPTOR :: setyt(double yt1) {
  (this->antptr+this->index)->tip_coord[1][0] = yt1;
}
 
void _RECEPTOR :: setzt(double zt1) {
  (this->antptr+this->index)->tip_coord[2][0] = zt1;
}

void _RECEPTOR :: set_theta(double theta){
  (this->antptr+this->index)->theta=theta;
}

void _RECEPTOR :: set_phi(double phi){
  (this->antptr+this->index)->phi=phi;
}

double _RECEPTOR :: get_theta(){
  return (this->antptr+this->index)->theta;
}

double _RECEPTOR :: get_phi(){
  return (this->antptr+this->index)->phi;
}

void _RECEPTOR :: setantigentip_position() {
  double normal_dir[3][1],rrx,rry,rrz;
  int tri_num = (this->antptr+this->index)->diffus_tri-1;                                           // get the triangle number
  if(tri_num != -1){
    normal_dir[0][0]=(this->triptr+tri_num)->fnor[0][0];
    normal_dir[1][0]=(this->triptr+tri_num)->fnor[1][0];                                            // choose the normal that is currently active
    normal_dir[2][0]=(this->triptr+tri_num)->fnor[2][0];
  }
  else{
    int ver_num = (this->antptr+this->index)->vertex-1;
    normal_dir[0][0]=(this->verptr+ver_num)->vnor[0][0];
    normal_dir[1][0]=(this->verptr+ver_num)->vnor[1][0];
    normal_dir[2][0]=(this->verptr+ver_num)->vnor[2][0];
  }
  //cout<<(this->antptr)<<" "<<(this->antptr+index)<<endl;
  //cout<<(this->antptr+index)->base_coord[0][0]<<" "<<(this->antptr+index)->base_coord[1][0]<<" "<<(this->antptr+index)->base_coord[2][0]<<endl;
  rrx = (this->antptr+index)->base_coord[0][0]+normal_dir[0][0]*this->ant_length;
  if(rrx > this->L) rrx -= this->L;
  if(rrx < 0) rrx += this->L;
  rry = (this->antptr+index)->base_coord[1][0]+normal_dir[1][0]*this->ant_length;
  if(rry > this->L) rry -= this->L;
  if(rry < 0) rry += this->L;
  rrz = (this->antptr+index)->base_coord[2][0]+normal_dir[2][0]*this->ant_length;
  (this->antptr+index)->tip_coord[0][0] = rrx;
  (this->antptr+index)->tip_coord[1][0] = rry;
  (this->antptr+index)->tip_coord[2][0] = rrz;
}

void _RECEPTOR :: setL(double L1) {
  this->L = L1;
}

bool _RECEPTOR :: getbonded(){
  return this->bonded;
}

void _RECEPTOR :: setbonded(bool temp){
  this->bonded = temp;
}

int _RECEPTOR :: getbondedtoves() {
  return this->bondedtoves;
}

void _RECEPTOR :: setbondedtoves(int temp) {
  this->bondedtoves = temp;
}

int _RECEPTOR :: getbondedtoab() {
  return this->bondedtoab;
}

void _RECEPTOR :: setbondedtoab(int temp) {
  this->bondedtoab = temp;
}

void _RECEPTOR :: antigen_translation () {                                                     // move an antigen to one of its nearest neighbour : Ram
  extern MTRand mtrand1;
  extern _membrane me;
  extern _potential p;
  extern _linklist link1;
  extern _selfavoidance sa;
  double In_energy=0.0,Fi_energy=0.0,delE=0.0;
  double baseold[3][1],tipold[3][1],trvec1[3][1],trvec2[3][1],tri_disp[3][1],normal_dir[3][1];
  double hhmat[3][3];
  int sel_nei,vertex_index,new_vertex;                                                 // vertex_index --> vertex number associated with antigen number index
  int tri_list[3],tri_list1[3],curr_tri;
  int i=0,n,boundary_tri_flag,tri_ant_index=-1,ver_ant_index=-1;
  bool move_accp_flag;
  double disp1,disp2;

  for(i=0;i<3;i++){
    hhmat[i][0]=hhmat[i][1]=hhmat[i][2]=0;
    baseold[i][0]=(this->antptr+index)->base_coord[i][0];                             // store the old base and tip coordinates
    tipold[i][0]=(this->antptr+index)->tip_coord[i][0];
  }

  boundary_tri_flag = 0;
  In_energy = p.calc(index,'c');                                                        // compute the energy before translation
  curr_tri = (this->antptr+index)->diffus_tri-1;                                              // the triangle the antigen is associated with
  vertex_index = (this->antptr+index)->vertex-1;                                              // vertex connected to the antigen
  
//->
  if (curr_tri != -1){                                                               // If the antigen is already diffusing on a triangle face 
    n=0;
    tri_list1[n]=vertex_index;					           // find the free vertices in curr_tri that includes new_vertex-1
    for (i=0;i<3;i++){
      if (((this->triptr+curr_tri)->vert[i]-1) != vertex_index)
	tri_list1[++n]=(this->triptr+curr_tri)->vert[i]-1;
      }
      sel_nei = mtrand1.randInt(2);                                                // randomly select a vertex within the associated triangle
      new_vertex = tri_list1[sel_nei];                   // the new vertex is chosen from the list

      if(new_vertex >= me.getnum_vertices()){
	new_vertex=(this->verptr+new_vertex)->imver-1;   // if the chosen vertex is an image, the antigen is moved to its parent vertex
	boundary_tri_flag=1;                             // boundary_tri_flag=1 will ensure that the antigen associated triangle is set to -1 on successful move
      }

      if(boundary_tri_flag==1){
	tri_disp[0][0]=tri_disp[1][0]=tri_disp[2][0]=0.0;         // The antigen is transfered to the image vertex only (Things are complicated if otherwise)
	normal_dir[0][0]=(this->verptr+new_vertex)->vnor[0][0];
	normal_dir[1][0]=(this->verptr+new_vertex)->vnor[1][0];
	normal_dir[2][0]=(this->verptr+new_vertex)->vnor[2][0];
      }
      else
      {
	n=0;
	tri_list[n] = new_vertex;                                             // find the free vertices in curr_tri that includes new_vertex-1
	for (i=0;i<3;i++){
	  if (tri_list1[i] != new_vertex){
	  tri_list[++n]=tri_list1[i];
	  }
	}
	disp1=(this->antptr+index)->disp1;
	disp2=(this->antptr+index)->disp2;
	// compute the edge vectors from new_vertex to its neighbors
	trvec1[0][0]=(this->verptr+tri_list[1])->vcoord[0][0]-(this->verptr+tri_list[0])->vcoord[0][0];  
	trvec1[1][0]=(this->verptr+tri_list[1])->vcoord[1][0]-(this->verptr+tri_list[0])->vcoord[1][0];
	trvec1[2][0]=(this->verptr+tri_list[1])->vcoord[2][0]-(this->verptr+tri_list[0])->vcoord[2][0];
	trvec2[0][0]=(this->verptr+tri_list[2])->vcoord[0][0]-(this->verptr+tri_list[0])->vcoord[0][0];
	trvec2[1][0]=(this->verptr+tri_list[2])->vcoord[1][0]-(this->verptr+tri_list[0])->vcoord[1][0];
	trvec2[2][0]=(this->verptr+tri_list[2])->vcoord[2][0]-(this->verptr+tri_list[0])->vcoord[2][0];
	tri_disp[0][0]=disp1*trvec1[0][0]+disp2*trvec2[0][0];                                         // displacement vector along the triangle curr_tri
	tri_disp[1][0]=disp1*trvec1[1][0]+disp2*trvec2[1][0];
	tri_disp[2][0]=disp1*trvec1[2][0]+disp2*trvec2[2][0];
	normal_dir[0][0]=(this->triptr+curr_tri)->fnor[0][0];    //triangle face normal is used in place of the vertex normal vector  to determine the tip pos 
	normal_dir[1][0]=(this->triptr+curr_tri)->fnor[1][0];
	normal_dir[2][0]=(this->triptr+curr_tri)->fnor[2][0];
      }
  }

  else{                                                                            // If the antigen is not associated with  a triangle
    sel_nei=mtrand1.randInt((this->verptr+vertex_index)->nonei-1);          // Choose a neighbouring vertex around the current antigen position
    new_vertex=(this->verptr+vertex_index)->vneipt[sel_nei]-1;
    if (new_vertex >= me.getnum_vertices()){
      new_vertex=(this->verptr+new_vertex)->imver-1;                    // if the chosen vertex is a phantom substitute it with its image 
      boundary_tri_flag=1;
    }

    if((this->verptr+new_vertex)->antigen_flag == 1) {                                    // do not allow for overlap of antigens
      new_vertex=vertex_index;
    }

    else {
      tri_disp[0][0]=tri_disp[1][0]=tri_disp[2][0]=0.0;                    // Initialize the displacement vector along the associated triangle
      normal_dir[0][0]=(this->verptr+new_vertex)->vnor[0][0]; 
      normal_dir[1][0]=(this->verptr+new_vertex)->vnor[1][0];
      normal_dir[2][0]=(this->verptr+new_vertex)->vnor[2][0];
    }
  }

  if ((new_vertex != vertex_index) && ((this->verptr+new_vertex)->nantigen < maxantigen_ver)){  
    (this->antptr+index)->base_coord[0][0]=(this->verptr+new_vertex)->vcoord[0][0]+0.5*tri_disp[0][0];   // new base coordinate
    (this->antptr+index)->base_coord[1][0]=(this->verptr+new_vertex)->vcoord[1][0]+0.5*tri_disp[1][0];                    
    (this->antptr+index)->base_coord[2][0]=(this->verptr+new_vertex)->vcoord[2][0]+0.5*tri_disp[2][0]; 
    // base coordinates do not obey PBC since they move on the mesh    
	
    if (this->bonded==1){                                                                 // If the antigen is bonded, restore the existing flexure
      double antigx_zframe = this->ant_length*sin((antptr+index)->theta)*cos((antptr+index)->phi);
      double antigy_zframe = this->ant_length*sin((antptr+index)->theta)*sin((antptr+index)->phi);
      double antigz_zframe = this->ant_length*cos((antptr+index)->theta);
      if ((curr_tri != -1) && (boundary_tri_flag != 1)){
	int triang_no = curr_tri+1;
	char ch='f';
	__module_curvcalc_MOD_compute_householdermatrix(&triang_no,&ch);
	hhmat[0][0] = (this->triptr+curr_tri)->HHM[0][0];
	hhmat[0][1] = (this->triptr+curr_tri)->HHM[0][1];
	hhmat[0][2] = (this->triptr+curr_tri)->HHM[0][2];
	hhmat[1][0] = (this->triptr+curr_tri)->HHM[1][0];
	hhmat[1][1] = (this->triptr+curr_tri)->HHM[1][1];
	hhmat[1][2] = (this->triptr+curr_tri)->HHM[1][2];
	hhmat[2][0] = (this->triptr+curr_tri)->HHM[2][0];
	hhmat[2][1] = (this->triptr+curr_tri)->HHM[2][1];
	hhmat[2][2] = (this->triptr+curr_tri)->HHM[2][2];
      }
      else{
	int vertex_num = new_vertex+1;
	char ch = 'n';
	__module_curvcalc_MOD_compute_householdermatrix(&vertex_num,&ch);
	hhmat[0][0] = (this->verptr+new_vertex)->HHM[0][0];
	hhmat[0][1] = (this->verptr+new_vertex)->HHM[0][1];
	hhmat[0][2] = (this->verptr+new_vertex)->HHM[0][2];
	hhmat[1][0] = (this->verptr+new_vertex)->HHM[1][0];
	hhmat[1][1] = (this->verptr+new_vertex)->HHM[1][1];
	hhmat[1][2] = (this->verptr+new_vertex)->HHM[1][2];
	hhmat[2][0] = (this->verptr+new_vertex)->HHM[2][0];
	hhmat[2][1] = (this->verptr+new_vertex)->HHM[2][1];
	hhmat[2][2] = (this->verptr+new_vertex)->HHM[2][2];
      }
	
      (this->antptr+index)->tip_coord[0][0] = hhmat[0][0]*antigx_zframe+hhmat[0][1]*antigy_zframe+hhmat[0][2]*antigz_zframe;    
      // transform the new antigen orientation along z to the normal using Householder trans
      (this->antptr+index)->tip_coord[1][0] = hhmat[1][0]*antigx_zframe+hhmat[1][1]*antigy_zframe+hhmat[1][2]*antigz_zframe; 
      (this->antptr+index)->tip_coord[2][0] = hhmat[2][0]*antigx_zframe+hhmat[2][1]*antigy_zframe+hhmat[2][2]*antigz_zframe; 
    }
    else{
      // Preserve the antigen to be along the normal direction even when moved
      (this->antptr+index)->tip_coord[0][0]=(this->antptr+index)->base_coord[0][0]+this->ant_length*normal_dir[0][0];                 
      (this->antptr+index)->tip_coord[1][0]=(this->antptr+index)->base_coord[1][0]+this->ant_length*normal_dir[1][0];     
      (this->antptr+index)->tip_coord[2][0]=(this->antptr+index)->base_coord[2][0]+this->ant_length*normal_dir[2][0];     
    }

    if ((this->antptr+index)->tip_coord[0][0] < 0.0) (this->antptr+index)->tip_coord[0][0] += this->L;          // periodic boundary along X
    if ((this->antptr+index)->tip_coord[0][0] > L) (this->antptr+index)->tip_coord[0][0] -= this->L;
    if ((this->antptr+index)->tip_coord[1][0] < 0.0) (this->antptr+index)->tip_coord[1][0] += this->L;         // periodic boundary along Y 
    if ((this->antptr+index)->tip_coord[1][0] > L)   (this->antptr+index)->tip_coord[1][0] -= this->L; 
    move_accp_flag=true;                                                                                          // move is accepted till this point

    if ((sa.overlap(index,'c','c')==true ) || (sa.overlap(index,'c','v')==true)) move_accp_flag=false;       // reject the move if it violates self-avoidance
    if ((move_accp_flag==true) && (bonded==1) && (does_bond_breaks()==true)) move_accp_flag=false;     // move is rejected is bond breakage is detected

    if(move_accp_flag == true){
      Fi_energy = p.calc(index,'c');                                              // compute the energy after translation
      delE = Fi_energy-In_energy;                                                                // Energy change after update
      if ((delE > 0) && (mtrand1.rand() > exp(-this->beta*delE))){                               // Metropolis scheme to accept the move 
	move_accp_flag=false;
      }
      else                                                                           // move is accepted and book keeping is performed in what follows 
      {
	(this->antptr+index)->vertex=new_vertex+1;                                       // the antigen is associated with a new vertex  
	for (i=0;i<(this->verptr+vertex_index)->nantigen;i++){
	  if ((this->verptr+vertex_index)->antigen_list[i] == index+1 )
	    ver_ant_index=i;
	}
	for (i=ver_ant_index;i<(this->verptr+vertex_index)->nantigen-1;i++){
	  (this->verptr+vertex_index)->antigen_list[i]=(verptr+vertex_index)->antigen_list[i+1];
	  (this->verptr+vertex_index)->antigen_list[i+1]=0;
	}
	(this->verptr+vertex_index)->nantigen -= 1;                      // Remove one antigen from the antigen count of the old vertex
	if ((this->verptr+vertex_index)->nantigen == 0)
	  (this->verptr+vertex_index)->antigen_flag=0;                   // set antigen flag for old vertex to zero
	
	(this->verptr+new_vertex)->antigen_list[(this->verptr+new_vertex)->nantigen] = index+1;    // Connect the selected antigen to the new vertex
	(this->verptr+new_vertex)->nantigen += 1;                                                // increment antigen count for new vertex by 1
	(this->verptr+new_vertex)->antigen_flag=1;                                               // set antigen flag for new vertex to one

	if(boundary_tri_flag == 1){
	  (this->antptr+index)->diffus_tri=0;                                                      // the antigen loses its diffusion triangle
	  if (curr_tri != -1){
	    for (i=0;i<(this->triptr+curr_tri)->nantigen;i++){                                                 
		if ((this->triptr+curr_tri)->antigen_list[i] == index+1) 
		  tri_ant_index=i;                                       // find the current position of antigen index in the list of curr_tri
	    }
	    for (i=tri_ant_index;i<(this->triptr+curr_tri)->nantigen-1;i++){
	      (this->triptr+curr_tri)->antigen_list[i]=(this->triptr+curr_tri)->antigen_list[i+1];  // transfer all contents from  R->L
	      (this->triptr+curr_tri)->antigen_list[i+1]=0;  
	    }
	    (this->triptr+curr_tri)->nantigen -= 1;   // reduce the number of antigen on the triangle by 1 since the triangle also loses an antigen
	  }
	}
	link1.calc(index,baseold[0][0],baseold[1][0],baseold[2][0],'c');                 // update the linkcell corresponding to the new position
      }
    }
	 
    if (move_accp_flag == false){                                                        // on rejection of the move restore the old coordinates
      (this->antptr+index)->base_coord[0][0]=baseold[0][0];
      (this->antptr+index)->base_coord[1][0]=baseold[1][0];
      (this->antptr+index)->base_coord[2][0]=baseold[2][0];                              // Restore the old coordinates
      (this->antptr+index)->tip_coord[0][0]=tipold[0][0];
      (this->antptr+index)->tip_coord[1][0]=tipold[1][0];
      (this->antptr+index)->tip_coord[2][0]=tipold[2][0];
    }
  } 
}

bool _RECEPTOR :: does_bond_breaks(double tempx,double tempy,double tempz){
  bool bond_break = false;
  extern _nanocarrier nc;
  extern _intparam intparam;
  if (this->bonded){
    double dx = this->getxt() + tempx - nc.getx_antibody(this->bondedtoves,this->bondedtoab);             // xcomp of new position
    double dy = this->getyt() + tempy - nc.gety_antibody(this->bondedtoves,this->bondedtoab);             // ycomp of new position
    double dz = this->getzt() + tempz - nc.getz_antibody(this->bondedtoves,this->bondedtoab);             // zcomp of new position
    dx = dx-round(dx/this->L) * this->L; dy = dy-round(dy/this->L)*this->L;
    double dist = dx*dx + dy*dy + dz*dz;
    double max_AB_ANT_distance = intparam.get_delrsq(this->bondedtoves, this->bondedtoab, this->index);
    if (dist > max_AB_ANT_distance)
      bond_break=true;                                               // if the movement stretches the bond beyond the reaction distance reject the move
  }
  return bond_break;
}

bool _RECEPTOR :: does_bond_breaks(){
  bool bond_break = false;
  extern _nanocarrier nc;
  extern _intparam intparam;
  if (this->bonded){
    double dx = this->getxt() - nc.getx_antibody(this->bondedtoves,this->bondedtoab);        // xcomp of new position
    double dy = this->getyt() - nc.gety_antibody(this->bondedtoves,this->bondedtoab);        // ycomp of new position
    double dz = this->getzt() - nc.getz_antibody(this->bondedtoves,this->bondedtoab);        // zcomp of new position
    dx = dx - round(dx/this->L)*this->L;
    dy = dy - round(dy/this->L)*this->L;
    double dist = dx*dx + dy*dy + dz*dz;
    double max_AB_ANT_distance = intparam.get_delrsq(this->bondedtoves, this->bondedtoab, this->index);
    if (dist > max_AB_ANT_distance)
      bond_break = true;                             	  // if the movement stretches the bond beyond the reaction distance reject the move
  }
  return bond_break;
}

void _receptor :: init() {
  this->set_array_sizes();
  this->c = new (nothrow) _RECEPTOR[this->num_antigens];
  if (this->c ==0) cout<<"\nError in System: memory could not be allocated";
  for(int i=0; i<this->num_antigens; i++)  
    this->c[i].init(i);
}

void _receptor :: init(ifstream &ifile){
	this->set_array_sizes();
	this->c = new (nothrow) _RECEPTOR[this->num_antigens];
	string temps;
	int i1,i2[3],i3[4];
	double fv[10],fv1[3];
	while (getline(ifile,temps)){
		if (temps.compare("<ANTIGEN>") == 0)
			break;
	}
  
	int n=0;
	while (n < this->num_antigens){						
		ifile>>i1>>fv[0]>>fv[1]>>fv[2]>>fv[3]>>fv[4]>>fv[5]>>
				   fv[6]>>fv[7]>>fv[8]>>fv[9]>>i2[0]>>i2[1]>>
				   i2[2]>>fv1[0]>>fv1[1]>>fv1[2]>>i3[0]>>i3[1]>>i3[2]>>i3[3];
		this->c[i1].init(i1,fv,i2,fv1,i3);
		n++;
	}
}

void _receptor :: set_array_sizes(){
  extern _datainput data;
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();
  this->num_antigens = data.getnum_members('c');        // number of antigens
  this->type_ant = data.get_types_of_antigens();
  this->num_type_ant = new(nothrow) int[this->type_ant];
  this->radius_ant = new(nothrow) double[this->type_ant];
  this->length_ant = new(nothrow) double[this->type_ant];
  this->kflex = new(nothrow) double[this->type_ant];	
  this->num_type_ant = data.get_number_of_antigen_types();
  this->radius_ant = data.get_antigen_radius();
  this->length_ant = data.get_antigen_length();
  this->kflex = data.get_antigen_flexure();
  this->ant_pattern = data.get_antigen_pattern();
}

_receptor :: _receptor(double L1, double H1, int num_antigens1){
  this->L = L1;
  this->H = H1;
  this->num_antigens = num_antigens1;
  this->num_memb_vert=__module_datastruct_MOD_nver;
  this->c = new (nothrow) _RECEPTOR[num_antigens];
  if (this->c == 0) cout<<"\nError in System: memory could not be allocated";
}

void _receptor :: reset_box_dimensions(){
  extern _datainput data;
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();	
  for(int i=0; i<this->num_antigens; i++)  
    c[i].reset_box_dimensions();
}

void _receptor :: dump_conf(ofstream &ofile){
  ofile<<"<ANTIGEN>"<<endl;
	for (int i=0; i< this->num_antigens; i++)   					
		c[i].dump_conf(ofile);
  ofile<<"</ANTIGEN>"<<endl<<endl;
}

void _receptor :: assign_antigen_types(){
  cout<<"test routine "<<endl;
}

void _receptor :: store_trajectory(int window){
	for (int i=0; i<this->num_antigens; i++)
		this->c[i].store_trajectory(window);
}

void _receptor :: dump_trajectory(int window){
	for (int i=0; i<this->num_antigens; i++)
		this->c[i].dump_trajectory(window);
}

_RECEPTOR _receptor :: getantigen(int j){
  return c[j];
}

int _receptor :: getnum_members(){
  return num_antigens;
}

double _receptor :: getL(){
  return L;
}

double _receptor :: getH(){
  return H;
}

void _receptor :: setantigen (_RECEPTOR c_temp, int i){
  c[i] = c_temp;
}

void _receptor :: setxc(int i, double x){
  c[i].setxc(x);
}

double _receptor :: getxc(int i){
  return c[i].getxc();
}

void _receptor :: setyc(int i, double y){
  c[i].setyc(y);
}

double _receptor :: getyc(int i){
  return c[i].getyc();
}

void _receptor :: setzc(int i, double z){
  c[i].setzc(z);
}

double _receptor :: getzc(int i){
  return c[i].getzc();
}

double _receptor :: getantigenxt(int i){                        // get the antigen tip positions
  return c[i].getxt();
}

double _receptor :: getantigenyt(int i){
  return c[i].getyt();
}

double _receptor :: getantigenzt(int i){
  return c[i].getzt();
}

double _receptor :: getxt(int i){                              // get the antigen tip positions
  return c[i].getxt();
}

double _receptor :: getyt(int i){
  return c[i].getyt();
}

double _receptor :: getzt(int i){
  return c[i].getzt();
}

void _receptor :: setantigentip_position(int i) {
  c[i].setantigentip_position();
}

void _receptor :: setantigenxt(int i, double xt1){            // set the antigen tip position
  c[i].setxt(xt1);
}

void _receptor :: setantigenyt(int i, double yt1){
  c[i].setyt(yt1);
}

void _receptor :: setantigenzt(int i, double zt1){
  c[i].setzt(zt1);
}

void _receptor :: setantigenrt(int i, double xt1, double yt1, double zt1){
  c[i].setxt(xt1);
  c[i].setyt(yt1);
  c[i].setzt(zt1);
}

void _receptor :: setradius(int i, double radius){
  c[i].setradius(radius);
}

double _receptor :: getradius(int i){
  return c[i].getradius();
}

double _receptor :: getlength(int i){
  return c[i].getlength();
}

int _receptor :: gettype(int i){
  return c[i].gettype();
}

double _receptor :: get_flexure(int i){
  return c[i].get_flexure();
}

void _receptor :: antigen_translation (int i) {
  c[i].antigen_translation();
}

bool _receptor :: getbonded(int i){
  return c[i].getbonded();
}

void _receptor :: setbonded (int i, bool temp){
  c[i].setbonded(temp);
}

int* _receptor :: getbondedto(int i){
  temp[0] = c[i].getbondedtoves();
  temp[1] = c[i].getbondedtoab();
  return temp;
}

void _receptor :: setbondedto(int i,  int* temp){
  c[i].setbondedtoves(temp[0]);
  c[i].setbondedtoab(temp[1]);
}

void _receptor :: reset_antigen_bond(int i){
  c[i].setbondedtoves(-1);
  c[i].setbondedtoab(-1);
}

bool _receptor ::does_bond_breaks(int i,double tempx,double tempy,double tempz){
  return c[i].does_bond_breaks(tempx,tempy,tempz);
}

bool _receptor ::does_bond_breaks(int i){
  return c[i].does_bond_breaks();
}

void _receptor :: set_theta_phi(int i,double theta, double phi){
  c[i].set_theta(theta);
  c[i].set_phi(phi);
}

void _receptor :: set_theta(int i,double theta){
  c[i].set_theta(theta);
}

void _receptor :: set_phi(int i,double phi){
  c[i].set_phi(phi);
}

double* _receptor :: get_theta_phi(int i){
  this->angle[0]=c[i].get_theta();
  this->angle[1]=c[i].get_phi();
  return angle;
}

double _receptor :: get_theta(int i){
  return c[i].get_theta();
}

double _receptor :: get_phi(int i){
  return c[i].get_phi();
}

void _receptor :: updatebonds(int i){        //To update the bonds
  int *temp;
  double dist;
  double rrx, rry, rrz;
  extern _nanocarrier nc;
  
  if (this->getbonded(i) == 1){          // To break the bond if the antigen goes beyond dist=chem_cutoffu
    temp = this->getbondedto(i);                         
    rrx = nc.getx_antibody(temp[0],temp[1])-this->getantigenxt(i);
    rry = nc.gety_antibody(temp[0],temp[1])-this->getantigenyt(i);
    rrz = nc.getz_antibody(temp[0],temp[1])-this->getantigenzt(i);
    rrx = rrx - this->L*round(rrx/this->L);
    rry = rry - this->L*round(rry/this->L);
    dist = rrx*rrx + rry*rry + rrz*rrz;

    if (dist > this->Dreaction*this->Dreaction){     //if (dist > chem_cutoffu*chem_cutoffu || dist < chem_cutoffl*chem_cutoffl)
      this->setbonded(i,0); 
      nc.setbonded(temp[0], temp[1],0);
      this->setantigentip_position(i);
      this->set_theta_phi(i,0.0,0.0);
      nc.setbondedto(temp[0],temp[1],-1);
      this->reset_antigen_bond(i);
    }
  }
}

void _receptor ::receptor_diffusion(){
  __module_mcsmoves_MOD_antigen_diffusion_on_triangle();
}

void _receptor ::receptor_hopping(){
  __module_mcsmoves_MOD_antigen_hopping_on_vertex();
}
