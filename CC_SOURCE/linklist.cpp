#include "_declarations.h"
#include "_linklist.h"
#include "_datainput.h"       // definition for data reading from input.in
#include "_nanocarrier.h"     // Lays down the datastructure for system variables
#include "_receptor.h"        // Lays down the datastructure for system variables
#include "_membrane.h"        // Lays down the datastructure for system variables
#include "_setup.h"
#include "_fortran_structures.h"

void linkedlist :: init(char ch1, double cellsize){
	extern _datainput data;
	extern _setup a;
	extern bool debug_mode;

	this->ch = ch1;
	this->num_members = a.getnum_members(ch1);
	this->L = data.getperiodic_box_length();
	this->H = data.getperiodic_box_height();
	this->cellIxy = 1.0/cellsize;
	this->cellIz = cellIxy;
	this->Mxy = ceil(cellIxy*L);                               // Mx=My=Mxy=L/getlink_cell_Length
	this->cellIxy = Mxy/L;                                     // We recompute the boxl since M is the nearest bigger integer
	this->Mz = ceil(cellIz*H);                                 // Mz=H/getlink_cell_Length
	this->cellIz = Mz/H;
	this->ncell = Mxy*Mxy*Mz;                                  // Number of cells of dimension (getlink_cell_Length)^3
	__module_datastruct_MOD_nlinkcells = this->ncell;          // pass the total number of linkcells to Fortran  

	if (debug_mode){
  	cout<<str_red<<"\n\t\t<LINKCELL>"<<" "<<ch<<str_black<<endl;
  	cout<<"\t --> (Lx=Ly, Lz): "<<1.0/cellIxy<<",  "<<1.0/cellIz<<" nm"<<endl;
  	cout<<"\t --> (Nx=Ny, Nz): "<<this->Mxy<<" "<<this->Mz<<endl;
  	cout<<"\t --> Nx*Ny*Nz: "<<ncell<<endl;
	}

 	//@ Initialize the _linklist
 	this->lscl = new (nothrow) int[this->num_members];                      
 	this->head = new (nothrow) int[this->ncell];        
 	this->cells_zdir = new(nothrow) int[this->Mz];          // lists all the cells in the  below a given cell (z-direction)
 	if (this->lscl == 0 || this->head==0){
  		cout<<"\n Error: memory could not be allocated\n";
  		exit(1);
 	}
 	fill_n(this->lscl,this->num_members,-1);
 	fill_n(this->head,this->ncell,-1);
 	fill_n(this->cells_zdir,this->Mz,-1);
}

void linkedlist :: reset_box_dimensions(){
 extern _datainput data;
 this->L = data.getperiodic_box_length();
 this->H = data.getperiodic_box_height();	
}

void linkedlist :: calc() {
 extern _membrane me;
 fill_n(this->lscl,this->num_members,-1);
 fill_n(this->head,this->ncell,-1);
 
 for (i=0;i<this->num_members;i++){         // Vector cell index to which this atom belongs
  icell = this->getcellnumber(i,ch);
  if(this->ch =='m'){
   me.setcellno(i,icell);                   // sets the cellnumber for membrane vertex that can be accessed by DTMC code
  }
  if (icell >= this->ncell){
	cout<<str_red;
	cout<<"In "<<__PRETTY_FUNCTION__<<endl<<endl;
	cout<<"Element "<<i<<" of type "<<this->ch<<" is out of bounds"<<endl<<endl;
  	cout<<"computed cell is "<<icell<<", while the total allocated cells are "<<this->ncell<<endl<<endl;
	cout<<"Possible reasons may include: "<<endl;
	cout<<"\t (i) Starting with a membrane file that is not generated for the chosen box dimensions"<<endl;
	cout<<"\n Exiting "<<endl<<str_black;
	exit(1);
  }
  this->lscl[i] = this->head[icell];       // Link to the previous occupant
  this->head[icell] = i;                   // The last one goes to the header
 }
}

void linkedlist :: calc(int sel_mem, double x, double y, double z) {   // Updates the link and header list when only one vesicle is displaced
  extern  _setup a;
  extern _membrane me;
  int cell_new, cell_old, mem; 
  cell_old = this->getcellnumber(x,y,z);
  cell_new = this->getcellnumber(sel_mem, ch);
  if (cell_old == cell_new) return;                     // If the vesicle did not move out of its previous linked cell, don't change link-list
  
  mem = this->head[cell_old];                           // To update the older list
  if (mem == sel_mem){                                  // If the displaced vesicle was in the header
   this->head[cell_old] = this->lscl[sel_mem];          // replace head of old cell with first link
  }
  else {
   do {
     
     if (this->lscl[mem] == sel_mem) {
      this->lscl[mem] = this->lscl[sel_mem];
      this->lscl[sel_mem]= -1; 
      break;
     }                                                   // remove sel_mem from the link list of cell_old
     mem =  this->lscl[mem];     
   } while (mem !=-1);
  }

  mem =  this->head[cell_new];                              // update 
  if (mem == -1){                                           // if the new list is empty then the selected vesicle is the header
   this->head[cell_new] = sel_mem;  
   this->lscl[sel_mem] = -1;
   if (ch == 'm') me.setcellno(sel_mem,cell_new);           // Pass the new cell number to Fortran
  }
  else {
  do {
   if ( this->lscl[mem] == sel_mem) {
   cout<<sel_mem<<" "<< this->ch<<" "<<x<<" "<<y<<" "<<z<<" "<<a.getxc(sel_mem,ch)<<" "<<a.getyc(sel_mem,ch)<<" "<<a.getzc(sel_mem,ch)<<endl;
   cout<<str_blue<<"Link cells before and after "<<cell_old<<" "<<cell_new<<endl;
   cout<<str_red<<"Error: New list already contain moved element."<<str_black; 
   exit(1);
  }

  if (this->lscl[mem] == -1) {
   this->lscl[mem] = sel_mem; 
   if (ch=='m') me.setcellno(sel_mem,cell_new);          // Pass the new cell number to Fortran
    this->lscl[sel_mem] = -1; 
    break;
  }
  
  mem = this->lscl[mem];
  } while (1);
 }
 return;
}

void linkedlist:: locate_allcells(int obj_no){
	int i,mem;
	for(i=0;i<this->ncell;i++){
		mem = this->head[i];
		if (mem == obj_no){
			cout<<"Object "<<obj_no<<" of type "<<ch<<" found in "<<i<<endl;
		}
		else{
		  do{
		   if (this->lscl[mem]==obj_no){
			cout<<"Object "<<obj_no<<" of type "<<ch<<" found in "<<i<<endl;
		   }
		   mem = this->lscl[mem];
		  } while (mem != -1);
		}
	}
}

void linkedlist :: printheadlist(){
	int i;
	for (i=0;i<this->ncell;i++){
	cout<<"Membrane Head list "<<i<<" points to "<<this->head[i]<<endl;
	}
}

void linkedlist :: printlscllist(){
	int i;
	for (i=0;i<this->num_members;i++){
		cout<<"Membrane lscl list "<<i<<" points to "<<lscl[i]<<endl;
		if (this->lscl[i]>10200) {
			cout<<"Found an instance of lscl being greater "<<endl;
			cin.get();
		}
	}
}

int* linkedlist ::  getlscl() {
	return this->lscl;
}

int* linkedlist ::  gethead() {
	return this->head;
}

double linkedlist ::  getcellIxy (){
	return this->cellIxy;
}

double linkedlist :: getcellIz (){
	return this->cellIz;
}

int linkedlist :: getMxy() {
	return this->Mxy;
}

int linkedlist :: getMz() {
	return this->Mz;
}

int* linkedlist ::  getcells(double x, double y, double z) {            // Returns the cell to which the vesicle belongs and also returns the 26 nearest cells
	int i,j,k,l,cell0x, cell0y, cell0z, cellx, celly, cellz;
	cell0x = (int)(x*this->cellIxy);                                // To calculate x,y,z coordinate of linked cell to which the selected_vesicle belongs
	cell0y = (int)(y*this->cellIxy);
	cell0z = (int)(z*this->cellIz);
	l = 0;
	for (i=-1; i<=1; i++){                                         // To calculate neighbouring linked-cells
	 cellz = cell0z + i;
	if (cellz < 0 || cellz >= this->Mz) cellz = cell0z;
	   for (j=-1;j<=1;j++){
	      celly = cell0y + j;
	      if(celly < 0) celly = celly + this->Mxy;                 // Periodic image of linked-cells
	      if(celly >= this->Mxy) celly = celly - this->Mxy;
	      for (k=-1;k<=1;k++){
	 	cellx = cell0x + k;
		if(cellx < 0) cellx = cellx + this->Mxy;
		if(cellx >= this->Mxy) cellx = cellx - this->Mxy;
		this->cell[l] = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
		if (this->cell[l] > this->ncell){
			cout<<"\n-----Error3: Neighbouring cell out of bound----"<<cell[l];
			exit(1);
		}
		l++;
	      }
	   }
	}
	return this->cell;
}

int linkedlist ::  getcell(double x, double y, double z) {                       // Returns the cell to which the vesicle belongs
	int icell;
	int cellx,celly,cellz;
	cellx = (int)(x*this->cellIxy);
	celly = (int)(y*this->cellIxy);
	cellz = (int)(z*this->cellIz);
	icell = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
	return icell;
}

int linkedlist ::  getcellnumber(double x, double y, double z) {                 // Returns the cell to which the vesicle belongs
	int icell;
	int cellx,celly,cellz;
	cellx = (int)(x*this->cellIxy);
	celly = (int)(y*this->cellIxy);
	cellz = (int)(z*this->cellIz);
	icell = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
	return icell;
}

int linkedlist :: getcellnumber(int i,char ch) {                                           // Returns the cell number for a vesicle or membrane i
	extern _setup a;
	int cellx,celly,cellz;
	cellx=(int)(a.getxc(i,ch)*this->cellIxy);
	celly=(int)(a.getyc(i,ch)*this->cellIxy);
	cellz=(int)(a.getzc(i,ch)*this->cellIz);
	icell = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
	return icell;
}

int* linkedlist ::  get2ringcells(int i,char ch) {	                // Returns the cell to which the vesicle belongs and also returns the 124 nearest cells
	extern _setup a;
	int l,cell0x, cell0y, cell0z, cellx, celly, cellz;
	int stack,block,row;
	cell0x = (int)(a.getxc(i,ch)*this->cellIxy);         // To calculate x,y,z coordinate of linked cell to which the selected_vesicle belongs
	cell0y = (int)(a.getyc(i,ch)*this->cellIxy);
	cell0z = (int)(a.getzc(i,ch)*this->cellIz);
	l = 0;
	for (stack=-2;stack<=2;stack++){	                         // To calculate neighbouring lnked-cells
		cellz = cell0z + stack;
		if (cellz < 0 || cellz >= this->Mz) cellz = cell0z;
		for (block=-2; block<=2; block++){
			celly = cell0y + block;
			if(celly < 0) celly = celly + this->Mxy;	          // Periodic image of linked-cells
			if(celly >= Mxy) celly = celly - this->Mxy;
			for (row=-2; row<=2; row++){
				cellx = cell0x + row;
				if(cellx < 0) cellx = cellx + this->Mxy;
				if(cellx >= Mxy) cellx = cellx - this->Mxy;
				this->cell2ring[l] = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
				if (this->cell2ring[l] > this->ncell){
					cout<<"\n-----Error2: Neighbouring cell out of bound----"<<cell2ring[l];
					exit(1);
				}
				l++;
			}
	     	}
	   }
	return this->cell2ring;
}

void linkedlist :: is_member_present_in_cell(int sel_mem){
	int mem;
	int flag=0;
        int sel_cell = this->getcellnumber(sel_mem,ch);
	mem = this->head[sel_cell];
	do {
	    if (mem == sel_mem){
		flag=1;
		break;
	     }
	     mem = this->lscl[mem];
	} while(mem !=-1);

	if (flag == 0){
		cout<<"The chosen member "<<sel_mem<<" of type "<<ch<<" is not found in the cell "<<sel_cell<<endl;
		cin.get();
	}
	
}

void linkedlist :: print_celllist(int sel_cell){
	int mem;
	mem = this->head[sel_cell];
	do {
		cout<<"Sel cell "<<sel_cell<<" has member "<<mem<<endl;
		mem = this->lscl[mem];
	} while(mem !=-1);

}


int* linkedlist ::  getcells_zdirection(int i,char ch) {	        	     // Returns the cell to which the vesicle belongs and also returns the 124 nearest cells
	extern _setup a;
	int l,cell0z, cellx, celly, cellz;
	for (l=0; l<this->Mz; l++)                                               // Initialize the cells array with -1 since the entries can be of varying length
		this->cells_zdir[l]=-1;

	cellx = (int)(a.getxc(i,ch)*this->cellIxy);		// To calculate x,y,z coordinate of linked cell to which the selected_vesicle belongs
	celly = (int)(a.getyc(i,ch)*this->cellIxy);
	cell0z = (int)(a.getzc(i,ch)*this->cellIz);
	l = 0;
	for (cellz=cell0z ; cellz>=0 ; cellz--){
		cells_zdir[l] = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
		l++;
	}
	return this->cells_zdir;
}


int* linkedlist ::  get1ringcells(int i,char ch) {		 // Returns the cell to which the vesicle belongs and also returns the 124 nearest cells
	extern _setup a;
	int l,cell0x, cell0y, cell0z, cellx, celly, cellz;
	int stack,block,row;
	cell0x = (int)(a.getxc(i,ch)*this->cellIxy);		 	  // To calculate x,y,z coordinate of linked cell to which the selected_vesicle belongs
	cell0y = (int)(a.getyc(i,ch)*this->cellIxy);
	cell0z = (int)(a.getzc(i,ch)*this->cellIz);
	l = 0;
	for (stack=-1;stack<=1;stack++){                   	                  			   // To calculate neighbouring lnked-cells
		cellz = cell0z + stack;
		if (cellz < 0 || cellz >= this->Mz) cellz = cell0z;
		for (block=-1 ; block<=1 ; block++){
			celly = cell0y + block;
			if(celly < 0) celly = celly + this->Mxy;	        	          	  // Periodic image of linked-cells
			if(celly >= this->Mxy) celly = celly - this->Mxy;
			for (row=-1 ; row<=1 ; row++){
				cellx = cell0x + row;
				if(cellx < 0) cellx = cellx + this->Mxy;
				if(cellx >= this->Mxy) cellx = cellx - this->Mxy;
				this->cell1ring[l] = cellx + celly*this->Mxy + cellz*this->Mxy*this->Mxy;
				if (this->cell1ring[l] > this->ncell){
					cout<<"\n-----Error1: Neighbouring cell out of bound----"<<this->cell1ring[l];
					exit(1);
				}
				l++;
			}
		 }
	   }
	return this->cell1ring;
}

void _linklist ::  init(){
  extern _nanocarrier nc;
  double minrad = nc.getsoft_radius(0);
  for (int i=1; i<nc.getnum_members(); i++){
    if (nc.getsoft_radius(i)>minrad) minrad = nc.getsoft_radius(i);
  }
  veslink.init('v',1*minrad);
  antigenlink.init('c',1*minrad);
  membvertlink.init('m',1*minrad);
  errorflag[0]=-1111111;
  derrorflag[0]=-111111.0;
}

void _linklist :: init(char ch){
  extern _nanocarrier nc;
  double minrad = nc.getsoft_radius(0);
  for (int i=1; i<nc.getnum_members(); i++){
    if (nc.getsoft_radius(i)>minrad) minrad = nc.getsoft_radius(i);
  }
  if (ch == 'c') antigenlink.init('c',2*minrad);
  if (ch == 'm') membvertlink.init('m',2*minrad);
  if (ch == 'v') veslink.init('v',2*minrad);
}


void _linklist :: reset_box_dimensions(){
		veslink.reset_box_dimensions();
		antigenlink.reset_box_dimensions();
		membvertlink.reset_box_dimensions();
}

void _linklist ::  calc(char ch) {
	if (ch == 'v')
		veslink.calc();
	if (ch == 'c')
		antigenlink.calc();
	if (ch=='m')
		membvertlink.calc();
}

void _linklist :: calc(int sel_mem, double x, double y, double z, char ch) {
	if (ch == 'v')
	{
		veslink.calc(sel_mem, x, y, z);
	}
	if (ch == 'c'){
		antigenlink.calc(sel_mem, x, y, z);
	}
	if (ch=='m'){
		membvertlink.calc(sel_mem,x,y,z);
	}
}

void _linklist:: locate_allcells(int obj_no,char ch){
	if (ch=='c'){
		antigenlink.locate_allcells(obj_no);
	}
	if (ch=='v'){
		veslink.locate_allcells(obj_no);
	}
	if (ch=='m'){
		membvertlink.locate_allcells(obj_no);
	}
}


void _linklist :: printlscllist(char ch){
	if (ch == 'v')
		return veslink.printlscllist();
	if (ch == 'c')
		return antigenlink.printlscllist();
	if (ch == 'm')
		return membvertlink.printlscllist();
}

void _linklist :: printheadlist(char ch){
	if (ch == 'v')
		return veslink.printheadlist();
	if (ch == 'c')
		return antigenlink.printheadlist();
	if (ch == 'm')
		return membvertlink.printheadlist();
}


int* _linklist ::  getlscl(char ch) {
	if (ch == 'v')
		return veslink.getlscl();
	if (ch == 'c')
		return antigenlink.getlscl();
	if (ch == 'm')
		return membvertlink.getlscl();
	throw_errormessage(__PRETTY_FUNCTION__,ch);                         // Throws a error message with an option to exit or continue (the same is done throughout)
	return this->errorflag;                                                   // Returns a predefined value if continued.  
}

int* _linklist :: gethead(char ch) {
	if (ch == 'v')
		return veslink.gethead();
	if (ch == 'c')
		return antigenlink.gethead();
	if (ch == 'm')
		return membvertlink.gethead();
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag;
}

double _linklist :: getcellIxy(char ch){
	if (ch == 'v')
		return veslink.getcellIxy();
	if (ch == 'c')
		return antigenlink.getcellIxy();
	if (ch == 'm')
		return membvertlink.getcellIxy();
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->derrorflag[0];
}

double _linklist :: getcellIz(char ch){
	if (ch == 'v')
		return veslink.getcellIz();
	if (ch == 'c')
		return antigenlink.getcellIz();
	if (ch == 'm')
		return membvertlink.getcellIz();
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->derrorflag[0];
}

int _linklist ::  getMxy(char ch){
	if (ch == 'v')
		return veslink.getMxy();
	if (ch == 'c')
		return antigenlink.getMxy();
	if (ch == 'm')
		return membvertlink.getMxy();
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag[0];
}

int _linklist :: getMz(char ch){
	if (ch == 'v')
		return veslink.getMz();
	if (ch == 'c')
		return antigenlink.getMz();
	if (ch == 'm')
		return membvertlink.getMz();
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag[0];
}

int* _linklist ::  getcells(double x, double y, double z, char ch) {
	if (ch == 'v')
	 return veslink.getcells(x, y, z);
	if (ch == 'c')
	 return antigenlink.getcells(x, y, z);
	if (ch == 'm')
	 return membvertlink.getcells(x,y,z);
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag;
}

int _linklist ::  getcell(double x, double y, double z, char ch) {
	if (ch == 'v')
		return veslink.getcell(x, y, z);
	if (ch == 'c')
		return antigenlink.getcell(x, y, z);
	if (ch == 'm')
		return membvertlink.getcell(x,y,z);
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag[0];
}


int* _linklist ::  get2ringcells(int i, char ch) {
	if (ch == 'v')
		return veslink.get2ringcells(i,ch);
	if (ch == 'c')
		return antigenlink.get2ringcells(i,ch);
	if (ch == 'm')
		return membvertlink.get2ringcells(i,ch);
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag;
}

int* _linklist ::  get1ringcells(int i, char ch) {
	if (ch == 'v')
		return veslink.get1ringcells(i,ch);
	if (ch == 'c')
		return antigenlink.get1ringcells(i,ch);
	if (ch == 'm')
		return membvertlink.get1ringcells(i,ch);
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag;
}
	
int* _linklist ::  getcells_zdirection(int i,char ch) {
	if (ch == 'v')
		return veslink.getcells_zdirection(i,ch);
	if (ch == 'c')
		return antigenlink.getcells_zdirection(i,ch);
	if (ch == 'm')
		return membvertlink.getcells_zdirection(i,ch);
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag;
}
	
	
int _linklist ::  getcellnumber(int i, char ch) {
	if (ch == 'v')
		return veslink.getcellnumber(i,ch);
	if (ch == 'c')
		return antigenlink.getcellnumber(i,ch);
	if (ch == 'm')
		return membvertlink.getcellnumber(i,ch);
	throw_errormessage(__PRETTY_FUNCTION__,ch);
	return this->errorflag[0];
}
	
void _linklist ::  is_member_present_in_cell(int i, char ch) {
	if (ch == 'v')
		return veslink.is_member_present_in_cell(i);
	if (ch == 'c')
		return antigenlink.is_member_present_in_cell(i);
	if (ch == 'm')
		return membvertlink.is_member_present_in_cell(i);
}
	
void _linklist :: print_celllist(int i, char ch) {
	if (ch == 'v')
		return veslink.print_celllist(i);
	if (ch == 'c')
		return antigenlink.print_celllist(i);
	if (ch == 'm')
		return membvertlink.print_celllist(i);
}
	
	
	
void _linklist :: throw_errormessage(string input, char ch){
		cout<<"\nError in function "<<input<<endl;
		cout<<"Wrong character input : passed value of ch is '"<<ch<<"'"<<endl;
		cout<<"enter 'c' to continue and 'e' to exit "<<endl;
		char action=cin.get();
		if (action=='e') exit(0);
}

