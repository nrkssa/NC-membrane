#include "_declarations.h"
#include "_constant.h"
#include "_parser.h"
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_membrane.h"
#include "_potential.h"
#include "_datainput.h"                            // definition for data reading from input.in
#include "MersenneTwister-ran2.h"                  //  header file for random number generator
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_linklist.h"
#include "_selfavoidance.h"
#include "_datawriter.h"
#include "_debugger.h"
#include "_interaction_parameters.h"

//@ Antibody datastructure

void _ANTIBODY :: init(int index){
     this->parse_abparameters(index);
}

void _ANTIBODY :: init(int nab, double* abl, double* abr, double* abc){
	extern bool debug_mode;
	this->nabtype = nab;
	this->abradius = new (nothrow) double[this->nabtype];
	this->ablength = new (nothrow) double[this->nabtype];
	this->ab_type_conc = new (nothrow) double[this->nabtype];
	
	for (int i=0; i<nab; i++){
		this->abradius[i] = abr[i];
		this->ablength[i] = abl[i];
		this->ab_type_conc[i] = abc[i];
	}
	
	if (debug_mode)
		this->print_message();
}

//@ read the antibody parameters enclosed between <antibody> and </antibody>
void _ANTIBODY :: parse_abparameters(int index){
	extern bool debug_mode;
	stringstream filename;
	_PARSER parse;
	char DELIMITER=',';
	filename<<"../PARAMETERS/nc-"<<index<<".ncinp";
	auto parsedata = parse.PARSEFILE(filename.str().c_str(),DELIMITER);
	
	std::array<int,MAX_FILE_ENTRY> parsval = parse.PARSEINT(parsedata,"ab_type",1);
	this->nabtype = parsval[0];
	this->abradius = new (nothrow) double[this->nabtype];
	this->ablength = new (nothrow) double[this->nabtype];
	this->ab_type_conc = new (nothrow) double[this->nabtype];

	std::array<double,MAX_FILE_ENTRY> parsval1 = parse.PARSEDOUBLE(parsedata,"ab_radius",this->nabtype);
	std::array<double,MAX_FILE_ENTRY> parsval2 = parse.PARSEDOUBLE(parsedata,"ab_length",this->nabtype);
	std::array<double,MAX_FILE_ENTRY> parsval3 = parse.PARSEDOUBLE(parsedata,"ab_type_conc",this->nabtype);

	for (int i=0; i<this->nabtype; i++){
	 this->abradius[i] = parsval1[i];
	 this->ablength[i] = parsval2[i];
	 this->ab_type_conc[i] = parsval3[i];  
	}
	
	if (debug_mode)
		this->print_message();
}

void _ANTIBODY :: print_message(){
	cout<<str_green<<"\n\t\t <ANTIBODY> \n"<<str_black<<endl;
	cout<<"\t\t\t --> Types of antibodies "<<this->nabtype<<endl;		
	cout<<"\t\t\t --> Radius of AB type "<<"\t";
	for (int i=0; i<this->nabtype; i++)
		cout<<i<<" "<<this->abradius[i]<<"\t";
	cout<<endl;
	cout<<"\t\t\t --> Length of AB type "<<"\t";
	for (int i=0; i<this->nabtype; i++)
		cout<<i<<" "<<this->ablength[i]<<"\t";
	cout<<endl;
}

int _ANTIBODY:: get_numtype(){
	return this->nabtype;
}

double _ANTIBODY:: get_ablength(int i){
	return this->ablength[i];
}

double _ANTIBODY :: get_abradius(int i){
	return this->abradius[i];
}

//@ Nanocarrier datastructure
void _NANOCARRIER :: init(int i){
  extern _datainput data;
  extern bool debug_mode;
  
  this->index = i;
  this->temperature = data.gettemperature();
  this->kBT = (kb*temperature);
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();
  
  if (debug_mode)
	  cout<<str_red<<"\n\t <NC#"<<this->index<<">"<<str_black<<endl;
  
  abclass.init(this->index);
  this->verptr = &__module_datastruct_MOD_ver;
  this->ncpointer= &__module_datastruct_MOD_nc_f + this->index;    //@ Point to the corresponding NC.
 
  this->xc = this->yc = this->zc = 0.0;
  this->anglex = this->angley = this->anglez = 0.0;
  this->multivalency = this->nbonded = this->nbroken = 0;         //@ counter to monitor the number of bonds
  this->adjust_interval = 1000;
  this->rstep_ub = this->rstep_b = 0.1;
  fill_n(this->trans_counter,4,0);
  fill_n(this->rot_counter,4,0);
  fill_n(this->bond_counters,3,0); 
  
  this->parse_nc_parameters(); 
  this->setnum_ab();          //@ Initialize all the required arrays for NC and ab 
  this->create_antibody();
 
  
  if (this->index == 0){
  cout<<this->getx_antibody(0)<<" "<<this->gety_antibody(0)<<" "<<this->getz_antibody(0)<<endl;
  cout<<this->getx_antibody(1)<<" "<<this->gety_antibody(1)<<" "<<this->getz_antibody(1)<<endl;
  }
  
  
  this->soft_radius = this->compute_max_soft_radius();
  this->tstep_ub = this->tstep_b = this->radius/100.0;
  this->update_fortran_pointers();
  
  if (debug_mode) cout<<endmarker<<endl;
}

void _NANOCARRIER :: init(int ncindex, ifstream &ifile) {
	extern bool debug_mode;
	string temps;
	stringstream tstring;
	bool foundflag = false;
	
	if (debug_mode)
		cout<<str_red<<"\n\t <NC#"<<ncindex<<">"<<str_black<<endl;
	
	this->index = ncindex;
	tstring<<"<NC#"<<this->index<<">";
	
	while (getline(ifile,temps)){
		if (temps.compare(tstring.str()) == 0){
			foundflag = true;
			break;
		}
	}
	
	int tind;
	ifile>>temps>>tind>>this->num_ab>>this->num_abtype;
	this->verptr = &__module_datastruct_MOD_ver;
	this->ncpointer= &__module_datastruct_MOD_nc_f + this->index;    //@ Point to NC.
	this->setnum_ab();
	
	double t1[this->num_abtype],t2[this->num_abtype],t3[this->num_abtype];
	this->tnum_ab_type = new(nothrow) int[this->num_abtype];
	this->ab_type_conc = new (nothrow) double[this->num_abtype];
	
	for (int i=0; i<this->num_abtype; i++){
		ifile>>tind>>t1[i]>>t2[i]>>t3[i]>>this->tnum_ab_type[i];
		this->ab_type_conc[i] = t3[i];
	}
	
	abclass.init(this->num_abtype,t1,t2,t3);
	
	ifile>>this->xc>>this->yc>>this->zc>>this->anglex>>this->angley>>this->anglez>>this->radius>>this->soft_radius;
	ifile>>this->L>>this->H>>this->kBT;
	ifile>>this->ncshape>>this->ab_arrangement;
	ifile>>this->rstart>>this->hstart;
	ifile>>this->bias_mode>>this->kbias>>this->biasref>>this->biasdir>>this->binsize>>this->lambda;
	ifile>>this->trans_counter[0]>>this->trans_counter[1]>>this->trans_counter[2]>>this->trans_counter[3];
	ifile>>this->rot_counter[0]>>this->rot_counter[1]>>this->rot_counter[2]>>this->rot_counter[3];
	ifile>>this->stepsize[0]>>this->stepsize[1]>>this->stepsize[2]>>this->stepsize[3];
	ifile>>this->multivalency>>this->nbonded>>this->nbroken>>this->bond_counters[0]>>this->bond_counters[1]>>this->bond_counters[2];
	ifile>>this->tstep_ub>>this->tstep_b>>this->rstep_ub>>this->rstep_b>>this->adjust_interval;
	ifile>>temps;
	for (int i=0; i<this->num_ab; i++){
		ifile>>this->x_antibody[i]>>this->y_antibody[i]>>this->z_antibody[i]>>this->xf_antibody[i]>>this->yf_antibody[i]>>this->zf_antibody[i]>>this->ab_type[i]>>this->ab_size[i]>>this->ab_radius[i]>>this->soft_rad_ab[i]>>this->bonded[i]>>this->bondedto[i];
	}
	
	this->update_fortran_pointers();
	
	if (!foundflag){
		cout<<"Could not find data for NC "<<this->index<<endl;
		cout<<"Exiting at "<<__PRETTY_FUNCTION__<<endl;
		exit(1);
	}
	
	if(debug_mode){
		cout<<"\n \t ==> Initialized NC#"<<this->index<<"  at ("<<this->xc<<", "<<this->yc<<", "<<this->zc<<")"<<endl;
		cout<<"   \t ==> Multivalency "<<this->multivalency<<endl;
		if (debug_mode) cout<<endmarker<<endl;
	}
	
	this->setxy_npb_antibody();
}

void _NANOCARRIER :: dump_conf(ofstream &ofile){
	ofile<<"<NC#"<<this->index<<">"<<endl;
	ofile<<"#NC "<<this->index<<" "<<this->num_ab<<" "<<this->num_abtype<<endl;
	for (int i=0; i<this->num_abtype; i++)
		ofile<<i<<" "<<abclass.get_ablength(i)<<" "<<abclass.get_abradius(i)<<" "<<this->ab_type_conc[i]<<" "<<this->tnum_ab_type[i]<<endl;
	ofile<<this->xc<<" "<<this->yc<<" "<<this->zc<<" "<<this->anglex<<" "<<this->angley<<" "<<this->anglez<<" "<<this->radius<<" "<<this->soft_radius<<endl;
	ofile<<this->L<<" "<<this->H<<" "<<this->kBT<<endl;
	ofile<<this->ncshape<<" "<<this->ab_arrangement<<" "<<endl;
	ofile<<this->rstart<<" "<<this->hstart<<endl;
	ofile<<this->bias_mode<<" "<<this->kbias<<" "<<this->biasref<<" "<<this->biasdir<<" "<<this->binsize<<" "<<this->lambda<<endl;
	ofile<<this->trans_counter[0]<<" "<<this->trans_counter[1]<<" "<<this->trans_counter[2]<<" "<<this->trans_counter[3]<<endl;
	ofile<<this->rot_counter[0]<<" "<<this->rot_counter[1]<<" "<<this->rot_counter[2]<<" "<<this->rot_counter[3]<<endl;
	ofile<<this->stepsize[0]<<" "<<this->stepsize[1]<<" "<<this->stepsize[2]<<" "<<this->stepsize[3]<<endl;
	ofile<<this->multivalency<<" "<<this->nbonded<<" "<<this->nbroken<<" "<<this->bond_counters[0]<<" "<<this->bond_counters[1]<<" "<<this->bond_counters[2]<<endl;
	ofile<<this->tstep_ub<<" "<<this->tstep_b<<" "<<this->rstep_ub<<" "<<this->rstep_b<<" "<<this->adjust_interval<<endl;
	
	ofile<<"<AB>"<<endl;
	for (int i=0; i<this->num_ab; i++){
		ofile<<this->x_antibody[i]<<" "<<this->y_antibody[i]<<" "<<this->z_antibody[i]<<" "<<this->xf_antibody[i]<<" "<<this->yf_antibody[i]<<" "<<this->zf_antibody[i]<<" "<<this->ab_type[i]<<" "<<this->ab_size[i]<<" "<<this->ab_radius[i]<<" "<<this->soft_rad_ab[i]<<" "<<this->bonded[i]<<" "<<this->bondedto[i]<<endl;
	}
	ofile<<"</AB>"<<endl;
	ofile<<"</NC#"<<this->index<<">"<<endl<<endl;
}

void _NANOCARRIER :: reset_box_dimensions(){
	extern _datainput data;
	this->L = data.getperiodic_box_length();
	this->H = data.getperiodic_box_height();	
}

void _NANOCARRIER :: print_address(){
  cout<<this->x_antibody<<" "<<&this->x_antibody[this->num_ab-1]<<endl;
  cout<<this->y_antibody<<" "<<&this->y_antibody[this->num_ab-1]<<endl;
  cout<<this->z_antibody<<" "<<&this->z_antibody[this->num_ab-1]<<endl;
}

void _NANOCARRIER :: setnum_ab() {
	this->x_antibody = new (nothrow) double[this->num_ab];
	this->y_antibody = new (nothrow) double[this->num_ab];
	this->z_antibody = new (nothrow) double[this->num_ab];
	this->xtold = new (nothrow) double[this->num_ab];
	this->ytold = new (nothrow) double[this->num_ab];
	this->ztold = new (nothrow) double[this->num_ab];
	this->xf_antibody = new (nothrow) double[this->num_ab];
	this->yf_antibody = new (nothrow) double[this->num_ab];
	this->zf_antibody = new (nothrow) double[this->num_ab];
	this->xbold = new (nothrow) double[this->num_ab];
	this->ybold = new (nothrow) double[this->num_ab];
	this->zbold = new (nothrow) double[this->num_ab];
	this->x_npb_antibody = new (nothrow) double[this->num_ab];
	this->y_npb_antibody = new (nothrow) double[this->num_ab];
	this->xf_npb_antibody = new (nothrow) double[this->num_ab];
	this->yf_npb_antibody = new (nothrow) double[this->num_ab];
	this->bonded = new (nothrow) bool[this->num_ab];
	this->bondedto = new (nothrow)  int[this->num_ab];
	if (this->x_antibody == 0 || this->y_antibody ==0 || this->z_antibody == 0 || this->bonded == 0 || this->bondedto == 0)
	 cout<<"\n Error in System: memory could not be allocated";
	if (this->xf_antibody == 0 || this->yf_antibody ==0 || this->zf_antibody == 0)
	 cout<<"\n Error in System: memory could not be allocated";
	
	this->soft_rad_ab = new (nothrow) double[this->num_ab];
	this->ab_size = new (nothrow) double[this->num_ab];
	this->ab_radius = new (nothrow) double[this->num_ab];
	this->ab_type = new (nothrow) int[this->num_ab];	
	if (this->soft_rad_ab ==0 || this->ab_size==0 || this->ab_radius ==0 || this->ab_type==0)
	 cout<<"\n Error in System: memory could not be allocated";

	fill_n(this->bonded,this->num_ab,false);
	fill_n(this->bondedto,this->num_ab,-1);
}

void _NANOCARRIER :: parse_nc_parameters(){
	stringstream filename;
	extern bool debug_mode;
	extern _datainput data;
	_PARSER parse;
	char DELIMITER = ',';
	
	runparam RP = data.get_runparameters();                              // get the run parameters from datainput.cpp
	
	if (debug_mode)
	  cout<<str_blue<<"\n \t\t =========> Parsing parameters for NC "<<this->index<<str_black<<endl<<endl;
	
	filename<<"../PARAMETERS/nc-"<<this->index<<".ncinp";
	auto parsedata = parse.PARSEFILE(filename.str().c_str(),DELIMITER);
	std::array<std::string,MAX_FILE_ENTRY> pars_str;
	std::array<double,MAX_FILE_ENTRY> pars_double;
	std::array<int,MAX_FILE_ENTRY> pars_int;

	//@ shape of the NC
	pars_str = parse.PARSESTRING(parsedata,"nc_shape",1);
	this->ncshape = pars_str[0];
	//@ particle radius
	pars_double = parse.PARSEDOUBLE(parsedata,"nc_radius",1);
	this->radius = pars_double[0];
	
	//@ total number of antibodies
	pars_int = parse.PARSEINT(parsedata,"nc_num_ab",1);
	this->num_ab = pars_int[0];
	
	//@ types of antibodies
	pars_int = parse.PARSEINT(parsedata,"ab_type",1);
	this->num_abtype = pars_int[0];
	this->ab_type_conc = new (nothrow) double[this->num_abtype];
	this->tnum_ab_type = new (nothrow) int[this->num_abtype];
	if (this->ab_type_conc ==0 || this->tnum_ab_type ==0){
		cout<<"Could not assign AB concentration arrays "<<endl;
		exit(1);
	}
	
	//@ concentration fo each type of antibodies
	pars_double = parse.PARSEDOUBLE(parsedata,"ab_type_conc",this->num_abtype);
	for (size_t i=0; i<(size_t) this->num_abtype; i++){
	  this->ab_type_conc[i] = pars_double[i];
	}
	
	pars_double = parse.PARSEDOUBLE(parsedata,"nc_rstart",1);
	this->rstart = pars_double[0];
	pars_double = parse.PARSEDOUBLE(parsedata,"nc_hstart",1);
	this->hstart = pars_double[0];
	
	
	//@ how to pattern the abs on the NC
	pars_str = parse.PARSESTRING(parsedata,"ab_mode",1);
        this->ab_arrangement = pars_str[0];
	
	//@ set the absolute number of each antibody based on their concentration
	int tempnum=0;
	for(int i=0; i<this->num_abtype-1; i++){
		this->tnum_ab_type[i] = round(this->ab_type_conc[i]*this->num_ab);
		tempnum += this->tnum_ab_type[i];
	}
	this->tnum_ab_type[this->num_abtype-1] = this->num_ab - tempnum;
	
	
	//@@ parse details of how the NC should be simulated
	//@ bias mode for the NC
	pars_str = parse.PARSESTRING(parsedata,"bias_mode",1);
	const char *temp=pars_str[0].c_str();
	this->bias_mode = temp[0];
	pars_double = parse.PARSEDOUBLE(parsedata,"binsize",1);
	double bsize = pars_double[0];
	pars_double = parse.PARSEDOUBLE(parsedata,"biasref",1);
	double bref = pars_double[0];
	pars_double = parse.PARSEDOUBLE(parsedata,"biasstrength",1);
	double bias_strength = pars_double[0];                         // kbias = nKT/(binsize^2)
	pars_str = parse.PARSESTRING(parsedata,"biasdir",1);
	string bdir = pars_str[0];
	pars_double = parse.PARSEDOUBLE(parsedata,"TI_lambda",1);
	double lambda_zero = pars_double[0];
	
	if (this->bias_mode == 'N'){
	  this->kbias = 0.0;
	  this->biasref = 0.0;
	  this->biasdir ="None";
	  this->binsize = 0.0;
	  this->lambda =1.0;
	}
	
	else if ((this->bias_mode == 'Z') || (this->bias_mode == 'H') || (this->bias_mode == 'T')){
	  this->binsize = bsize;
	  this->biasdir = bdir;
	  
	  this->kbias = 2.0*bias_strength*this->kBT/pow(this->binsize,2);         // k=2*kbt/(dx**2)
	  
	  if ((this->bias_mode == 'Z') || (this->bias_mode == 'H')){
	    if (this->biasdir.compare("approach") == 0){
	      this->biasref = bref - RP.win_start*bsize;
	    }
	    else if (this->biasdir.compare("recede") == 0){
	      this->biasref = bref + RP.win_start*this->binsize;
	    }
	    else {
	      cout<<"Incorrect bias direction "<<bdir<< " for NC "<<this->index<<endl;
	      cout<<"Options are 'approach'/'recede'"<<endl;
	      cout<<"Exiting - Restart with modified "<<filename.str()<<endl;
	      exit(1);
	     }
	    this->lambda = 1.0;
	  }
	  
	  else if (this->bias_mode == 'T'){
	   this->biasref = bref;
	   
	   if (RP.total_proc == 1){
	     this->lambda = lambda_zero;
	   }
	   else{
	    double dlambda = (1.0-lambda_zero)/(RP.total_proc-1);
	    this->lambda = lambda_zero + RP.curr_proc*dlambda;
	   }
	  }
	}
	
	if (debug_mode){
	//@ print the parsed value to the screen
	cout<<"\t\t Nanoparticle number: "<<this->index+1<<endl;
	cout<<"\t\t Shape: "<<this->ncshape<<endl;
	cout<<"\t\t Radius: "<<this->radius<<endl;	
	cout<<"\t\t Number of AB: "<<this->num_ab<<endl;
	cout<<"\t\t Types of AB: "<<this->num_abtype<<endl;
	cout<<"\t\t Individual AB numbers: ";
	for (int i=0; i<this->num_abtype; i++)
	  cout<<this->tnum_ab_type[i]<<"\t";
	cout<<endl;
	
	cout<<"\t\t Bias_mode: "<<this->bias_mode<<endl;
	cout<<"\t\t Bias_strength: "<<this->kbias<<endl;
	cout<<"\t\t Bias_ref: "<<this->biasref<<endl;
	cout<<"\t\t Bias binsize: "<<this->binsize<<endl;
	cout<<"\t\t Bias direction: "<<this->biasdir<<endl;
	cout<<"\t\t Lambda value: "<<this->lambda<<endl;
	cout<<str_blue<<"\n\t\t =========> Completed Parsing "<<this->index<<str_black<<endl<<endl;
	}
}

void _NANOCARRIER :: create_antibody(){
	int i;
	stringstream ncfile;
	bool *tempabfill;
	double *xpos,*ypos,*zpos,sradius;
	extern MTRand mtrand1;
	
	tempabfill = new (nothrow) bool[this->num_ab];
	xpos = new(nothrow) double[this->num_ab];
	ypos = new(nothrow) double[this->num_ab];
	zpos = new(nothrow) double[this->num_ab];
	

	if ((xpos ==0) || (ypos==0) || (zpos==0) || (tempabfill==0)){
	  cout<<"Error in alllocation of xpos, ypos, and zpos | line number nanocarrier.cpp:140"<<endl;
	  exit(1);
	}

	ncfile<<"../PARAMETERS/"<<this->ncshape<<this->num_ab<<".in";                      // read from the data file
	ifstream infile (ncfile.str().c_str(),ios::in);                         // read from the data file
	if (!infile) {
	  cout << "Cannot open sphere data file."<<endl;
	  exit(1);
	}
	
	for (i=0;i<num_ab;i++) {                                                // read the discrete coordinates for the sphere
	 infile >> xpos[i];
	 infile >> ypos[i];
	 infile >> zpos[i];
	}
	
	infile.close();
	fill_n(tempabfill,this->num_ab,false);                               //@ set tempabfill to false 
	this->shift_to_COM_frame(this->num_ab,xpos,ypos,zpos);               //@ shift coordinates to center of mass

        //@ Only one type of antibody 
	
	if (this->num_abtype == 1){
	  for (int i=0; i < this->num_ab; i++){
	    sradius = this->radius + this->abclass.get_ablength(0);                  //@ compute the soft radius for each antibody
	    this->soft_rad_ab[i] = sradius;                                     //@ assign soft radius
	    this->ab_size[i] = this->abclass.get_ablength(0);
	    this->ab_radius[i] = this->abclass.get_abradius(0);
	    this->ab_type[i] = 0;                                               //@ assign antibody type 
	    this->set_base_tip_antibody(i,xpos[i],ypos[i],zpos[i]);             //@ set the tip and base position 
	  }
	}

	//@ uniform distribution of multiple types
	else if ((this->ab_arrangement.compare("random") == 0) && (this->num_abtype >1)){
	  for (int i=0; i < this->num_abtype; i++){
	    int ii = 0;
	    while (ii < this->tnum_ab_type[i]){
	      int iii = mtrand1.randInt(this->num_ab-1);                               //@ choose a random vertex
	      if (!tempabfill[iii]){                                             //@ if not assigned already
		sradius = this->radius + this->abclass.get_ablength(i);                  //@ compute the soft radius for each antibody
		this->soft_rad_ab[iii] = sradius;                                     //@ assign soft radius
		this->ab_size[iii] = this->abclass.get_ablength(i);
		this->ab_radius[iii] = this->abclass.get_abradius(i);
		this->ab_type[iii] = i;                                               //@ assign antibody type 
		this->set_base_tip_antibody(iii,xpos[iii],ypos[iii],zpos[iii]);  //@ set the tip and base position 
		tempabfill[iii]=true;                                            //@ set flag to assigned
		ii += 1;
	        }
	     }
	  }
	}

	//@ Dipole like patterns with two types of ABs
	else if ((this->ab_arrangement.compare("dipole")==0) && (this->num_abtype ==2)){
	  int ii =0;
	  for (int i=0; i<this->num_ab; i++){
	    if (xpos[i]<0.0) {
		sradius = this->radius + this->abclass.get_ablength(0);                  //@ compute the soft radius for each antibody
		this->ab_type[i] = 0;
		this->soft_rad_ab[i]=sradius;
		this->ab_size[i] = this->abclass.get_ablength(0);
		this->ab_radius[i] = this->abclass.get_abradius(0);
		this->ab_type[i] =0;
		ii += 1;
	    }
	    else {
		sradius = this->radius + this->abclass.get_ablength(1);                  //@ compute the soft radius for each antibody
		this->ab_type[i] = 1;  
		this->soft_rad_ab[i]=sradius;
		this->ab_size[i] = this->abclass.get_ablength(1);
		this->ab_radius[i] = this->abclass.get_abradius(1);
		this->ab_type[i] =1;
	    }
	this->set_base_tip_antibody(i,xpos[i],ypos[i],zpos[i]);
	  }
	  this->tnum_ab_type[0] = ii;
	  this->tnum_ab_type[1] = this->num_ab-ii;                             //@ if in case the number of 0 and 1 are not exactly equal
	}
	
	else if ((this->ab_arrangement.compare("angle")==0) && (this->num_abtype ==2)){
		int ii =0;
		double anglezero = (180-45)*PI/180;
		for (int i=0; i<this->num_ab; i++){
			double ttheta = acos(zpos[i]);
			if (ttheta > anglezero){
				sradius = this->radius + this->abclass.get_ablength(0);                  //@ compute the soft radius for each antibody
				this->ab_type[i] = 0;
				this->soft_rad_ab[i]=sradius;
				this->ab_size[i] = this->abclass.get_ablength(0);
				this->ab_radius[i] = this->abclass.get_abradius(0);
				ii += 1;
			}
			else {
				sradius = this->radius + this->abclass.get_ablength(1);                  //@ compute the soft radius for each antibody
				this->ab_type[i] = 1;  
				this->soft_rad_ab[i]=sradius;
				this->ab_size[i] = this->abclass.get_ablength(1);
				this->ab_radius[i] = this->abclass.get_abradius(1);
			}
			this->set_base_tip_antibody(i,xpos[i],ypos[i],zpos[i]);
		}
		this->tnum_ab_type[0] = ii;
		this->tnum_ab_type[1] = this->num_ab-ii;                             //@ if in case the number of 0 and 1 are not exactly equal
	}

	
	//@ quadrupole like patterns with 2 types of antibodies
	else if ((this->ab_arrangement.compare("quadrapole")==0) && (this->num_abtype ==2)){
	  int ii =0;
	  for (int i=0; i<this->num_ab; i++){
	    double tphi = atan2(ypos[i],xpos[i])+PI;
	    double ttheta = acos(zpos[i]);
	      if (cos(ttheta)-0.7*cos(2*tphi) < 0){
	        sradius = this->radius + this->abclass.get_ablength(0);                  //@ compute the soft radius for each antibody
		this->ab_type[i] = 0;
		this->soft_rad_ab[i]=sradius;
		this->ab_size[i] = this->abclass.get_ablength(0);
		this->ab_radius[i] = this->abclass.get_abradius(0);
		this->ab_type[i] =0;
		ii += 1;
	    }
	    else {
		sradius = this->radius + this->abclass.get_ablength(1);                  //@ compute the soft radius for each antibody
		this->ab_type[i] = 1;  
		this->soft_rad_ab[i]=sradius;
		this->ab_size[i] = this->abclass.get_ablength(1);
		this->ab_radius[i] = this->abclass.get_abradius(1);
		this->ab_type[i] =1;
	    }
	    this->set_base_tip_antibody(i,xpos[i],ypos[i],zpos[i]);
	  }
	  this->tnum_ab_type[0] = ii;
	  this->tnum_ab_type[1] = this->num_ab-ii;                             //@ if in case the number of 0 and 1 are not exactly equal
	}

	else if ((this->ab_arrangement.compare("quadrapole-lower")==0) && (this->num_abtype ==2)){
	  int ii =0;
	  for (int i=0; i<this->num_ab; i++){
	    double tphi = atan2(ypos[i],xpos[i])+PI;
	    double ttheta = acos(zpos[i]);
	      if (cos(ttheta) < 0.0 && cos(2*tphi)< 0.0){
	        sradius = this->radius + this->abclass.get_ablength(0);                  //@ compute the soft radius for each antibody
		this->ab_type[i] = 0;
		this->soft_rad_ab[i]=sradius;
		this->ab_size[i] = this->abclass.get_ablength(0);
		this->ab_radius[i] = this->abclass.get_abradius(0);
		this->ab_type[i] =0;
		ii += 1;
	    }
	    else {
		sradius = this->radius + this->abclass.get_ablength(1);                  //@ compute the soft radius for each antibody
		this->ab_type[i] = 1;  
		this->soft_rad_ab[i]=sradius;
		this->ab_size[i] = this->abclass.get_ablength(1);
		this->ab_radius[i] = this->abclass.get_abradius(1);
		this->ab_type[i] =1;
	    }
	    this->set_base_tip_antibody(i,xpos[i],ypos[i],zpos[i]);
	  }
	  this->tnum_ab_type[0] = ii;
	  this->tnum_ab_type[1] = this->num_ab-ii;                             //@ if in case the number of 0 and 1 are not exactly equal
	}
	//@ other types
	else{
	  cout<< str_red;
	  cout<< " The specified pattern `"<<str_blue<<this->ab_arrangement<<str_red<<"` for NC# "<<this->index<< " is not defined"<<endl;
	  cout<< " Continue after adding this pattern type to:  `nanocarrier.cpp:_NANOCARRIER::create_antibody()`"<<endl;
	  cout<< " Receiving exit signal at `nanocarrier.cpp:422`"<<endl;
	  cout<< str_black;
	  exit(1);
	}
	delete[] xpos;
	delete[] ypos;
	delete[] zpos;
	delete[] tempabfill;
}

void _NANOCARRIER :: update_fortran_pointers(){
  (this->ncpointer)->kbias = (this->kbias/this->kBT);             //@ k*l^2/kBT
  (this->ncpointer)->biasref = this->biasref;
  (this->ncpointer)->lambda = this->lambda;
  (this->ncpointer)->multivalency = this->multivalency;
  (this->ncpointer)->bias_mode = this->bias_mode;
  (this->ncpointer)->radius = this->radius;
  (this->ncpointer)->soft_radius = this->soft_radius;
  (this->ncpointer)->coord[0][0] = this->xc;
  (this->ncpointer)->coord[1][0] = this->yc;
  (this->ncpointer)->coord[2][0] = this->zc;
  (this->ncpointer)->shadow_cutoff = this->soft_radius*2.0;
  (this->ncpointer)->shadow_cutoffsq = pow((this->ncpointer)->shadow_cutoff,2);
}

void _NANOCARRIER :: write_NC_configuration(string filename){
    ofstream test1;
    test1.open(filename,ios::out);
    for (int i=0; i<this->num_ab; i++){
      test1<<this->x_antibody[i]<<" "<<this->y_antibody[i]<<" "<<this->z_antibody[i]<<endl;
    }
    test1.close();
}

void _NANOCARRIER :: shift_to_COM_frame(int nab, double *xpos, double *ypos, double *zpos){
  double xsum=0.0, ysum=0.0, zsum=0.0;
  for (int i=0; i<nab; i++){
    xsum += xpos[i];
    ysum += ypos[i];
    zsum += zpos[i];
  }
  
  xsum = xsum/(double) nab;
  ysum = ysum/(double) nab;
  zsum = zsum/(double) nab;
   
  for (int i=0; i<nab; i++){
    xpos[i] -= xsum;
    ypos[i] -= ysum;
    zpos[i] -= zsum;
  }  
}

void _NANOCARRIER :: set_base_tip_antibody(int i,double xpos, double ypos, double zpos){
    this->xf_antibody[i] = this->shift_PBC(this->xc + this->radius*xpos);                //base
    this->x_antibody[i] = this->shift_PBC(this->xc + this->soft_rad_ab[i]*xpos);         //tip
    this->yf_antibody[i] = this->shift_PBC(this->yc + this->radius*ypos);                //base
    this->y_antibody[i] = this->shift_PBC(this->yc + this->soft_rad_ab[i]*ypos);         //tip
    this->zf_antibody[i]= this->zc + this->radius*zpos;                                  //base
    this->z_antibody[i]= this->zc + this->soft_rad_ab[i]*zpos;                           //tip    
}

double _NANOCARRIER :: compute_max_soft_radius(){
  double maxrad = 0.0;
  for (int i=0; i<this->num_abtype; i++){
    double maxrad1 = this->radius + this->abclass.get_ablength(i);
    if (maxrad1 >maxrad)
      maxrad = maxrad1;
  }
  return maxrad;
}

void _NANOCARRIER :: shift_nc_position(double xpos, double ypos, double zpos){
  // shift the NC position from (xc,yc,zc) to (xpos,ypos,zpos)
  for (int i=0; i<this->num_ab; i++){    
     this->x_antibody[i] = this->shift_PBC(this->x_antibody[i]+xpos-this->xc);
     this->xf_antibody[i] = this->shift_PBC(this->xf_antibody[i]+xpos-this->xc);
     this->y_antibody[i] = this->shift_PBC(this->y_antibody[i]+ypos-this->yc);
     this->yf_antibody[i] = this->shift_PBC(this->yf_antibody[i]+ypos-this->yc);
     this->z_antibody[i] += zpos-this->zc;
     this->zf_antibody[i] += zpos-this->zc;
  }
  this->setrc(this->shift_PBC(xpos), this->shift_PBC(ypos), zpos);
}

double _NANOCARRIER :: shift_PBC(double rval){
    if (rval > this->L) rval -= this->L;
    if (rval < 0) rval += this->L;
    return rval;
}

double _NANOCARRIER :: apply_PBC(double rval){
    rval -= round(rval/this->L) * this->L;
    return rval;
}

void _NANOCARRIER :: store_configuration(){
	this->xcold = this->xc;
	this->ycold = this->yc;
	this->zcold = this->zc;
	this->anglexold = this->anglex;
	this->angleyold = this->angley;
	this->anglezold = this->anglez;
	for (int i=0;i<this->num_ab; i++){
		this->xtold[i]=this->x_antibody[i];
		this->ytold[i]=this->y_antibody[i];
		this->ztold[i]=this->z_antibody[i];
		this->xbold[i]=this->xf_antibody[i];
		this->ybold[i]=this->yf_antibody[i];
		this->zbold[i]=this->zf_antibody[i];
	}
}

void _NANOCARRIER :: restore_configuration(){
	this->xc = this->xcold;
	this->yc = this->ycold;
	this->zc = this->zcold;
	this->anglex = this->anglexold;
	this->angley = this->angleyold;
	this->anglez = this->anglezold;
	for (int i=0;i<this->num_ab; i++){
		this->x_antibody[i] = this->xtold[i];
		this->y_antibody[i] = this->ytold[i];
		this->z_antibody[i] = this->ztold[i];
		this->xf_antibody[i] = this->xbold[i];
		this->yf_antibody[i] = this->ybold[i];
		this->zf_antibody[i] = this->zbold[i];
	}
}

void _NANOCARRIER :: compute_ves_AB_distance(){
  for (int i=0;i<this->num_ab;i++) {                              
   double dx = this->apply_PBC(this->xc - this->x_antibody[i]);
   double dy = this->apply_PBC(this->yc - this->y_antibody[i]);
   double dz = this->zc - this->z_antibody[i];    
   cout<<str_green<<i<<" "<<pow(pow(dx,2)+pow(dy,2)+pow(dz,2),0.5)<<" "<<this->soft_rad_ab[i]<<str_black<<endl;
  }
}

double _NANOCARRIER :: compute_bond_length(int j, int k){
  extern _receptor re;
  double dx = this->apply_PBC(re.getxt(k) - this->x_antibody[j]);
  double dy = this->apply_PBC(re.getyt(k) - this->y_antibody[j]);
  double dz = re.getzt(k) - this->z_antibody[j];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

double _NANOCARRIER :: get_shadowr_dist(){
  return (this->ncpointer)->meandist;
}

double* _NANOCARRIER :: get_shadowr(){
	this->meanr[0] = this->ncpointer->meanr[0][0];
	this->meanr[1] = this->ncpointer->meanr[1][0];
	this->meanr[2] = this->ncpointer->meanr[2][0];
	return this->meanr;
}

char _NANOCARRIER :: get_bias_mode(){
      return this->bias_mode;
}

int _NANOCARRIER :: get_num_abtype(){
	return this->num_abtype;
}

int* _NANOCARRIER :: get_ab_types(){
  return this->ab_type;
}

int _NANOCARRIER :: get_ab_type(int k){
  return this->ab_type[k];
}

double _NANOCARRIER :: get_shadowcurv(){
	return (this->ncpointer)->meancurv;
}

int _NANOCARRIER :: get_nshadowvert(){
	return (this->ncpointer)->nshadow_vert;
}

void _NANOCARRIER :: reset_euler_angles(){
	this->anglex = this->anglez = 0.0;
	this->angley = PI/2.0;
}

void _NANOCARRIER :: setL (double L1){
	this->L = L1;
}

void _NANOCARRIER :: setH (double H1){
	this->H = H1;
}

void _NANOCARRIER :: setradius(double radius1, double soft_radius1) {
  this->radius = radius1;
  this->soft_radius = soft_radius1;
  (this->ncpointer)->radius = this->radius;
  (this->ncpointer)->soft_radius = this->soft_radius;
}

void _NANOCARRIER :: setrc(double xc1, double yc1, double zc1){
  this->xc = xc1;
  this->yc = yc1;
  this->zc = zc1;  
  (this->ncpointer)->coord[0][0] = this->xc;
  (this->ncpointer)->coord[1][0] = this->yc;
  (this->ncpointer)->coord[2][0] = this->zc;
}

void _NANOCARRIER :: setxc(double xc1) {
      this->xc = xc1;
      (this->ncpointer)->coord[0][0] = xc1;
}

void _NANOCARRIER :: setyc(double yc1) {
  this->yc = yc1;
  (this->ncpointer)->coord[1][0] = yc1;
}

void _NANOCARRIER :: setzc(double zc1) {
  this->zc = zc1;
  (this->ncpointer)->coord[2][0] = zc1;
}

void _NANOCARRIER :: setanglex(double ax1){
	this->anglex = ax1;
}

void _NANOCARRIER :: setangley(double ay1){
	this->angley = ay1;
}

void _NANOCARRIER :: setanglez(double az1){
	this->anglez = az1;
}

double _NANOCARRIER :: getradius() {
	return this->radius;
}

double _NANOCARRIER :: getxc() {
	return this->xc;
}

double _NANOCARRIER :: getyc() {
	return this->yc;
}

double _NANOCARRIER :: getzc() {
	return this->zc;
}

double _NANOCARRIER :: getanglex(){
	return this->anglex;
}

double _NANOCARRIER :: getangley(){
	return this->angley;
}

double _NANOCARRIER :: getanglez(){
        return this->anglez;
}

int _NANOCARRIER :: getnum_ab() {
	return this->num_ab;
}

double _NANOCARRIER :: getx_antibody(int k) {
	return this->x_antibody[k];
}

double _NANOCARRIER :: gety_antibody(int k) {
	return this->y_antibody[k];
}

double _NANOCARRIER :: getz_antibody(int k) {
	return this->z_antibody[k];
}

double _NANOCARRIER :: getxf_antibody(int k) {
	return this->xf_antibody[k];
}

double _NANOCARRIER :: getyf_antibody(int k) {
	return this->yf_antibody[k];
}

double _NANOCARRIER :: getzf_antibody(int k) {
	return this->zf_antibody[k];
}

double _NANOCARRIER :: getxf_npb_antibody(int k) {
	return this->xf_npb_antibody[k];
}

double _NANOCARRIER :: getyf_npb_antibody(int k) {
	return this->yf_npb_antibody[k];
}

double _NANOCARRIER :: getx_npb_antibody(int k) {
	return this->x_npb_antibody[k];
}

double _NANOCARRIER :: gety_npb_antibody(int k) {
	return this->y_npb_antibody[k];
}

double _NANOCARRIER :: getsoft_radius(){
	return this->soft_radius;
}

void _NANOCARRIER :: setxy_npb_antibody(){
	int i;
	double dx,dy,dz;
	cout.precision(dbl::digits10);
	
	for (i=0; i<this->num_ab; i++){
	  this->x_npb_antibody[i] = this->x_antibody[i] + round((this->xc-this->x_antibody[i])/this->L)*this->L;
	  this->xf_npb_antibody[i] = this->xf_antibody[i] + round((this->xc-this->xf_antibody[i])/this->L)*this->L;
	  this->y_npb_antibody[i] = this->y_antibody[i] + round((this->yc-this->y_antibody[i])/this->L)*this->L;
	  this->yf_npb_antibody[i] = this->yf_antibody[i] + round((this->yc-this->yf_antibody[i])/this->L)*this->L;
	
	  dx = this->xf_npb_antibody[i]-this->xc; 
	  dy = this->yf_npb_antibody[i]-this->yc;
	  dz = this->zf_antibody[i]-this->zc;
	  
	  double comp_radius = pow((dx*dx+dy*dy+dz*dz)/(this->radius*this->radius),2);
	  
	  // The input file is not a perfect sphere. So we allow for some amount of deviation from a sphere
	  if(comp_radius < 0.97 || comp_radius>1.03){                  
	       cout<<"Ratio of radius computed "<<comp_radius<<endl;
	       cout<<"Radius computed "<<dx*dx+dy*dy+dz*dz<<" "<<this->radius*this->radius<<endl;
	       cout<<"Shifting shows error for Antibody "<<i<<endl;
	       cout<<"PBC positions for NC center "<<this->xc<<" "<<this->yc<<" "<<this->zc<<endl;
	       cout<<"PBC position for antibody "<<this->x_antibody[i]<<" "<<this->y_antibody[i]<<" "<<this->z_antibody[i]<<endl;
	       cout<<"PBC position for antibody "<<this->xf_antibody[i]<<" "<<this->yf_antibody[i]<<" "<<this->zf_antibody[i]<<endl;
	       cout<<"Non PBC position for antibody "<<this->xf_npb_antibody[i]<<" "<<this->yf_npb_antibody[i]<<" "<<this->zf_antibody[i]<<endl;
	       cin.get();
	       }
	}
	return;
}

bool _NANOCARRIER :: getbonded(int k){
	return this->bonded[k];
}

void _NANOCARRIER :: setbonded(int k, bool temp){
  this->bonded[k] = temp;
}

int _NANOCARRIER :: getbondedto(int k){
	return this->bondedto[k];
}

void _NANOCARRIER :: setbondedto(int k,  int temp){
	this->bondedto[k] = temp;
}

void _NANOCARRIER :: translate (double dx, double dy, double dz) { 
  // shift the NC position from (xc,yc,zc) to (xc+dx,yc+dy,zc+dz)
  for (int i=0; i<this->num_ab; i++){    
     this->x_antibody[i] = this->shift_PBC(this->x_antibody[i]+dx);
     this->xf_antibody[i] = this->shift_PBC(this->xf_antibody[i]+dx);
     this->y_antibody[i] = this->shift_PBC(this->y_antibody[i]+dy);
     this->yf_antibody[i] = this->shift_PBC(this->yf_antibody[i]+dy);
     this->z_antibody[i] += dz; 
     this->zf_antibody[i] += dz;
  }
  
  this->setrc(this->shift_PBC(this->xc+dx), this->shift_PBC(this->yc+dy), this->zc+dz); 
} 

void _NANOCARRIER :: translate_OMP (double dx, double dy, double dz) { 
	// shift the NC position from (xc,yc,zc) to (xc+dx,yc+dy,zc+dz)
	
  #pragma omp parallel for
	for (int i=0; i<this->num_ab; i++){    
		this->x_antibody[i] = this->shift_PBC(this->x_antibody[i] + dx);
		this->xf_antibody[i] = this->shift_PBC(this->xf_antibody[i] + dx);
		this->y_antibody[i] = this->shift_PBC(this->y_antibody[i] + dy);
		this->yf_antibody[i] = this->shift_PBC(this->yf_antibody[i] + dy);
		this->z_antibody[i] += dz; 
		this->zf_antibody[i] += dz;
	}
	this->setrc(this->shift_PBC(this->xc+dx), this->shift_PBC(this->yc+dy), this->zc+dz);
} 

void _NANOCARRIER :: rotate (double phi, double theta, double psi,bool transpose) {   // For MD, use quaternion and not eulerian angles
	double cosphi, sinphi, costheta, sintheta, cospsi, sinpsi;
	double tempx, tempy, tempz, dx, dy, dz;
	double rot_mat[3][3];
	int i;

	this->anglex += phi;
	this->angley += theta;
	this->anglez += psi;
	
	if (this->anglex < 0) this->anglex += TWOPI;
	if (this->anglex >= TWOPI) this->anglex -= TWOPI;
	if (this->angley < 0) this->angley += PI;
	if (this->angley >= PI) this->angley -= PI;
	if (this->anglez < 0) this->anglez += TWOPI;
	if (this->anglez >= TWOPI) this->anglez -= TWOPI;
	
	//cout<<str_red<<omp_get_thread_num()<<" "<<this->index<<" "<<str_blue<<this->anglex<<" "<<this->angley<<" "<<this->anglez<<str_black<<endl;
	
	cosphi = cos(phi);
	cospsi = cos(psi);							 	 			  //temporary variables
	costheta = cos(theta);
	sinphi = sin(phi);
	sinpsi = sin(psi);
	sintheta = sin(theta);
	
	if (transpose){
		rot_mat[0][0] = cosphi*cospsi - sinphi*costheta*sinpsi;
		rot_mat[0][1] = sinphi*cospsi + cosphi*costheta*sinpsi;
		rot_mat[0][2] = sintheta*sinpsi;		                    //Rotation matrix (ZXZ rotation: see http://en.wikipedia.org/wiki/Euler_angles)
		rot_mat[1][0] = -cosphi*sinpsi - sinphi*costheta*cospsi;
		rot_mat[1][1] = cosphi*costheta*cospsi -sinphi*sinpsi;
		rot_mat[1][2] = sintheta*cospsi;
		rot_mat[2][0] = sinphi*sintheta;
		rot_mat[2][1] = -cosphi*sintheta;
		rot_mat[2][2] = costheta;
	}
	else{
		rot_mat[0][0] = cosphi*cospsi - sinphi*costheta*sinpsi;
		rot_mat[1][0] = sinphi*cospsi + cosphi*costheta*sinpsi;
		rot_mat[2][0] = sintheta*sinpsi;		                    //Rotation matrix (ZXZ rotation: see http://en.wikipedia.org/wiki/Euler_angles)
		rot_mat[0][1] = -cosphi*sinpsi - sinphi*costheta*cospsi;
		rot_mat[1][1] = cosphi*costheta*cospsi -sinphi*sinpsi;
		rot_mat[2][1] = sintheta*cospsi;
		rot_mat[0][2] = sinphi*sintheta;
		rot_mat[1][2] = -cosphi*sintheta;
		rot_mat[2][2] = costheta;
	}

	for (i=0;i<num_ab;i++){
		dx = this->x_antibody[i] - this->xc;
		dy = this->y_antibody[i] - this->yc;
		dz = this->z_antibody[i] - this->zc;
		if(dx > this->soft_rad_ab[i]) {
			if(abs(dx-this->soft_rad_ab[i]) > COMPARE_TOL){                // double comparison
		      		this->x_antibody[i] -= this->L;
					dx  = this->x_antibody[i]-this->xc;
			}
		}
		if(dx < -this->soft_rad_ab[i]){
			if(abs(dx+this->soft_rad_ab[i]) > COMPARE_TOL){               // double comparison 
			       	this->x_antibody[i] += this->L;
					dx = this->x_antibody[i]-this->xc;
			}
		}
		if(dy > this->soft_rad_ab[i]){
			if(abs(dy-this->soft_rad_ab[i]) > COMPARE_TOL){               // double comparison 
			      	this->y_antibody[i] -= this->L;
					dy = this->y_antibody[i]-this->yc;
			}
		}

		if(dy < -this->soft_rad_ab[i]){
			if(abs(dy+this->soft_rad_ab[i]) > COMPARE_TOL){
		       	this->y_antibody[i] += this->L;
				dy = this->y_antibody[i]-this->yc;
			}
		}
	
		tempx = rot_mat[0][0]*dx + rot_mat[0][1]*dy + rot_mat[0][2]*dz;
		tempy = rot_mat[1][0]*dx + rot_mat[1][1]*dy + rot_mat[1][2]*dz;
		tempz = rot_mat[2][0]*dx + rot_mat[2][1]*dy + rot_mat[2][2]*dz;
	
		this->x_antibody[i] = this->shift_PBC(tempx + this->xc);
		this->y_antibody[i] = this->shift_PBC(tempy + this->yc);
		this->z_antibody[i] = tempz + this->zc;
	}

	for (i=0;i<this->num_ab;i++){                                                                        		// to rotate fixed points
		dx = this->xf_antibody[i] - this->xc; 
		dy = this->yf_antibody[i] - this->yc; 
		dz = this->zf_antibody[i] - this->zc;

		if(dx > this->radius) {
			if(abs(dx-this->radius) > COMPARE_TOL){         
			    this->xf_antibody[i] -= this->L;
				dx=this->xf_antibody[i]-this->xc;
			}
		}
		if(dx < -this->radius){
			if(abs(dx+this->radius) > COMPARE_TOL){
				this->xf_antibody[i] += this->L;
				dx=this->xf_antibody[i]-this->xc;
			}
		}
		if(dy > this->radius){
			if(abs(dy-this->radius) > COMPARE_TOL){      
				this->yf_antibody[i] -= this->L;
				dy = this->yf_antibody[i]-this->yc;
			}
		}

		if(dy < -this->radius){
			if(abs(dy+this->radius) > COMPARE_TOL){                  
				this->yf_antibody[i] += this->L;
				dy = this->yf_antibody[i]-this->yc;
			}
		}

		tempx = rot_mat[0][0]*dx + rot_mat[0][1]*dy + rot_mat[0][2]*dz;
		tempy = rot_mat[1][0]*dx + rot_mat[1][1]*dy + rot_mat[1][2]*dz;
		tempz = rot_mat[2][0]*dx + rot_mat[2][1]*dy + rot_mat[2][2]*dz;

		this->xf_antibody[i] = tempx+this->xc;
		this->yf_antibody[i] = tempy+this->yc;
		this->zf_antibody[i] = tempz+this->zc;
 
		if (this->xf_antibody[i] < 0) this->xf_antibody[i] = this->xf_antibody[i] + this->L;
		if (this->yf_antibody[i] < 0) this->yf_antibody[i] = this->yf_antibody[i] + this->L;
		if (this->xf_antibody[i] > this->L) this->xf_antibody[i] = this->xf_antibody[i] - this->L;
		if (this->yf_antibody[i] > this->L) this->yf_antibody[i] = this->yf_antibody[i] - this->L;
	}
}

void _NANOCARRIER :: rotate_mp (double phi, double theta, double psi,bool transpose) {   // For MD, use quaternion and not eulerian angles
	double cosphi, sinphi, costheta, sintheta, cospsi, sinpsi;
	double rot_mat[3][3];

	this->anglex += phi;
	this->angley += theta;
	this->anglez += psi;
	
	if (this->anglex < 0) this->anglex += TWOPI;
	if (this->anglex >= TWOPI) this->anglex -= TWOPI;
	if (this->angley < 0) this->angley += PI;
	if (this->angley >= PI) this->angley -= PI;
	if (this->anglez < 0) this->anglez += TWOPI;
	if (this->anglez >= TWOPI) this->anglez -= TWOPI;
	
	cosphi = cos(phi);
	cospsi = cos(psi);							 	 			  //temporary variables
	costheta = cos(theta);
	sinphi = sin(phi);
	sinpsi = sin(psi);
	sintheta = sin(theta);
	
	if (transpose){
		rot_mat[0][0] = cosphi*cospsi - sinphi*costheta*sinpsi;
		rot_mat[0][1] = sinphi*cospsi + cosphi*costheta*sinpsi;
		rot_mat[0][2] = sintheta*sinpsi;		                    //Rotation matrix (ZXZ rotation: see http://en.wikipedia.org/wiki/Euler_angles)
		rot_mat[1][0] = -cosphi*sinpsi - sinphi*costheta*cospsi;
		rot_mat[1][1] = cosphi*costheta*cospsi -sinphi*sinpsi;
		rot_mat[1][2] = sintheta*cospsi;
		rot_mat[2][0] = sinphi*sintheta;
		rot_mat[2][1] = -cosphi*sintheta;
		rot_mat[2][2] = costheta;
	}
	else{
		rot_mat[0][0] = cosphi*cospsi - sinphi*costheta*sinpsi;
		rot_mat[1][0] = sinphi*cospsi + cosphi*costheta*sinpsi;
		rot_mat[2][0] = sintheta*sinpsi;		                    //Rotation matrix (ZXZ rotation: see http://en.wikipedia.org/wiki/Euler_angles)
		rot_mat[0][1] = -cosphi*sinpsi - sinphi*costheta*cospsi;
		rot_mat[1][1] = cosphi*costheta*cospsi -sinphi*sinpsi;
		rot_mat[2][1] = sintheta*cospsi;
		rot_mat[0][2] = sinphi*sintheta;
		rot_mat[1][2] = -cosphi*sintheta;
		rot_mat[2][2] = costheta;
	}

	#pragma omp for
	for (int i=0;i<num_ab;i++){
		double dx, dy, dz, tempx, tempy, tempz;
		double dx1, dy1, dz1, tempx1, tempy1, tempz1;
		dx = this->x_antibody[i] - this->xc;
		dy = this->y_antibody[i] - this->yc;
		dz = this->z_antibody[i] - this->zc;
		if(dx > this->soft_rad_ab[i]){
			if(abs(dx-this->soft_rad_ab[i]) > COMPARE_TOL){                // double comparison
		      	this->x_antibody[i] -= this->L;
				dx  = this->x_antibody[i]-this->xc;
			}
		}
		if(dx < -this->soft_rad_ab[i]){
			if(abs(dx+this->soft_rad_ab[i]) > COMPARE_TOL){               // double comparison 
			       	this->x_antibody[i] += this->L;
				dx = this->x_antibody[i]-this->xc;
			}
		}
		if(dy > this->soft_rad_ab[i]){
			if(abs(dy-this->soft_rad_ab[i]) > COMPARE_TOL){               // double comparison 
			      	this->y_antibody[i] -= this->L;
				dy = this->y_antibody[i]-this->yc;
			}
		}

		if(dy < -this->soft_rad_ab[i]){
			if(abs(dy+this->soft_rad_ab[i]) > COMPARE_TOL){
		       	this->y_antibody[i] += this->L;
			dy = this->y_antibody[i]-this->yc;
			}
		}
	
		tempx = rot_mat[0][0]*dx + rot_mat[0][1]*dy + rot_mat[0][2]*dz;
		tempy = rot_mat[1][0]*dx + rot_mat[1][1]*dy + rot_mat[1][2]*dz;
		tempz = rot_mat[2][0]*dx + rot_mat[2][1]*dy + rot_mat[2][2]*dz;
	
		this->x_antibody[i] = tempx + xc;
		this->y_antibody[i] = tempy + yc;
		this->z_antibody[i] = tempz + zc;
	
		if (this->x_antibody[i] < 0) this->x_antibody[i] = this->x_antibody[i] + this->L;
		if (this->y_antibody[i] < 0) this->y_antibody[i] = this->y_antibody[i] + this->L;
		if (this->x_antibody[i] > this->L) this->x_antibody[i] = this->x_antibody[i] - this->L;
		if (this->y_antibody[i] > this->L) this->y_antibody[i] = this->y_antibody[i] - this->L;

		dx1 = this->xf_antibody[i]-this->xc; 
		dy1 = this->yf_antibody[i]-this->yc; 
		dz1 = this->zf_antibody[i]-this->zc;

		if(dx1 > this->radius) {
			if(abs(dx1-this->radius) > COMPARE_TOL){         
			    this->xf_antibody[i] -= this->L;
				dx1 = this->xf_antibody[i]-this->xc;
			}
		}
		
		if(dx1 < -this->radius){
			if(abs(dx1+this->radius) > COMPARE_TOL){
				this->xf_antibody[i] += this->L;
				dx = this->xf_antibody[i]-this->xc;
			}
		}
		
		if(dy1 > this->radius){
			if(abs(dy1-this->radius) > COMPARE_TOL){      
				this->yf_antibody[i] -= this->L;
				dy1 = this->yf_antibody[i]-this->yc;
			}
		}

		if(dy1 < -this->radius){
			if(abs(dy1 + this->radius) > COMPARE_TOL){                  
				this->yf_antibody[i] += this->L;
				dy1 = this->yf_antibody[i]-this->yc;
			}
		}

		tempx1 = rot_mat[0][0]*dx1 + rot_mat[0][1]*dy1 + rot_mat[0][2]*dz1;
		tempy1 = rot_mat[1][0]*dx1 + rot_mat[1][1]*dy1 + rot_mat[1][2]*dz1;
		tempz1 = rot_mat[2][0]*dx1 + rot_mat[2][1]*dy1 + rot_mat[2][2]*dz1;

		this->xf_antibody[i] = tempx1 + this->xc;
		this->yf_antibody[i] = tempy1 + this->yc;
		this->zf_antibody[i] = tempz1 + this->zc;
 
		if (this->xf_antibody[i] < 0) this->xf_antibody[i] = this->xf_antibody[i] + this->L;
		if (this->yf_antibody[i] < 0) this->yf_antibody[i] = this->yf_antibody[i] + this->L;
		if (this->xf_antibody[i] > this->L) this->xf_antibody[i] = this->xf_antibody[i] - this->L;
		if (this->yf_antibody[i] > this->L) this->yf_antibody[i] = this->yf_antibody[i] - this->L;
	}
	#pragma omp barrier
}

bool _NANOCARRIER :: does_bond_breaks(double tempx, double tempy, double tempz){
   bool bond_break = false;
   extern _receptor re;
   extern _intparam intparam;
   for (int ab=0; ab<this->num_ab; ab++){
    if (this->bonded[ab]){
     int ant_no = this->bondedto[ab];
     double max_AB_ANT_distance = intparam.get_delrsq(re.gettype(ant_no), this->ab_type[ab]);
     double dx = this->x_antibody[ab]+tempx - re.getantigenxt(ant_no);   
     double dy = this->y_antibody[ab]+tempy - re.getantigenyt(ant_no);
     double dz = this->z_antibody[ab]+tempz - re.getantigenzt(ant_no);
     dx = dx-round(dx/this->L)*this->L; 
	 dy=dy - round(dy/this->L)*this->L;
     double dist = dx*dx + dy*dy + dz*dz;
     if (dist > max_AB_ANT_distance){
      bond_break=true;                          // if the movement stretches the bond beyond the reaction distance reject the move
      return bond_break;
     }
    }
   }
   return bond_break;
}

bool _NANOCARRIER :: does_bond_breaks(){
   bool bond_break = false;
   extern _receptor re;
   extern _intparam intparam;
   for (int ab=0; ab<this->num_ab;ab++){
    if (this->bonded[ab]){
     int ant_no = this->bondedto[ab];
     double max_AB_ANT_distance = 
			intparam.get_delrsq(re.gettype(ant_no), this->ab_type[ab]);
     double dx = this->x_antibody[ab] - re.getantigenxt(ant_no);                 // check for the distance between the antigen tip and antibody
     double dy = this->y_antibody[ab] - re.getantigenyt(ant_no);
     double dz = this->z_antibody[ab] - re.getantigenzt(ant_no);
     dx = dx - round(dx/this->L) * this->L; 
     dy = dy - round(dy/this->L) * this->L;
     double dist=dx*dx+dy*dy+dz*dz;
     
     if (dist > max_AB_ANT_distance){
     bond_break=true;                                        // if the movement stretches the bond beyond the reaction distance reject the move
     return bond_break;
     }
    }
   }
   return bond_break;
}


bool _NANOCARRIER :: Does_Bond_Breaks(){
	bool bond_break = false;
	extern _receptor re;
	extern _intparam intparam;
	for (int ab=0; ab<this->num_ab;ab++){
		if (this->bonded[ab]){
			int ant_no = this->bondedto[ab];
			double max_AB_ANT_distance = intparam.get_delrsq(re.gettype(ant_no), this->ab_type[ab]);
			double dx = this->x_antibody[ab] - re.getantigenxt(ant_no);                 // check for the distance between the antigen tip and antibody
			double dy = this->y_antibody[ab] - re.getantigenyt(ant_no);
			double dz = this->z_antibody[ab] - re.getantigenzt(ant_no);
			dx = dx - round(dx/this->L) * this->L; 
			dy = dy - round(dy/this->L) * this->L;
			double dist=dx*dx+dy*dy+dz*dz;
			if (dist > max_AB_ANT_distance){
				bond_break=true;                                        // if the movement stretches the bond beyond the reaction distance reject the move
				cout<<this->index<<" "<<ab<<" "<<ant_no<<" "<<dist<<" "<<max_AB_ANT_distance<<endl;
				return bond_break;
			}
		}
	}
	return bond_break;
}

double _NANOCARRIER :: get_rstart(){
  return this->rstart;
}

double _NANOCARRIER :: get_hstart(){
  return this->hstart;
}

double _NANOCARRIER :: get_kbias(){
  return this->kbias;
}

double _NANOCARRIER :: get_biasref(){
  return this->biasref;
}

double _NANOCARRIER :: get_lambda(){
  return this->lambda;
}

string _NANOCARRIER :: get_biasdir(){
  return this->biasdir;
}

double _NANOCARRIER :: get_biasbinsize(){
  return this->binsize;
}

double* _NANOCARRIER :: get_biasparameters(){
  this->biasparm[0] = this->kbias;
  this->biasparm[1] = this->biasref;
  this->biasparm[2]= this->binsize;
  this->biasparm[3] = this->lambda;
  return this->biasparm;
}

void _NANOCARRIER :: set_biasref(double refval){
  this->biasref = refval;
  (this->ncpointer)->biasref = this->biasref;
  
}

void _NANOCARRIER :: shift_biaswindow(){
  if (this->biasdir.compare("approach") == 0)
    this->biasref -= this->binsize;
  else if (this->biasdir.compare("recede") == 0)
    this->biasref += this->binsize;
}

//@ Functions for counters to handle MC move statistics

void  _NANOCARRIER :: reset_counters(string cname){
 if (cname.compare("trans") == 0){
  fill_n(this->trans_counter,4,0);
 }
 else if (cname.compare("rot") == 0){
  fill_n(this->rot_counter,4,0);
 }
 else{
   cout <<"Invalid entry for reset counters in nanocarrier.cpp"<<endl;
   cout <<"Options are 'trans'/'rot'";
   exit(1);
 }
}

void  _NANOCARRIER :: increment_counters(string cname, bool acceptance){
  if (cname.compare("trans") == 0){
    if (this->multivalency == 0){
      this->trans_counter[1] += 1;
      if (acceptance) this->trans_counter[0] += 1;
      if (this->trans_counter[1] == this->adjust_interval) this->reset_step_size(1);
    }
    
    else{
      this->trans_counter[3] += 1;
      if (acceptance) this->trans_counter[2] += 1;
      if (this->trans_counter[3] == this->adjust_interval) this->reset_step_size(2);
    }  
  }
  
  else if (cname.compare("rot") == 0){
    if (this->multivalency == 0){
      this->rot_counter[1] += 1;
      if (acceptance) this->rot_counter[0] += 1;
      if (this->rot_counter[1] == this->adjust_interval) this->reset_step_size(3);
    }
    
    else{
      this->rot_counter[3] += 1;
      if (acceptance) this->rot_counter[2] += 1;
      if (this->rot_counter[3] == this->adjust_interval) this->reset_step_size(4);
    }  
  }
  
  else{
   cout <<"Invalid entry for increment counters in nanocarrier.cpp"<<endl;
   cout <<"Options are 'trans'/'rot'";
   exit(1);
 }
}

int* _NANOCARRIER :: get_counter(string cname){
  if (cname.compare("trans")==0)
    return this->trans_counter;
  else if (cname.compare("rot")==0)
    return this->rot_counter;
  else{
   cout <<"Invalid entry for reset counters in nanocarrier.cpp"<<endl;
   cout <<"Options are 'trans'/'rot'";
   exit(1);
 }
}

void _NANOCARRIER :: reset_step_size(int mode){ 
  switch (mode){
    case 1:{            // trans/unbonded
      double dratioub = (double) this->trans_counter[0]/(double) this->trans_counter[1];
      if (dratioub < 0.5)
       tstep_ub = max(0.95*tstep_ub,this->radius/10000.0);
      else
       tstep_ub = min(1.05*tstep_ub, this->radius/100.0);
      this->trans_counter[0] = this->trans_counter[1] = 0;
      break;
    }
  
    case 2:{           // trans/bonded
      double dratiob = (double) this->trans_counter[2]/(double) this->trans_counter[3];   
      if (dratiob <0.5)
       tstep_b = max(0.95*tstep_b,this->radius/10000.0);
      else
       tstep_b = min(1.05*tstep_b,this->radius/100.0);
      this->trans_counter[2] = this->trans_counter[3] = 0;
      break;
    }
  
    case 3:{             // rot/unbonded
     double dratiorub = (double) this->rot_counter[0]/(double) this->rot_counter[1];
     if (dratiorub < 0.5)
      rstep_ub = max(0.95*rstep_ub,PI/20000.0);
     else 
      rstep_ub = min(1.05*rstep_ub,PI/100.0);
     this->rot_counter[0] = this->rot_counter[1] = 0;
     break;
    }
  
    case 4:{           // rot/bonded
     double dratiorb = (double) this->rot_counter[2]/(double) this->rot_counter[3]; 
     if (dratiorb <0.5)
       rstep_b = max(0.95*rstep_b,PI/20000.0);
     else
       rstep_b = min(1.05*rstep_b,PI/100.0);
     this->rot_counter[2] = this->rot_counter[3] = 0;
     break;
    }
  
    default:{
     cout <<"Invalid entry for reset counters in nanocarrier.cpp"<<endl;
     cout <<"Options are 'trans'/'rot'";
     exit(1);
     break;
    }
  } 
}

double _NANOCARRIER :: get_step_size(string cname){
  double stepsize = 0.0;
  if (cname.compare("trans") == 0){
    if (this->multivalency == 0){
      stepsize = tstep_ub;
    }
    else{
      stepsize = tstep_b;
    }
  }
  
  else if (cname.compare("rot") == 0){
    if (this->multivalency == 0){
      stepsize = rstep_ub;
    }
    else{
      stepsize = rstep_b;
    }
  }
  
  else{
   cout <<"Invalid entry for reset counters in nanocarrier.cpp"<<endl;
   cout <<"Options are 'trans'/'rot'";
   exit(1);
  }
  
  return stepsize;
}

double* _NANOCARRIER :: get_mc_stepsizes(){
  this->stepsize[0] = this->tstep_ub;
  this->stepsize[1] = this->tstep_b;
  this->stepsize[2] = this->rstep_ub;  
  this->stepsize[3] = this->rstep_b;
  return this->stepsize;
}

//@ Functions to handle multivalency

int _NANOCARRIER :: get_multivalency(){
  return this->multivalency;
}

void _NANOCARRIER :: set_multivalency(int mval){
  this->multivalency = mval;
  (this->ncpointer)->multivalency = this->multivalency;
}

void _NANOCARRIER :: update_bond_counters(int dmval, int dnbon, int dnbr){
  this->multivalency += dmval;
  (this->ncpointer)->multivalency = this->multivalency;
  this->nbonded += dnbon;
  this->nbroken += dnbr;
}

int* _NANOCARRIER :: get_bond_counters(){
  bond_counters[0] = this->multivalency;
  bond_counters[1] = this->nbonded;
  bond_counters[2] = this->nbroken;
  return bond_counters;
}

void _NANOCARRIER :: reset_bond_counters(){
  this->nbonded = 0;
  this->nbroken =0;
}

void _nanocarrier :: init() {
  set_array_sizes();	
  this->v = new (nothrow) _NANOCARRIER[num_vesicles];
  for (int i=0; i<this->num_vesicles;i++) 
     this->v[i].init(i);
}

void _nanocarrier :: init(ifstream &ifile){
  set_array_sizes();	
  this->v = new (nothrow) _NANOCARRIER[num_vesicles];
  if (this->v == 0) {
	  cout<<"\n Memory could not be allocated "<<__PRETTY_FUNCTION__<<endl;
	  exit(1);
  }
  for (int i=0; i<this->num_vesicles;i++){
     this->v[i].init(i,ifile);
	 ifile.clear();
	 ifile.seekg(0,ios::beg);
  }
}

void _nanocarrier :: set_array_sizes(){
	extern _datainput data;
	this->L = data.getperiodic_box_length();
	this->H = data.getperiodic_box_height();
	this->temperature = data.gettemperature();
	this->kBT = (kb*temperature);
	this->beta = 1.0/this->kBT;
	this->num_vesicles = data.getnum_members('v');            // number of vesicles
	this->num_antigens = data.getnum_members('c');
	
	this->thetat = new(nothrow) double[this->num_samp];
	this->phit = new(nothrow) double[this->num_samp];
	if ((this->thetat == 0) || (this->phit == 0)) cout<<"\n Memory could not be allocated - Pos 1";
	
	this->xpn = new(nothrow) double[this->num_samp];
	this->ypn = new(nothrow) double[this->num_samp];
	this->zpn = new(nothrow) double[this->num_samp];
	if ((this->xpn == 0) || (this->ypn == 0) || (this->zpn == 0)) cout<<"\n Memory could not be allocated - Pos 2";
	
	this->b = new(nothrow) int[this->num_samp];
	this->weight = new(nothrow) double[this->num_samp];
	this->Eo = new(nothrow) double[this->num_samp];
	
	//Pointers to Fortran data
	this->verptr = &__module_datastruct_MOD_ver;
	this->triptr = &__module_datastruct_MOD_tri;
	this->antptr = &__module_datastruct_MOD_antig;
}

_nanocarrier :: _nanocarrier(double L1, double H1, int num_vesicles1){
  this->L = L1;
  this->H = H1;
  this->num_vesicles = num_vesicles1;
  v = new (nothrow) _NANOCARRIER[this->num_vesicles];
  if (v == 0)
    cout<<"\nError in System: memory could not be allocated";
}

void _nanocarrier :: reset_box_dimensions(){
	extern _datainput data;
	this->L = data.getperiodic_box_length();
	this->H = data.getperiodic_box_height();
	for (int i=0; i<this->num_vesicles;i++) 
		v[i].reset_box_dimensions();
}

void _nanocarrier :: dump_conf(int i, ofstream &ofile){
	v[i].dump_conf(ofile);
}

_NANOCARRIER _nanocarrier :: getvesicle(int j){
  return v[j];
}

void _nanocarrier :: print_address(int i){
  v[i].print_address();
}

int _nanocarrier :: getnum_members(){
  return num_vesicles;
}

double _nanocarrier :: getL(){
  return L;
}

double _nanocarrier :: getH(){
  return H;
}

void _nanocarrier :: setvesicle (_NANOCARRIER v_temp, int i){
  v[i] = v_temp;
}

void _nanocarrier :: write_NC_configuration(int  i, string filename){
  v[i].write_NC_configuration(filename);
}

//@ NC center of mass
void _nanocarrier :: store_configuration(int i){
	this->v[i].store_configuration();
}

void _nanocarrier :: restore_configuration(int i){
	this->v[i].restore_configuration();
}

void _nanocarrier :: setxc(int i, double x){
  v[i].setxc(x);
}

void _nanocarrier :: setyc(int i, double y){
  v[i].setyc(y);
}

void _nanocarrier :: setzc(int i, double z){
  v[i].setzc(z);
}

double _nanocarrier :: getxc(int i){
  return v[i].getxc();
}

double _nanocarrier :: getyc(int i){
  return v[i].getyc();
}

double _nanocarrier :: getzc(int i){
  return v[i].getzc();
}

//@ NC angles
void _nanocarrier :: setanglex(int i, double x){
  v[i].setanglex(x);
}

void _nanocarrier :: setangley(int i, double y){
  v[i].setangley(y);
}

void _nanocarrier :: setanglez(int i, double z){
  v[i].setanglez(z);
}

double _nanocarrier :: getanglex(int i){
  return v[i].getanglex();
}

double _nanocarrier :: getangley(int i){
  return v[i].getangley();
}

double _nanocarrier :: getanglez(int i){
  return v[i].getanglez();
}

void _nanocarrier :: reset_euler_angles(int i){
  v[i].reset_euler_angles();
}

//@ NC dimensions

void _nanocarrier :: setradius(int i, double radius, double soft_radius1){
  v[i].setradius(radius, soft_radius1);
}

double _nanocarrier :: getradius(int i){
  return v[i].getradius();
}

double _nanocarrier :: getsoft_radius(int i){
  return v[i].getsoft_radius();
}

void _nanocarrier :: setnum_ab(int i, int n){
  v[i].setnum_ab();
}

int _nanocarrier :: getnum_ab(int i){
  return v[i].getnum_ab();
}

int _nanocarrier :: get_num_abtype(int i){
  return v[i].get_num_abtype();
}

int* _nanocarrier :: get_ab_types(int i){
  return v[i].get_ab_types();
}

int _nanocarrier :: get_ab_type(int i, int k){
  return v[i].get_ab_type(k);
}

//@ Solid body motions

void _nanocarrier :: shift_nc_position(int i, double xpos, double ypos, double zpos){
  v[i].shift_nc_position(xpos,ypos,zpos);
}

void _nanocarrier :: translate(int i, double tempx, double tempy, double tempz){
  v[i].translate(tempx, tempy, tempz);
}

void _nanocarrier :: rotate(int i, double theta, double psi, double phi,bool transpose){
  v[i].rotate(theta, psi, phi,transpose);
}

//@ Antibody base and tip positions

void _nanocarrier :: setxy_npb_antibody(int ves_num){
  v[ves_num].setxy_npb_antibody();
}

double _nanocarrier :: getx_antibody(int i, int k) {
  return v[i].getx_antibody(k);
}

double _nanocarrier :: gety_antibody(int i, int k) {
  return v[i].gety_antibody(k);
}

double _nanocarrier :: getz_antibody(int i, int k) {
  return v[i].getz_antibody(k);
}

double _nanocarrier :: getxf_antibody(int i, int k) {
  return v[i].getxf_antibody(k);
}

double _nanocarrier :: getyf_antibody(int i, int k) {
  return v[i].getyf_antibody(k);
}

double _nanocarrier :: getzf_antibody(int i, int k) {
  return v[i].getzf_antibody(k);
}

double _nanocarrier :: getx_npb_antibody(int i, int k) {
  return v[i].getx_npb_antibody(k);
}

double _nanocarrier :: gety_npb_antibody(int i, int k) {
  return v[i].gety_npb_antibody(k);
}

double _nanocarrier :: getxf_npb_antibody(int i, int k) {
  return v[i].getxf_npb_antibody(k);
}

double _nanocarrier :: getyf_npb_antibody(int i, int k) {
  return v[i].getyf_npb_antibody(k);
}


//@ bonded details
bool _nanocarrier :: getbonded(int i, int k){
  return v[i].getbonded(k);
}

bool _nanocarrier :: getbonded(int i){                    
  int vnum_ab = getnum_ab(i);
  for (int j=0;j<vnum_ab;j++) {
    if (getbonded(i,j)==1) return 1;
  }
  return 0;
}

void _nanocarrier :: setbonded(int i, int k, bool temp){
  v[i].setbonded(k, temp);
}

int _nanocarrier :: getbondedto(int i, int k){
  return v[i].getbondedto(k);
}

void _nanocarrier :: setbondedto(int i, int k,  int temp){
  v[i].setbondedto(k, temp);
}

bool _nanocarrier :: does_bond_breaks(int i,double tempx,double tempy,double tempz){
  return v[i].does_bond_breaks(tempx,tempy,tempz);
}

bool _nanocarrier :: does_bond_breaks(int i){
  return v[i].does_bond_breaks();
}

bool _nanocarrier :: Does_Bond_Breaks(int i){
	return v[i].Does_Bond_Breaks();
}

void _nanocarrier :: updatebonds(int i){                             //To update the bonds
  int k;
  int tem;
  double dist, rrx, rry, rrz;
  extern _receptor re;
  extern _intparam intparam;

  for (k=0; k<this->getnum_ab(i); k++) {
    if (this->getbonded(i,k) == 1){
      tem = this->getbondedto(i,k);
      rrx = this->getx_antibody(i,k) - re.getantigenxt(tem);
      rry = this->gety_antibody(i,k) - re.getantigenyt(tem);                        //calculate the real distance between antibody and antigen tip:Jin
      rrz = this->getz_antibody(i,k) - re.getantigenzt(tem);
      rrx = rrx - this->L*round(rrx/this->L);
      rry = rry - this->L*round(rry/this->L);
      dist = rrx*rrx + rry*rry + rrz*rrz;
      double Dreaction = intparam.get_delrsq(i, k, tem);
      if (dist > Dreaction){                                        // if (dist > chem_cutoffu*chem_cutoffu || dist<chem_cutoffl*chem_cutoffl)
       this->setbonded(i,k, 0);                                     // antibody k on vesicle i is unbonded  
       re.setbonded(tem,0);                                         // antigen tem is also unbonded
       this->setbondedto(i,k,-1);
       re.reset_antigen_bond(tem);
       re.setantigentip_position(tem);                             // automatically sets the antigen tip to be along the vertex normal  
       re.set_theta_phi(tem,0.0,0.0);
      }
    }
  } 
}

void _nanocarrier :: reset_antab_bonds(int i){
  for (int k=0; k<this->getnum_ab(i); k++) 
       this->setbonded(i,k, 0);                                     // antibody k on vesicle i is unbonded 
}

//@ distance calculations

void _nanocarrier :: compute_ves_AB_distance(int i){
  v[i].compute_ves_AB_distance();
}

double _nanocarrier :: compute_bond_length(int i, int j, int k){
  return v[i].compute_bond_length(j,k);
}

//@ deals with the NC shadow on the membrane

double _nanocarrier :: get_shadowr_dist(int i){
  return v[i].get_shadowr_dist();
}

double* _nanocarrier :: get_shadowr(int i){
  return v[i].get_shadowr();
}

double _nanocarrier :: get_shadowcurv(int i){
  return v[i].get_shadowcurv();
}

int _nanocarrier :: get_nshadowvert(int i){
  return v[i].get_nshadowvert();
}

void _nanocarrier :: store_restore_shadow(int flag, int ncnum){
  __module_datastruct_MOD_store_restore_shadowvertices(&flag,&ncnum);
}

void _nanocarrier :: compute_shadow(int ncnum){
  __module_datastruct_MOD_compute_nc_shadow_vertices(&ncnum);
}

//@ handle all bias parameters
char _nanocarrier:: get_bias_mode(int i){
  return v[i].get_bias_mode();
}

double _nanocarrier :: get_rstart(int i){
  return v[i].get_rstart();
}

double _nanocarrier :: get_hstart(int i){
  return v[i].get_hstart();
}

double _nanocarrier :: get_kbias(int i){
  return v[i].get_kbias();
}

double _nanocarrier :: get_biasref(int i){
  return v[i].get_biasref();
}

double _nanocarrier :: get_lambda(int i){
  return v[i].get_lambda();
}

string _nanocarrier :: get_biasdir(int i){
  return v[i].get_biasdir();
}

double _nanocarrier :: get_biasbinsize(int i){
  return v[i].get_biasbinsize();
}

double* _nanocarrier :: get_biasparameters(int i){
  return v[i].get_biasparameters();
}

void _nanocarrier :: set_biasref(int i,double refval){
  v[i].set_biasref(refval);
}

void _nanocarrier :: shift_biaswindow(){
  for (int i=0; i<this->num_vesicles; i++){
    if ( (this->get_bias_mode(i) =='Z') || (this->get_bias_mode(i) =='H') ){
      cout<<"Before -- NC "<<i<<this->get_bias_mode(i)<<" "<<this->get_biasref(i)<<endl;
      v[i].shift_biaswindow();
      cout<<"After -- NC "<<i<<this->get_bias_mode(i)<<" "<<this->get_biasref(i)<<endl;
    }
  }
}

//@ Function to handle NC counters

void _nanocarrier :: reset_counters(int i, string cname, bool allflag){
  if (allflag){
    for(int j=0; j<this->num_vesicles; j++)
      v[j].reset_counters(cname);
  }
  else{
    v[i].reset_counters(cname);
  }
}

void _nanocarrier :: increment_counters(int i, string cname, bool acceptance){
  v[i].increment_counters(cname, acceptance);
}

int* _nanocarrier :: get_counter(int i, string cname){
  return v[i].get_counter(cname);
}

double _nanocarrier :: get_step_size(int i, string cname){
  return v[i].get_step_size(cname);
}

double* _nanocarrier :: get_mc_stepsizes(int i){
  return v[i].get_mc_stepsizes();
}

//@ Functions to handle multivalency

int _nanocarrier :: get_multivalency(int i){
  return v[i].get_multivalency();
}

void _nanocarrier :: set_multivalency(int i, int mval){
  v[i].set_multivalency(mval);
}

void _nanocarrier :: update_bond_counters(int i, int dmval, int dnbon, int dnbr){
  v[i].update_bond_counters(dmval, dnbon, dnbr);
}

int* _nanocarrier :: get_bond_counters(int i){
  return v[i].get_bond_counters();
}

void _nanocarrier :: reset_bond_counters(int i){
  v[i].reset_bond_counters();
}

//@ Function to handles various Monte Carlo moves

void _nanocarrier :: translation(){
  extern _potential p;
  extern _receptor re;
  extern _linklist link1;
  extern MTRand mtrand1;
  extern _selfavoidance sa;

  double  *rmean_old, *rmean_new, bias_ener,meanH_old,dE;
  int flag,sel_member;
  bool bond_break = false, accpflag = false;
  
  bias_ener = meanH_old = 0.0;

  if (this->num_vesicles > 1)                                               //  Choose a vesicle randomly if there are multiple NC
    sel_member = mtrand1.randInt(num_vesicles-1);
  else
    sel_member = 0;

  if (sel_member > this->num_vesicles){                                      // Check just for safety
    cout <<"\n Fatal Error in selecting vesicles"; 
    exit(1); 
  }

  double stepsize = this->get_step_size(sel_member,"trans");                 // step size is different for bonded and unbonded vesicles
  double tempx = mtrand1.rand(-1.0,1.0)*stepsize;                            // random disp. of max length dl along x          
  double tempy = mtrand1.rand(-1.0,1.0)*stepsize;                            // random disp. of max length dl along y
  double tempz = mtrand1.rand(-1.0,1.0)*stepsize;                            // disp along z direction
  
  
  if(this->getzc(sel_member) + tempz + this->getsoft_radius(sel_member) >= this->H) tempz=0.0;  // NO PBC along z

  if (this->getbonded(sel_member))
    bond_break = this->does_bond_breaks(sel_member,tempx,tempy,tempz);         // check if a  bond breaks after movement
	
    
  if (bond_break == false){                                                     // perform the move only if the bond remains intact 
	 if(!sa.overlap(sel_member,'v','m',tempx,tempy,tempz) && 
	   !sa.overlap(sel_member,'v','c',tempx,tempy,tempz) &&
	   !sa.overlap(sel_member,'v','v',tempx,tempy,tempz)){                         // after self avoidance is valid
      double energy_old = p.calc(sel_member,'v');                                  // energy for old configuration ( Antigen-AB + glycocalyx)   
      double rold[3] = {this->getxc(sel_member),this->getyc(sel_member),this->getzc(sel_member)};
      rmean_old = this->get_shadowr(sel_member);
      meanH_old = this->get_shadowcurv(sel_member);
	  this->store_configuration(sel_member);
      this->translate(sel_member, tempx, tempy, tempz);                         
      double rnew[3] = {this->getxc(sel_member),this->getyc(sel_member),this->getzc(sel_member)};
      flag=0;
      this->store_restore_shadow(flag,sel_member+1);
      this->compute_shadow(sel_member+1);
   
      //@ if a bias has been applied
      char bmode = this->get_bias_mode(sel_member);
      if (bmode == 'H')                                                        // bias the height under the nanocarrier  
        bias_ener = p.compute_biasing_potential(sel_member,meanH_old,this->get_shadowcurv(sel_member));

      if ((bmode == 'Z') || (bmode=='T')){                                     // bias the height under the nanocarrier   
        rmean_new = this->get_shadowr(sel_member);
        bias_ener = p.compute_biasing_potential(sel_member, rold, rmean_old, rnew, rmean_new); 
      }

      //@ final energy
      double energy_new = p.calc(sel_member,'v');                               // Total Antigen-Antibody energy + glycocalyx interaction
      dE = energy_new - energy_old + bias_ener;
      //@ Metropolis Scheme to accept or reject moves
      if ((dE > 0.0) && (mtrand1.rand() > exp(-this->beta*dE))){               // move rejected
        //this->translate(sel_member, -tempx, -tempy, -tempz);                     // translate back to original position
		this->restore_configuration(sel_member);
        flag=1;
        this->store_restore_shadow(flag, sel_member+1);
      } 
      else{                                                                                                // move accepted
        link1.calc(sel_member,this->getxc(sel_member),this->getyc(sel_member),this->getzc(sel_member),'v'); //update the link cell
        accpflag = true;
      }
    }
  }
  // log the results in the counter
  this->increment_counters(sel_member,"trans",accpflag);
}

void _nanocarrier :: rotation(){
	extern MTRand mtrand1;
	extern _potential p;
	int sel_member;
	bool accpflag = false;

	//@
	if (this->num_vesicles >= 2)
		sel_member = mtrand1.randInt(this->num_vesicles-1);
	else
		sel_member = 0;
  
	this->store_configuration(sel_member);

	//@
	double energy_old = p.calc(sel_member,'v');                                                           // energy for old configuration
	double angle_max = this->get_step_size(sel_member,"rot");                                             // rotational magnitude
    
	//random increment of the Euler angles
	double tempx = mtrand1.rand(-1.0,1.0)*angle_max;
	double tempy = mtrand1.rand(-1.0,1.0)*angle_max;
	double tempz = mtrand1.rand(-1.0,1.0)*angle_max;

	this->rotate(sel_member, tempx, tempy, tempz,false);
	bool bond_breaks = this->does_bond_breaks(sel_member);

	if(!bond_breaks){
		double energy_new = p.calc(sel_member,'v');
		if (mtrand1.rand() > exp(-this->beta*(energy_new-energy_old))) {
			this->restore_configuration(sel_member);
		}
		else {                                          //we break the old bonds only if the bead was initially within cutoff distance
			accpflag = true;
		}
	}
  
	else {
		this->restore_configuration(sel_member);
	}
  
	this->increment_counters(sel_member,"rot",accpflag);
}

void _nanocarrier :: bond_formation(){
  extern _linklist link1;
  extern _receptor re;
  extern _membrane me;
  extern MTRand mtrand1;
  extern _intparam intparam;
  extern _datainput data;
  //extern _datawriter w;
  
  int i, j, sel_ves, sel_ab, sel_antigen,temp_bond[2];
  int N_sel, bt, *membcell,*lscl,*head, mem, current_multivalency;
  double Eot,Ebond_old,distance,theta,phi;
  double xpnt,ypnt,zpnt;
  double dist,sumw,wn,wo;
  double antImx,antImy,antImz,rrx,rry,rrz;
  int sel_ves_list[this->num_vesicles], sel_antigen_list[this->num_antigens];
  double max_nc_memb_distance, max_nc_anttip_distance, max_AB_ANT_distance;
  
  i = j = sel_ves = sel_ab = sel_antigen = -1;                     // Initialize various variables
  fill_n(sel_ves_list,this->num_vesicles,-1);                      // Initialize the array to -1
  
  lscl = link1.getlscl('m');                                       //antigenlink.lscl
  head = link1.gethead('m');                                       //antigenlink.head 
  
  int kv=-1; 
  
  double max_antlen = data.get_antigen_maxlength();
  double max_rxnlen = intparam.get_max_reaction_length();

  for(i=0; i<this->num_vesicles; i++){
	// do not check distance if already bound
	if (this->get_multivalency(i)>0){
        sel_ves_list[++kv] = i;
	}
	else {
		max_nc_memb_distance = pow(this->getsoft_radius(i) + max_antlen + max_rxnlen, 2);
		membcell = link1.get2ringcells(i,'v');
		for(j=0; j<125; j++){
      		mem = head[membcell[j]];
      		if (mem == -1) continue;
      		do{
       			rrx = this->getxc(i) - me.getxc(mem);
       			rry = this->getyc(i) - me.getyc(mem);
       			rrx = rrx - round(rrx/this->L)*this->L;           // vesicle membrane obey periodic boundary condition
       			rry = rry - round(rry/this->L)*this->L;
       			rrz = this->getzc(i) - me.getzc(mem);
       			distance = rrx*rrx + rry*rry + rrz*rrz;
       			if(distance <= max_nc_memb_distance){
        			sel_ves_list[++kv] = i;
        			break;                                        // breaks the do while loop since the vesicle is already selected
       			}
       			mem = lscl[mem];
      		} while (mem != -1);
     		if (kv!=-1 && sel_ves_list[kv]==i) break;            // breaks the for(j=...) loop since the vesicle i is already accounted for                                              
    	}
   	}
  }

  //select a NC randomly from the list
  if (kv==-1) return;
  if (kv==0)
    sel_ves = sel_ves_list[0];
  else
    sel_ves = sel_ves_list[mtrand1.randInt(kv)];

  double x_selves = this->getxc(sel_ves);
  double y_selves = this->getyc(sel_ves);
  double z_selves = this->getzc(sel_ves);
  current_multivalency = this->get_multivalency(sel_ves);
  max_nc_anttip_distance = pow(this->getsoft_radius(sel_ves)+2.0*intparam.get_max_reaction_length(),2);

//-------------------------------------

  membcell = link1.get2ringcells(sel_ves,'v');                              // 2 cell neighbourhood around the chosen vesicle
  lscl = link1.getlscl('c');                                                // antigenlink.lscl ( link list for the antigens)
  head = link1.gethead('c');                                                // antigenlink.head ( head list for antigen in a cell)
  fill_n(sel_antigen_list,this->num_antigens,-1);                           // Initialize an array sel_antigen_list, with size num_antigens to -1    
  int kant = 0;

  for(j=0; j<125; j++){
    mem = head[membcell[j]];
    if (mem == -1) continue;
     do{
      rrx = x_selves - re.getxt(mem);
      rry = y_selves - re.getyt(mem);
      rrz = z_selves - re.getzt(mem);
      rrx -= round(rrx/this->L)*this->L;
      rry -= round(rry/this->L)*this->L;
      distance = rrx*rrx + rry*rry + rrz*rrz;                               // compute distance between antigen tip and vesicle
      if(distance <= max_nc_anttip_distance){                               // choose all antigens within a cutoff
       sel_antigen_list[kant] = mem;                                        // Make a list of all unbonded antigens
       kant++;
      }
      mem = lscl[mem];
     } while(mem !=-1);
  }

  if (kant == 0) return;
  
  int num_ab = this->getnum_ab(sel_ves);                                            // Number of antibodies in the chosen vesicle
  int sel_antibody_list[num_ab*kant][2];                                          
  double x_selant, y_selant, z_selant;

  int kab = -1;    
  for (int k1=0; k1<kant; k1++){
    x_selant = re.getxt(sel_antigen_list[k1]);                                        // tip position for the selected antigen
    y_selant = re.getyt(sel_antigen_list[k1]);
    z_selant = re.getzt(sel_antigen_list[k1]);
    max_AB_ANT_distance = 2.0*pow(re.getlength(sel_antigen_list[k1]),2)*onemoneosqrt2;
    
  //------------------------------------
    for (j=0; j<num_ab; j++){ 
     rrx = x_selant - this->getx_antibody(sel_ves, j);
     rry = y_selant - this->gety_antibody(sel_ves, j);
     rrz = z_selant - this->getz_antibody(sel_ves, j);
     rrx -= this->L*round(rrx/this->L);
     rry -= this->L*round(rry/this->L);
     dist = (rrx*rrx + rry*rry + rrz*rrz);
     if (dist <= max_AB_ANT_distance){
      kab++;
      sel_antibody_list[kab][0]=j;
      sel_antibody_list[kab][1]=sel_antigen_list[k1];
     }
    }
  }

  if (kab == -1) return;
  if (kab > 0){
    int k11 = mtrand1.randInt(kab);
    sel_ab = sel_antibody_list[k11][0];                             // Randomly choose one pair of antigen and antibody
    sel_antigen = sel_antibody_list[k11][1];                 
  }
  else {
    sel_ab = sel_antibody_list[0][0];                               // Randomly choose one pair of antigen and antibody
    sel_antigen = sel_antibody_list[0][1];                 
  }
  
  if (this->getbonded(sel_ves,sel_ab) && (this->getbondedto(sel_ves,sel_ab) != sel_antigen)) return; // if bound, make sure it is bonded to the selected antigen
  if ((re.getbonded(sel_antigen) == 1) && (this->getbondedto(sel_ves,sel_ab) != sel_antigen) ) return;   // if bound, make sure it is bonded to the selected antigen

  x_selant = re.getxt(sel_antigen);                                  // tip position for the selected antigen
  y_selant = re.getyt(sel_antigen);
  z_selant = re.getzt(sel_antigen);

  double x_selab = this->getx_antibody(sel_ves,sel_ab);              // Tip of the antibody(sel_ab) in vesicle(sel_vesicle)
  double y_selab = this->gety_antibody(sel_ves,sel_ab);
  double z_selab = this->getz_antibody(sel_ves,sel_ab);

  double x_selantb = re.getxc(sel_antigen);                          // base of the antigen (sel_antigen)    
  double y_selantb = re.getyc(sel_antigen);                 
  double z_selantb = re.getzc(sel_antigen);      
  
//------------ Get the rotation matrix for the triangle or for the vertex depending on the location of the antigen ------------------
  tmatrix AMAT = this->construct_householder_matrix(sel_antigen);     // Get normal and rotation matrix for the selected antigen
  float dist_flex,dist_bond_spring;
  float flexed_x,flexed_y,flexed_z;
  
  //@ Parameters for the selected antigen
  double ant_size = re.getlength(sel_antigen);
  double stheta = re.get_theta(sel_antigen);
  int typant = re.gettype(sel_antigen);
  int typab = this->get_ab_type(sel_ves, sel_ab);
  double Dreaction = pow(intparam.get_delr(typant,typab),2);          // Square of the max reaction distance
//--------------------------------------------------------------------------------------------------------------------------------------------------

  sumw = 0.0;
  for (i=0 ; i<this->num_samp ; i++){                                 // generate num_samp antigen orientations
	  
  this->thetat[i] = acos(oneosqrt2 + onemoneosqrt2*mtrand1.rand());  // Pick up a polar angle between 0 and pi/4
  this->phit[i] = mtrand1.rand()*TWOPI;                              // random angle between 0 and 2*pi
  
  antImx =  ant_size*sin(this->thetat[i])*cos(this->phit[i]);          // x,y and z component assuming that the unflexed antigen is along the z direction
  antImy =  ant_size*sin(this->thetat[i])*sin(this->phit[i]);                       
  antImz =  ant_size*cos(this->thetat[i]);

  flexed_x = AMAT.hhm[0][0]*antImx + AMAT.hhm[0][1]*antImy + AMAT.hhm[0][2]*antImz;    
  flexed_y = AMAT.hhm[1][0]*antImx + AMAT.hhm[1][1]*antImy + AMAT.hhm[1][2]*antImz;   // transform the new orientation along z to the normal using Householder transf.
  flexed_z = AMAT.hhm[2][0]*antImx + AMAT.hhm[2][1]*antImy + AMAT.hhm[2][2]*antImz; 
  
  xpn[i] = x_selantb + flexed_x; 
  ypn[i] = y_selantb + flexed_y;
  zpn[i] = z_selantb + flexed_z;                                            // position of antigen tip in the Cartesian coordinates
  
  this->b[i] = mtrand1.randInt(1);              // random bond breakage (0=break,1=bond)
  dist_bond_spring = 0.0;

  if (this->b[i] !=0){
    rrx = xpn[i] - x_selab;                                // calculate the distance between antigen tip and antibody 
    rry = ypn[i] - y_selab;
    rrz = zpn[i] - z_selab;
    rrx -= this->L*round(rrx/this->L);                    // PBC along x and y   
    rry -= this->L*round(rry/this->L);
    dist_bond_spring = (rrx*rrx + rry*rry + rrz*rrz);     // computed with the frame of reference at (0,0,0)	
    if (dist_bond_spring < Dreaction){
	   dist_flex= ant_size*sin(this->thetat[i]);           // computed from the frame of reference at the vertex or triangle ( distance along the normal direction)
	   this->Eo[i] = this->calc_energy(sel_ves,sel_ab,sel_antigen,dist_bond_spring,dist_flex*dist_flex,this->b[i]);     // calculate the energy (flexure+ reaction)
	   this->weight[i] = exp(-this->beta*Eo[i]);                                                 // Boltzmann weight of the new configuration  
	}
	else {
      this->b[i] = 0;                                     // mark b=0 if the antigen and antibody are far apart
	}
  }
  
  if (this->b[i]==0){                                     //unbonded bond energy is zero
	  this->Eo[i]=0.0;
	  this->weight[i]=1.0;
  }
   sumw += this->weight[i];                                                                // Rosenbluth weight factor
  }

  wn = sumw;                                                                                // Rosenbluth factor for the new configurations 
  N_sel = this->select_from_Rosenbluth(weight,sumw);                                        // Select one configuration amongst the generated configs
  xpn[N_sel] -= round(xpn[N_sel]/this->L)*this->L;                                          // apply  PBC only for the chosen set
  ypn[N_sel] -= round(ypn[N_sel]/this->L)*this->L;

  //@ select num_samp-1 samples for the old configuration
  sumw = 0.0;
  if (re.getbonded(sel_antigen)==1){                                                        // If already bonded
   rrx = x_selant - x_selab;                                                                 //calculate the old bond energy
   rry = y_selant - y_selab;
   rrz = z_selant - z_selab;
   rrx -=  this->L*round(rrx/this->L);
   rry -=  this->L*round(rry/this->L);
   dist_bond_spring = (rrx*rrx + rry*rry + rrz*rrz);
   dist_flex= ant_size*(sin(stheta));
   Ebond_old = this->calc_energy(sel_ves,sel_ab, sel_antigen, dist_bond_spring,dist_flex*dist_flex,1);
   sumw += exp(-beta*Ebond_old);                                                             // Rosenbluth factor for the old configuration
  } 
  else {
	  sumw=1.0;
  }

  for (i=1; i<this->num_samp; i++) {                                                            // Generate another num_samp-1 configurations
    theta = acos(oneosqrt2+onemoneosqrt2*mtrand1.rand());                         // Pick up a polar angle between 0 and pi/4
    phi = mtrand1.rand()*TWOPI;

    antImx =  ant_size*sin(theta)*cos(phi);
    antImy =  ant_size*sin(theta)*sin(phi);                  // generate num_samp random configuration  by assuming that the antigen is along z
    antImz =  ant_size*cos(theta);

    flexed_x = AMAT.hhm[0][0]*antImx + AMAT.hhm[0][1]*antImy + AMAT.hhm[0][2]*antImz;    
    flexed_y = AMAT.hhm[1][0]*antImx + AMAT.hhm[1][1]*antImy + AMAT.hhm[1][2]*antImz;  // transform the new antigen orientation along z to the normal using Householder
    flexed_z = AMAT.hhm[2][0]*antImx + AMAT.hhm[2][1]*antImy + AMAT.hhm[2][2]*antImz;
	
    xpnt = x_selantb+flexed_x ;
	ypnt = y_selantb+flexed_y ;
	zpnt = z_selantb+flexed_z;
	
	if ( bt == 1 ){
		rrx = xpnt - x_selab;
		rrx -= round(rrx/this->L)*this->L;
		rry = ypnt - y_selab;                                           // calculate the distance between antigen tip and antibody
		rry -= round(rry/this->L)*this->L;
		rrz = zpnt - z_selab;
		dist_bond_spring = (rrx*rrx + rry*rry + rrz*rrz);
		if (dist_bond_spring < Dreaction){
			dist_flex= ant_size*(sin(theta));                        // computed from the frame of reference at the vertex or triangle
			Eot = this->calc_energy(sel_ves,sel_ab, sel_antigen, dist_bond_spring,dist_flex*dist_flex,bt);                             // calculate the energy
			sumw += exp(-this->beta*Eot);
		}
		else {
			bt=0;
		}
	}
	
	if (bt == 0) sumw += 1.0 ;
  }
  wo = sumw;
  
  if (wo == 0.0) wo = 1.0e-20;                                                            // Metropolis criteria
  
  if (mtrand1.rand() < wn/wo) {                                                           // accept the move, update the bond information
    if ((re.getbonded(sel_antigen)) && (this->getbondedto(sel_ves,sel_ab) == sel_antigen) && (b[N_sel] == 0)){                  // bond break
      this->setbonded(sel_ves,sel_ab,0);
      this->setbondedto(sel_ves,sel_ab,-1);
      re.setbonded(sel_antigen,0);
      re.reset_antigen_bond(sel_antigen);
      rrx = x_selantb + ant_size*AMAT.normal[0][0]; 
	  rrx -= round(rrx/this->L)*this->L;
      rry = y_selantb + ant_size*AMAT.normal[1][0];                                       
	  rry -= round(rry/this->L)*this->L;	  
      rrz = z_selantb + ant_size*AMAT.normal[2][0];    
      re.setantigenxt(sel_antigen,rrx);                                  // set antigen tip to point along the vertex/triangle normal
      re.setantigenyt(sel_antigen,rry);
      re.setantigenzt(sel_antigen,rrz);
      re.set_theta_phi(sel_antigen,0.0,0.0);                             // reset theta and phi values for the antigen to zero 
      this->update_bond_counters(sel_ves,-1,0,1);
    }

    if (b[N_sel] == 1) {
      if (current_multivalency == 0) {
       this->reset_euler_angles(sel_ves);                               // reset the Euler angles to (0,pi,0) when the first bond is formed
      }
      
      if (re.getbonded(sel_antigen) == 0) {
        this->update_bond_counters(sel_ves,1,1,0);
      }
      else {
        this->update_bond_counters(sel_ves,0,1,0);
      }
      
      this->setbonded(sel_ves,sel_ab,1);
      re.setbonded(sel_antigen,1);
      temp_bond[0] = sel_ves;
      temp_bond[1] = sel_ab;
      re.setbondedto(sel_antigen, temp_bond);
      this->setbondedto(sel_ves,sel_ab,sel_antigen);
      re.setantigenxt(sel_antigen,xpn[N_sel]);
      re.setantigenyt(sel_antigen,ypn[N_sel]);                      // update the antigen tip position ( This will also update the Fortran structure)
      re.setantigenzt(sel_antigen,zpn[N_sel]);  
      re.set_theta_phi(sel_antigen,thetat[N_sel],phit[N_sel]);      // set the current value of theta and phi for the antigens 
    }
  }
}

tmatrix _nanocarrier :: construct_householder_matrix(int sel_antigen){
  tmatrix antmatrix;
  int fantvert = (this->antptr+sel_antigen)->vertex;                // membrane vertex linked to the chose antigen (result is returned from fortran)
  int fanttri = (this->antptr+sel_antigen)->diffus_tri;             // membrane vertex linked to the chose antigen (result is returned from fortran)
  int antvert = fantvert-1;                                         // antvert=fantvert-1 is for use in C++ to retreive the values 
  int anttri = fanttri-1;
  
  if (anttri != -1){
   char marker='f';
   __module_curvcalc_MOD_compute_householdermatrix(&fanttri,&marker);   // compute the householder matrix for the face anttri ('f' is marker for triangle)
   for (int i=0; i<3; i++){
    antmatrix.normal[i][0]=(this->triptr+anttri)->fnor[i][0];            // use the face normal
    for (int j=0; j<3; j++){
     antmatrix.hhm[i][j]=(this->triptr+anttri)->HHM[i][j];
    }
   }
  }
  else
  {
    char marker='n';
    __module_curvcalc_MOD_compute_householdermatrix(&fantvert,&marker); // compute the householder matrix for the vertex antvert ('n' is marker for node)
    for (int i=0; i<3; i++){
     antmatrix.normal[i][0]=(this->verptr+antvert)->vnor[i][0];            // use the face normal
     for (int j=0; j<3; j++){
      antmatrix.hhm[i][j]=(this->verptr+antvert)->HHM[i][j];
     }
    }
  }
  return antmatrix;
}

double _nanocarrier :: calc_energy(int sel_ves, int sel_ab, int sel_ant, double bdistsq, double fdistsq, int bound){
  double energy = 0.0;
  extern _receptor re;
  extern _intparam intparam; 
  if (bound == 1) {
  	double kflex = re.get_flexure(sel_ant);
  	double* bbpar = intparam.get_bellbond_parameters(sel_ves,sel_ab,sel_ant);
  	double lambda = this->get_lambda(sel_ves);
  	energy = kflex*fdistsq*base_length_sq + lambda*(bbpar[0] + 0.5*bbpar[1]*bdistsq*base_length_sq);
  }
  return energy;
}

int _nanocarrier :: select_from_Rosenbluth(double *w, double sumw){
  extern MTRand mtrand1;
  double ws,cumw;
  int n;
  ws = mtrand1.rand(0.0,1.0)*sumw;
  cumw = w[0];
  n = 0;
  while (cumw < ws) {
   n++;
   cumw += w[n];
  }
  return n;
}

