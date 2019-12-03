#include "_declarations.h"
#include "_datainput.h"
#include "_fortran_structures.h"
#include "_parser.h"
#include "_nanocarrier.h"
#include "_interaction_parameters.h"
#include "_receptor.h"

void _datainput :: init(){
  extern bool debug_mode;
  this->read_antigen_parameters();
  this->read_system_geometry();
  this->periodic_box_length = double(this->N+1) * this->memb_init_bond_length;                         // nm 
  this->periodic_box_height = this->periodic_box_length * this->periodic_box_length_height_ratio;      // wrt to the periodic_box_length
  this->membrane_init_z = this->periodic_box_height * 0.25;                                   // Initial height at which the fluctuating membrane is placed.
  kBT = kb * this->temperature;                                                                                                      
   __module_datastruct_MOD_blen = this->periodic_box_length/(double)(this->N+1);             // scale box size with initial blength for setting up the membrane
  this->check_periodicbox_dimensions();
  this->initialize_fortran_variables();
  
  if (debug_mode){
	  cout<<str_red<<"\n \t <SYSTEM PARAMTERS>: READ FROM PARAMETERS FILE \t \t "<<endl<<str_black;
	  cout<<"\n \t  --> Periodic_box_length (nm) "<<this->periodic_box_length;
	  cout<<"\n\t  -->  Periodic_box_height (nm) "<<this->periodic_box_height;
	  cout<<"\n \t  --> Number of nanocarriers "<<this->num_nanocarriers;
	  cout<<"\n \t  --> Temperature in deg K "<<this->temperature;
	  cout<<"\n \t  --> Number of grid points for membrane in each direction "<<N;
	  cout<<"\n \t  --> Initial Membrane bond length (nm) "<<this->memb_init_bond_length;
	  cout<<"\n \t  --> OSD printing frequency "<<this->sampling_interval<<str_black<<endl;
	  cout<<endmarker<<endl;
  }
}

void _datainput :: read_antigen_parameters(){
  extern bool debug_mode;
  _PARSER parse;
  char DELIMITER=' ';
  string filename = "../PARAMETERS/antigen_parameters.inp";
  
  //@ total number of antigens
  auto pd1 = parse.PARSEFILE(filename,DELIMITER,"total");
  this->num_antigens = std::atoi(pd1[0][1].c_str());
  
  
   //@ how to pattern the antigens
  auto pd11 = parse.PARSEFILE(filename,DELIMITER,"pattern");
  this->pattern_ant = pd11[0][1].c_str();
  
  
  //@ total types of antigens
  auto pd2 = parse.PARSEFILE(filename,DELIMITER,"type");
  this->type_ant = std::atoi(pd2[0][1].c_str());
  
  
  //@Allocate antigen parameter arrays
  this->conc_ant = new(nothrow) double[this->type_ant];
  this->num_type_ant = new(nothrow) int[this->type_ant];
  this->radius_ant = new(nothrow) double[this->type_ant];
  this->length_ant = new(nothrow) double[this->type_ant];
  this->kflex = new(nothrow) double[this->type_ant];
  
  //@ parse the concentration
  int tempsum=0;
  int tint = -1;
  double tdoub = 0.0;
  
  auto pd3 = parse.PARSEFILE(filename,DELIMITER,"conc");
  for (int i=0; i<this->type_ant; i++){
    tint = std::atoi(pd3[i][0].c_str());
    tdoub = std::atof(pd3[i][1].c_str());
    this->conc_ant[tint] = tdoub;
    this->num_type_ant[tint] = round(tdoub*this->num_antigens);
    if (i<this->type_ant-1) tempsum += this->num_type_ant[tint];
  }
  this->num_type_ant[tint] = this->num_antigens - tempsum;  //@ to take any rounding problems into account
  
  
  //@ parse the radius
  auto pd4 = parse.PARSEFILE(filename,DELIMITER,"radius");
  for (int i=0; i<this->type_ant; i++){
    tint = std::atoi(pd4[i][0].c_str());
    tdoub = std::atof(pd4[i][1].c_str());
    this->radius_ant[tint] = tdoub;
  }
  
  //@ parse the receptor lengths
  auto pd5 = parse.PARSEFILE(filename,DELIMITER,"length");
  for (int i=0; i<this->type_ant; i++){
    tint = std::atoi(pd5[i][0].c_str());
    tdoub = std::atof(pd5[i][1].c_str());
    this->length_ant[tint] = tdoub;
  }
  
  this->min_antlen = this->max_antlen = this->length_ant[0];
  for (int i=1; i<this->type_ant; i++){
    if (this->length_ant[i] < this->min_antlen) this->min_antlen = this->length_ant[i];
    if (this->length_ant[i] > this->max_antlen) this->max_antlen = this->length_ant[i];
  }
    
  //@ parse the receptor flexure
  auto pd6 = parse.PARSEFILE(filename,DELIMITER,"flexure");
  for (int i=0; i<this->type_ant; i++){
    tint = std::atoi(pd6[i][0].c_str());
    tdoub = std::atof(pd6[i][1].c_str());
    this->kflex[tint] = (2*tdoub*1e-30)/pow(this->length_ant[tint]*base_length,3);           // supplied flex is in pN. nm^2 and length is in nm
    // See eqn.S2 in BPJ, 101, 2, 319-326, 2011
  }
  if (debug_mode) print_message_ant();
}

void _datainput :: print_message_ant(){
	cout<<str_red<<"\n\t <ANTIGEN>"<<endl<<str_black;
	cout<<"\t --> Total number of antigens = "<<this->num_antigens<<endl;
	cout<<"\t --> Antigens will be patterned in "<<this->pattern_ant<<endl;
	cout<<"\t --> Number of antigen types = "<<this->type_ant<<endl;
	cout<<"\t --> Concentration of antigens: ";
	for(int i=0; i<this->type_ant; i++){
		cout<<"  ANT"<<i<<" = "<<this->conc_ant[i];
		if (i<this->type_ant-1) cout<<" && ";
	}
	cout<<endl;
	cout<<"\t --> Number of antigens in each type: ";
	for(int i=0; i<this->type_ant; i++){
		cout<<"  ANT"<<i<<" = "<<this->num_type_ant[i];
		if (i<this->type_ant-1) cout<<" && ";
	}
	cout<<endl;    
	cout<<"\t --> Radius of antigens: ";
	for(int i=0; i<this->type_ant; i++){
		cout<<"  ANT"<<i<<" = "<<this->radius_ant[i]<<" nm ";
		if (i<this->type_ant-1) cout<<" && ";
	}
	cout<<endl;
	cout<<"\t --> Length of antigens: ";
	for(int i=0; i<this->type_ant; i++){
		cout<<"  ANT"<<i<<" = "<<this->length_ant[i]<<" nm ";
		if (i<this->type_ant-1) cout<<" && ";
	}
	cout<<endl;
	cout<<"\t --> Flexural rigidity of antigens: ";
	for(int i=0; i<this->type_ant; i++){
		cout<<"  ANT"<<i<<" = "<<this->kflex[i]<<" pN. nm^2 ";
		if (i<this->type_ant-1) cout<<" && ";
	}
	cout<<endl;
	cout<<endmarker<<endl;
}

void _datainput :: dump_antigen_parameters(ofstream &ofile){
	ofile<<"<ANTPARAM>"<<endl;
	ofile<<this->num_antigens<<" "<<this->pattern_ant<<" "<<this->type_ant<<endl;
	for (int i=0; i<this->type_ant; i++)
		ofile<<i<<" "<<this->conc_ant[i]<<" "<<this->num_type_ant[i]<<" "<<this->radius_ant[i]<<" "<<this->length_ant[i]<<" "<<this->kflex[i]<<endl;
	ofile<<"</ANTPARAM>"<<endl<<endl;
}

void _datainput :: read_antigen_parameters(ifstream &ifile){
	extern bool debug_mode;
	string temps;
	while (getline(ifile,temps)){
		if (temps.compare("<ANTPARAM>")==0){
			break;
		}
	}
	
	ifile>>this->num_antigens>>this->pattern_ant>>this->type_ant;
	
	//@Allocate antigen parameter arrays
	this->conc_ant = new(nothrow) double[this->type_ant];
	this->num_type_ant = new(nothrow) int[this->type_ant];
	this->radius_ant = new(nothrow) double[this->type_ant];
	this->length_ant = new(nothrow) double[this->type_ant];
	this->kflex = new(nothrow) double[this->type_ant];
	
	int tind;
	for (int i=0; i<this->type_ant; i++)
		ifile>>tind>>this->conc_ant[i]>>this->num_type_ant[i]>>this->radius_ant[i]>>this->length_ant[i]>>this->kflex[i];
	
	if (debug_mode) this->print_message_ant();
}

void _datainput :: read_system_geometry(){
  _PARSER parse;
  std::string filename = "../PARAMETERS/system-parameters.inp";
  char DELIMITER=',';
  
  //@ boxsize and membrane parameters
  this->periodic_box_length_height_ratio = parse.PARSEDOUBLE(filename, DELIMITER, "pbchoverl");
  this->num_nanocarriers = parse.PARSEINT(filename, DELIMITER, "num_nanocarriers");
  this->N = parse.PARSEINT(filename, DELIMITER, "memb_lattice_size");
  this->memb_init_bond_length = parse.PARSEDOUBLE(filename, DELIMITER, "memb_lattice_space");
  
  //@temperature
  this->temperature = parse.PARSEDOUBLE(filename, DELIMITER, "temperature");
  
  //@ runmode specifications
  this->run_system = parse.PARSESTRING(filename, DELIMITER, "run_system");
  this->mpi_mode = parse.PARSESTRING(filename, DELIMITER, "mpi_mode");
  
  this->win_start = parse.PARSEINT(filename, DELIMITER, "win_start");
  this->win_end = parse.PARSEINT(filename, DELIMITER, "win_end");
  
  //@system run parameters
  this->num_moves = parse.PARSEINT(filename, DELIMITER, "num_mc_moves");
  this->num_dtmc_moves = parse.PARSEINT(filename, DELIMITER, "num_dtmc_moves");
  this->sampling_interval = parse.PARSEINT(filename, DELIMITER, "sampling_interval");
  this->osd_print_interval = parse.PARSEINT(filename, DELIMITER, "osd_print_interval");
  this->window_antab_data_interval = parse.PARSEINT(filename, DELIMITER, "window_antab_data_interval");
  
  //@ Intervals to dump datafiles
  this->window_conf_interval = parse.PARSEINT(filename, DELIMITER, "window_conf_interval");
  this->antigen_memb_write_interval = parse.PARSEINT(filename, DELIMITER, "antigen_memb_write_interval");
  this->antigen_traj_interval = parse.PARSEINT(filename, DELIMITER, "antigen_traj_interval");
  
  //@ Binsizes
  this->binsize_rdf = parse.PARSEDOUBLE(filename, DELIMITER, "binsize_rdf");
  
  this->setup_runparameters();
  
}

void _datainput :: dump_system_geometry(ofstream &ofile){
	ofile<<"<SYSTEM>"<<endl;
	ofile<<"pbchoverl "<<this->periodic_box_length_height_ratio<<endl;
	ofile<<"boxlength "<<this->periodic_box_length<<endl;
	ofile<<"boxheight "<<this->periodic_box_height<<endl;
	ofile<<"num_nanocarriers "<<this->num_nanocarriers<<endl;
	ofile<<"memb_lattice_size "<<this->N<<endl;
	ofile<<"memb_lattice_space "<<this->memb_init_bond_length<<endl;
	ofile<<"temperature "<<this->temperature<<endl;
	ofile<<"run_system "<<this->run_system<<endl;
	ofile<<"mpi_mode "<<this->mpi_mode<<endl;
	ofile<<"win_start "<<this->win_start<<endl;
	ofile<<"win_end "<<this->win_end<<endl;
	ofile<<"num_mc_moves "<<this->num_moves<<endl;
	ofile<<"num_dtmc_moves "<<this->num_dtmc_moves<<endl;
	ofile<<"sampling_interval "<<this->sampling_interval<<endl;
	ofile<<"osd_print_interval "<<this->osd_print_interval<<endl;
	ofile<<"window_antab_data_interval "<<this->window_antab_data_interval<<endl;
	ofile<<"window_conf_interval "<<this->window_conf_interval<<endl;
	ofile<<"antigen_memb_write_interval "<<this->antigen_memb_write_interval<<endl;
	ofile<<"antigen_conf_interval "<<this->antigen_traj_interval<<endl;
	ofile<<"binsize_rdf "<<this->binsize_rdf<<endl;
	ofile<<"</SYSTEM>"<<endl<<endl;
}

void _datainput :: read_system_geometry(ifstream &ifile){
	extern bool debug_mode;
	string temps;
	while (getline(ifile,temps)){
		if (temps.compare("<SYSTEM>") == 0)
			break;
	}
	ifile>>temps>>this->periodic_box_length_height_ratio;
	ifile>>temps>>this->periodic_box_length;
	ifile>>temps>>this->periodic_box_height;
	ifile>>temps>>this->num_nanocarriers;
	ifile>>temps>>this->N;
	ifile>>temps>>this->memb_init_bond_length;
	ifile>>temps>>this->temperature;
	ifile>>temps>>this->run_system;
	ifile>>temps>>this->mpi_mode;
	ifile>>temps>>this->win_start;
	ifile>>temps>>this->win_end;
	ifile>>temps>>this->num_moves;
	ifile>>temps>>this->num_dtmc_moves;
	ifile>>temps>>this->sampling_interval;
	ifile>>temps>>this->osd_print_interval;
	ifile>>temps>>this->window_antab_data_interval;
	ifile>>temps>>this->window_conf_interval;
	ifile>>temps>>this->antigen_memb_write_interval;
	ifile>>temps>>this->antigen_traj_interval;
	ifile>>temps>>this->binsize_rdf;
	
	kBT = kb * this->temperature;                                                                                                      
	if (debug_mode){
		cout<<str_red<<"\n \t <SYSTEM PARAMTERS>: READ FROM RESTART FILE \t \t "<<endl<<str_black;
		cout<<"\n \t  --> Periodic_box_length (nm) "<<this->periodic_box_length;
		cout<<"\n\t  -->  Periodic_box_height (nm) "<<this->periodic_box_height;
		cout<<"\n \t  --> Number of nanocarriers "<<this->num_nanocarriers;
		cout<<"\n \t  --> Temperature in deg K "<<this->temperature;
		cout<<"\n \t  --> Number of grid points for membrane in each direction "<<N;
		cout<<"\n \t  --> Initial Membrane bond length (nm) "<<this->memb_init_bond_length;
		cout<<"\n \t  --> OSD printing frequency "<<this->sampling_interval<<str_black<<endl;
		cout<<endmarker<<endl;
	}
}

double _datainput :: get_max_ncsize(){
	stringstream filename;
	_PARSER parse;
	char DELIMITER = ',';
	std::array<double,MAX_FILE_ENTRY> pars_double;
	double maxrad=0.0;
	
	for (int i=0; i<this->num_nanocarriers; i++){
		filename<<"../PARAMETERS/nc-"<<i<<".ncinp";
		auto parsedata = parse.PARSEFILE(filename.str().c_str(),DELIMITER);
		pars_double = parse.PARSEDOUBLE(parsedata,"nc_radius",1);
		if (i==0){
			maxrad = pars_double[0];
		}
		else {
			if (maxrad < pars_double[0])
				maxrad = pars_double[0];
		}
	}
	return maxrad;
}

//@
void _datainput :: check_periodicbox_dimensions(){
  extern bool debug_mode;
  extern _receptor re;
  extern _intparam intparam;
  
   int boxflag=0;
  //@Max NC radius
   double sr1 = this->get_max_ncsize();
   
   if (debug_mode) cout<<str_red<<"\n \t <BOXINFO> "<<str_black<<endl;
   
  if (this->membrane_init_z + 3*sr1 + 25.0 > this->periodic_box_height){       // Lengths required in init.cpp for placing vesicles in the box
   boxflag=1;
   if (debug_mode)
   		cout<<"\t\t\t\t  Current box height "<<periodic_box_height<<" not in the prescribed limit"<<endl;
  }
  
  //@ reset box size
  if (boxflag==1){
   do{
    this->periodic_box_height = 1.05 * this->periodic_box_height;
     this->membrane_init_z = this->periodic_box_height * 0.25;
     cout<<"\t\t\t\t\t --> Retrying with box height "<< this->periodic_box_height<<endl;
     if (this->membrane_init_z +  3 * sr1 + 25.0 < this->periodic_box_height)                       // Lengths required in init.cpp for placing vesicles in the box
       boxflag=0;
    } while(boxflag==1);
  
   this->periodic_box_length_height_ratio = this->periodic_box_height/this->periodic_box_length;
   
   if (debug_mode)
      cout <<"\n \t\t\t\t Rescaling average bond distance to "<<__module_datastruct_MOD_blen<<endl;
  }
  
  if (debug_mode){
    cout<<"\t\t\t\t Final Box length and height are "<<periodic_box_length<<" "<<periodic_box_height<<endl;
	cout<<endmarker<<endl;
  }
}

void _datainput :: resetperiodic_box_dimensions(){
	extern _nanocarrier nc;
	nc.reset_box_dimensions();
}

//@
void _datainput :: setup_runparameters(){
  extern int processor, nprocessors;

  //@ details of the processor number  
  this->rparam.curr_proc = processor;
  this->rparam.total_proc = nprocessors;
  
  
  this->rparam.system = this->run_system;
  
  //@ how many types of MC moves are allowed (36 if all moves are present)
  if (this->run_system.compare("nc_membrane")==0 ){
    this->rparam.nmcss = 0.0;
    this->rparam.nmcse = 36.0;
    this->rparam.DTMCsteps = this->num_dtmc_moves;
    this->rparam.DTMCflag = true;
    this->rparam.sampling_cutoff = this->num_moves/2;
  }
  else if (this->run_system.compare("nc_planar")==0 ){
    this->rparam.nmcss = 0.0;
    this->rparam.nmcse = 24.0;  
    this->rparam.DTMCsteps = 0;
    this->rparam.DTMCflag = false;
    this->rparam.sampling_cutoff = this->num_moves/2;
  }
  else if (this->run_system.compare("membrane_DTMC")==0 ) {
    this->rparam.nmcss = 24.0;
    this->rparam.nmcse = 36.0;  
    this->rparam.DTMCsteps = this->num_dtmc_moves;    
    this->rparam.DTMCflag = true;
    this->rparam.sampling_cutoff = this->num_moves/2;
  }
  
  this->rparam.mpi_mode = this->mpi_mode;
  //@ winstart and win_end
  
  if (this->mpi_mode.compare("serial")==0){
    //@ a single processor loops over all the PMF windows
    this->rparam.win_start = this->win_start;
    this->rparam.win_end = this->win_end;
	this->rparam.ensno=processor;
  }
  
  else if (this->mpi_mode.compare("parallel")==0){
    //@ each window in a PMF calculation is given to one processor
    //@ the reference value of a window will be supplied by rref and href in nanocarrier.cpp
    this->rparam.win_start = this->win_start + processor;
    this->rparam.win_end = this->win_start + processor + 1;
	this->rparam.ensno=0;
  }
  
  this->rparam.step_start = 0;
  this->rparam.step_end = this->num_moves;
  this->rparam.samp_freq = this->sampling_interval;
  this->rparam.osdp_freq = this->osd_print_interval;
  this->rparam.antab_freq = this->window_antab_data_interval;
  this->rparam.conf_freq = this->window_conf_interval;
  this->rparam.ant_write_freq = this->antigen_memb_write_interval;
  this->rparam.ant_traj_freq = this->antigen_traj_interval;
}

void  _datainput :: dump_run_parameters(ofstream &ofile){
	ofile<<"<RUNPARAMETERS>"<<endl;
	ofile<<"proc "<<this->rparam.curr_proc<<" "<<this->rparam.total_proc<<endl;
	ofile<<"system "<<this->rparam.system<<endl;
	ofile<<"mcs "<<this->rparam.nmcss<<" "<<this->rparam.nmcse<<" "<<this->rparam.DTMCflag<<" "<<this->rparam.DTMCsteps<<" "<<this->rparam.sampling_cutoff<<endl;
	ofile<<"mpi_mode "<<this->rparam.mpi_mode<<endl;
	ofile<<"wininfo "<<this->rparam.win_start<<" "<<this->rparam.win_end<<" "<<this->rparam.ensno<<endl;
	ofile<<"step "<<this->rparam.step_start<<" "<<this->rparam.step_end<<" "<<this->rparam.samp_freq<<" "<<this->rparam.osdp_freq<<" "<<this->rparam.antab_freq<<" "<<this->rparam.conf_freq<<" "<<this->rparam.ant_write_freq<<" "<<this->rparam.ant_traj_freq<<endl;
	ofile<<"</RUNPARAMETERS>"<<endl<<endl;
}

void  _datainput :: read_run_parameters(int ensno, int window, ifstream &ifile){
	extern bool debug_mode;
	string temps;
	while (getline(ifile,temps)){
		if (temps.compare("<RUNPARAMETERS>") == 0)
			break;
	}
	ifile>>temps>>this->rparam.curr_proc>>this->rparam.total_proc;
	ifile>>temps>>this->rparam.system;
	ifile>>temps>>this->rparam.nmcss>>this->rparam.nmcse>>this->rparam.DTMCflag>>this->rparam.DTMCsteps>>this->rparam.sampling_cutoff;
	ifile>>temps>>this->rparam.mpi_mode;
	ifile>>temps>>this->rparam.win_start>>this->rparam.win_end>>this->rparam.ensno;
	ifile>>temps>>this->rparam.step_start>>this->rparam.step_end>>this->rparam.samp_freq>>this->rparam.osdp_freq>>this->rparam.antab_freq>>this->rparam.conf_freq>>this->rparam.ant_write_freq>>rparam.ant_traj_freq;
	
	this->rparam.ensno = ensno;
	this->rparam.win_start = window+1;  //@ Set up for running the next window
	if (debug_mode) this->print_runparam_message();
}

void _datainput :: print_runparam_message(){
	cout<<str_red<<"\n\t\t<RUNPARAMETERS>\t\t"<<endl<<str_black;
	cout<<"\t --> proc: "<<this->rparam.curr_proc<<" "<<this->rparam.total_proc<<endl;
	cout<<"\t --> system: "<<this->rparam.system<<endl;
	cout<<"\t --> mcs: "<<this->rparam.nmcss<<" "<<this->rparam.nmcse<<" "<<this->rparam.DTMCflag<<" "<<this->rparam.DTMCsteps<<" "<<this->rparam.sampling_cutoff<<endl;
	cout<<"\t --> mpi_mode: "<<this->rparam.mpi_mode<<endl;
	cout<<"\t --> wininfo: "<<this->rparam.win_start<<" "<<this->rparam.win_end<<" "<<this->rparam.ensno<<endl;
	cout<<"\t --> step: "<<this->rparam.step_start<<" "<<this->rparam.step_end<<" "<<this->rparam.samp_freq<<" "<<this->rparam.osdp_freq<<" "<<this->rparam.antab_freq<<" "<<this->rparam.conf_freq<<" "<<this->rparam.ant_write_freq<<" "<<this->rparam.ant_traj_freq<<endl<<endl;
}

runparam _datainput :: get_runparameters(){
  return this->rparam;
}

//@
int _datainput :: get_types_of_antigens(){
  return this->type_ant;
}

int* _datainput :: get_number_of_antigen_types(){
  return this->num_type_ant;
}

double* _datainput :: get_antigen_radius(){
  return this->radius_ant;
}

double _datainput :: get_antigen_radius(int antnum){
  return this->radius_ant[antnum];
}

double* _datainput :: get_antigen_length(){
  return this->length_ant;
}

double* _datainput :: get_antigen_flexure(){
  return this->kflex;
}

string _datainput :: get_antigen_pattern(){
  return this->pattern_ant;
}

double _datainput :: get_antigen_maxlength(){
  return this->max_antlen;
}

double _datainput :: get_antigen_minlength(){
  return this->min_antlen;
}

//@
double _datainput ::  getperiodic_box_length(){
  return periodic_box_length;
}

void _datainput ::  resetperiodic_box_dimensions(double boxlength){                    // reset box dimensions if the fortran code adjusts link lengths
  periodic_box_length=boxlength;
  periodic_box_height = periodic_box_length*periodic_box_length_height_ratio;                                           
}

double _datainput ::  getperiodic_box_height(){
  return periodic_box_height;
}

int _datainput ::  getnum_members(char ch){
 if (ch == 'v')
   return num_nanocarriers;
 if (ch =='c')
   return num_antigens;
 if (ch == 'm')
   return __module_datastruct_MOD_nver;
 return 0.0;
}

void _datainput ::  setnum_members(char ch, int nval){
	if (ch == 'v')
		num_nanocarriers = nval;
	if (ch =='c')
		num_antigens = nval;
	if (ch == 'm')
		__module_datastruct_MOD_nver = nval;
}

int _datainput ::  getnum_moves(){
  return num_moves;
}

int _datainput ::  getnum_eq_moves(){
  return num_eq_moves;
}

double _datainput :: gettemperature(){
  return temperature;
}

int _datainput ::  getN(){
  return N;
}

void _datainput ::  setN(int N_new){
  N=N_new;
}

int _datainput :: getsampling_interval(){
  return this->sampling_interval;
}

double _datainput :: getmembrane_init_z(){
  return this->membrane_init_z;
}

double _datainput :: get_kBT(){
	return this->kBT;
}

double _datainput :: get_binsize_rdf(){
	return this->binsize_rdf;
}


int _datainput :: get_window_conf_interval(){
	return this->window_conf_interval;
}

int _datainput :: get_window_antab_data_interval(){
	return this->window_antab_data_interval;
}

int _datainput :: get_antigen_memb_write_interval(){
	return this->antigen_memb_write_interval;
}

int _datainput :: get_antigen_traj_interval(){
	return this->antigen_traj_interval;
}

int _datainput :: get_antigen_traj_samples(){
	return (this->num_moves/2/this->antigen_traj_interval);
}

int _datainput :: get_osd_print_interval(){
  return this->osd_print_interval;
}

void _datainput :: initialize_fortran_variables(){
  __module_datastruct_MOD_blen = memb_init_bond_length;                        // scale box size with initial blength for setting up the membrane
  __module_datastruct_MOD_gsize = getN()+1;                                    // Pass the grid size to the fortran subroutine (0 -> N) for fortran
  __module_datastruct_MOD_num_antigens = getnum_members('c');                  // Pass the number of antigen to the fortran
  __module_datastruct_MOD_num_nanocarrier = getnum_members('v');               // Pass the number of nanocarriers to the fortran
  __module_datastruct_MOD_periodic_box_length =getperiodic_box_length();       // Pass the box length
  __module_datastruct_MOD_periodic_box_height =getperiodic_box_height();       // Pass the box height
  __module_datastruct_MOD_membrane_init_z = getmembrane_init_z();              // Initial height at which the membrane is initialized
  __module_datastruct_MOD_beta = (double) 300.0/gettemperature();
}
