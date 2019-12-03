#include "_declarations.h"
#include "_restart.h"
#include "_datainput.h"        // definition for data reading from input.in
#include "_setup.h"            // Lays down the datastructure for system variables
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_membrane.h"
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_interaction_parameters.h"

void _restart :: write_restart_conf(int ens, int window, int step){
	ofstream ofile;
	stringstream st1;
	
	st1<<"../SYS_STATE/DUMP/restart-ENS-"<<ens<<"-window-"<<window<<"-frame-"<<step<<".dump"; 
	ofile.open(st1.str(),ios::binary|ios::out);
	ofile<<"ENSEMBLE, WINDOW, STEP "<<ens<<" "<<window<<" "<<step<<endl;
	this->_write_restart_conf(ofile);
	ofile.close();
}

void _restart :: write_restart_conf(int ens, int window){
	ofstream ofile;
	stringstream st1;
	
	st1<<"../SYS_STATE/DUMP/restart-ENS-"<<ens<<"-window-"<<window<<".dump"; 
	ofile.open(st1.str(),ios::binary|ios::out);
	ofile<<"ENSEMBLE, WINDOW, STEP "<<ens<<" "<<window<<" "<<0<<endl;
	this->_write_restart_conf(ofile);
	ofile.close();
}

void _restart :: write_restart_conf(int window){
	ofstream ofile;
	stringstream st1;
	
	st1<<"../SYS_STATE/DUMP/restart-window-"<<window<<".dump"; 
	ofile.open(st1.str(),ios::binary|ios::out);
	ofile<<"ENSEMBLE, WINDOW, STEP "<<0<<" "<<window<<" "<<0<<endl;
	this->_write_restart_conf(ofile);
	ofile.close();
}

void _restart :: _write_restart_conf(ofstream &ofile){
  extern _nanocarrier nc;
  extern _receptor re;
  extern _intparam intparam;
  extern _datainput data;
  
  this->num_ves = nc.getnum_members();
  this->num_ant = re.getnum_members();
  
  ofile<<"TOTALNC "<<this->num_ves<<endl;
  ofile<<"TOTAL_ANTIGEN "<<this->num_ant<<endl<<endl;
  
  data.dump_system_geometry(ofile);
  data.dump_run_parameters(ofile);
  data.dump_antigen_parameters(ofile);
  intparam.dump_parameters(ofile);
    
  ofile<<"<NCDATA>"<<endl;
  for (i=0; i< this->num_ves; i++){                                // write the vesicle position
	  nc.dump_conf(i,ofile);
  }
  ofile<<"</NCDATA>"<<endl<<endl;
  re.dump_conf(ofile);
}

void _restart :: read_restart_conf(int procnum){
	extern _nanocarrier nc;
	extern _receptor re;
	extern _datainput data;
	extern _intparam intparam;
	extern _membrane me;
	
	int tnum_nc,tnum_antigen,ens,window,step;
	ifstream ifile;
	string temps;	
	stringstream st1;
	
	st1<<"../SYS_STATE/RESTART/restart-conf-"<<procnum<<".dump"; 
	ifile.open(st1.str(),ios::binary|ios::in);	
	
	ifile>>temps>>temps>>temps>>ens>>window>>step;
	ifile>>temps>>tnum_nc;	
	ifile>>temps>>tnum_antigen;
	
	if ((tnum_nc != data.getnum_members('v')) || (tnum_antigen != data.getnum_members('c'))){
	}
	
	data.read_system_geometry(ifile);
	data.read_run_parameters(ens,window,ifile);
	data.read_antigen_parameters(ifile);
	intparam.read_parameters(ifile);
	nc.init(ifile);
	re.init(ifile);
	ifile.close();
	me.init();
	
	for(int i=0; i<num_ves; i++){
		nc.compute_shadow(i);
	}
	
}