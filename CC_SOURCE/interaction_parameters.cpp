#include "_declarations.h"
#include "_interaction_parameters.h"
#include "_datainput.h"
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_parser.h"
#include "_constant.h"

void _intparam :: init(){
  extern _nanocarrier nc;
  extern _datainput data;  
  this->nabtype = nc.get_num_abtype(0);
  this->nanttype = data.get_types_of_antigens();
  this->L = data.getperiodic_box_length();
  cout<<"Coming here 1.1"<<endl;
  this->set_array_sizes();
  this->read_interaction_parameters();
  this->set_minmax_reaction_length();
}

void _intparam :: dump_parameters(ofstream &ofile){
	ofile<<"<INTPARAM>"<<endl;
	ofile<<this->nabtype<<" "<<this->nanttype<<" "<<this->L<<endl;
	ofile<<this->min_delr<<" "<<this->max_delr<<endl;
	for (int i=0; i<this->nabtype; i++){
		for (int j=0; j<this->nanttype; j++){
			ofile<<i<<" "<<j<<" "<<this->delg[i][j]<<" "<<this->kbond[i][j]<<" "<<this->delr[i][j]<<endl;
		}
	}
	ofile<<"</INTPARAM>"<<endl<<endl;
}

void _intparam :: read_parameters(ifstream &ifile){
	string temps;
	while (getline(ifile,temps)){
		if (temps.compare("<INTPARAM>") == 0)
			break;
	}

	ifile>>this->nabtype>>this->nanttype>>this->L;
	this->set_array_sizes();
	ifile>>this->min_delr>>this->max_delr;
	int a,b;
	double c,d,e;
  cout<<"Coming here 1.1"<<endl;
	for (int i=0; i<this->nabtype; i++){
		for (int j=0; j<this->nanttype; j++){
			ifile>>a>>b>>c>>d>>e;
      cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl;
			this->delg[a][b] = c;
			this->kbond[a][b] = d;
			this->delr[a][b] = e;
			this->delrsq[a][b] = e*e;
		}
	}
}

void _intparam :: set_array_sizes(){
	this->delg = new(nothrow) double*[this->nanttype];
	this->kbond = new(nothrow) double*[this->nanttype];
	this->delr = new(nothrow) double*[this->nanttype];
	this->delrsq = new(nothrow) double*[this->nanttype];
	
	for(int i=0; i<this->nanttype; i++){
		this->delg[i] = new(nothrow) double[this->nabtype];
		this->kbond[i] = new(nothrow) double[this->nabtype];
		this->delr[i] = new(nothrow) double[this->nabtype];
		this->delrsq[i] = new(nothrow) double[this->nabtype];
		fill_n(this->delg[i],this->nabtype,0.0);
		fill_n(this->kbond[i],this->nabtype,0.0);
		fill_n(this->delr[i],this->nabtype,0.0);
	}
}


void _intparam :: reset_box_dimensions(){
  extern _datainput data;
  this->L = data.getperiodic_box_length();
}

void _intparam :: read_interaction_parameters(){
  _PARSER parse;
  char DELIMITER=',';
  extern bool debug_mode;
  //@ parse for the spring constant
  

  
  auto parsekbond = parse.PARSEFILE("../PARAMETERS/int_param.inp",DELIMITER,"kbond");
  for (std::size_t i=0; i<(std::size_t) (this->nanttype*this->nabtype); i++){
      int ind1 = std::atoi(parsekbond[i][0].c_str());
      int ind2 = std::atoi(parsekbond[i][1].c_str());
      double kbo = std::atof(parsekbond[i][2].c_str());
      kbond[ind1][ind2] = kbo;
  } 
  
  //@ parse for delg
  auto parsedelg = parse.PARSEFILE("../PARAMETERS/int_param.inp",DELIMITER,"delg");  
  for (std::size_t i=0; i<(std::size_t) (this->nanttype*this->nabtype); i++){
	  int ind1 = std::atoi(parsedelg[i][0].c_str());
	  int ind2 = std::atoi(parsedelg[i][1].c_str());
	  double kbo = std::atof(parsedelg[i][2].c_str());
	  cout<<ind1<<" "<<ind2<<" "<<kbo<<endl;
	  delg[ind1][ind2] = kbo;
  } 
  
 
  
  //@ compute delr from delg and kbond
  for (int i=0; i<this->nanttype; i++){
    for (int j=0; j<this->nabtype; j++){
      if (abs(delg[i][j])>0.0){
		delr[i][j] = sqrt(2.0*abs(delg[i][j])/kbond[i][j])/base_length;                    // delr is given in the base length of the code
      }
      else{
      	delr[i][j] = 0.0;
      }
      delrsq[i][j] = pow(delr[i][j],2);          // square of delr to minimize the number of calculations 
	  cout<<i<<" "<<j<<" "<<delrsq[i][j]<<endl;
	  
	  
    }
  }  
  
  if (debug_mode){
  cout<<str_red<<"\n\t <INTERACTION PARAMETERS> \n"<<str_black<<endl;
  cout<<str_blue<<"\t \t --> Spring constants (k) \n"<<str_black<<endl;
  for (int i=0; i<this->nanttype; i++){
    for (int j=0; j<this->nabtype; j++){
      cout<<"\t \t \t ANT"<<i<<" "<<"\t AB"<<j<<" \t"<<kbond[i][j]<<" N/m"<<endl;
    }
  }
  cout<<endl;
  
  cout<<str_blue<<"\t\t --> Activation energy for AB-ANT pairs (delg) \n"<<str_black<<endl;    
  for (int i=0; i<this->nanttype; i++){
    for (int j=0; j<this->nabtype; j++){
      cout<<"\t \t \t ANT"<<i<<" "<<"\t AB"<<j<<" \t"<<delg[i][j]<<" J"<<endl;
    }
  }
  cout<<endl;
  
  cout<<str_blue<<"\t \t --> Cutoff distance for receptor-ligand bonds \n"<<str_black<<endl;
  for (int i=0; i<this->nanttype; i++){
    for (int j=0; j<this->nabtype; j++){
      cout<<"\t \t \t ANT"<<i<<" "<<"\t AB"<<j<<" \t"<<delr[i][j]<<" nm"<<endl;
    }
  }
  cout<<endl;
  cout<<endmarker<<endl;
  }
}

void _intparam :: set_minmax_reaction_length(){
    this->min_delr = this->max_delr = delr[0][0];
    for (int i=0; i<this->nanttype; i++){
     for (int j=0; j<this->nabtype; j++){
      if (delr[i][j] < this->min_delr) this->min_delr = delr[i][j];
      if (delr[i][j] > this->max_delr) this->max_delr = delr[i][j];     
     }
    } 
}

double _intparam :: get_kbond(int i, int j){
  return this->kbond[i][j];
}

double _intparam :: get_kbond(int ncnum, int abnum, int antnum){
  extern _nanocarrier nc;
  extern _receptor re;
  int i = re.gettype(antnum);
  int j = nc.get_ab_type(ncnum, abnum);
  return this->kbond[i][j];
}

double _intparam :: get_delg(int i, int j){
  return this->delg[i][j];
}

double _intparam :: get_delg(int ncnum, int abnum, int antnum){
  extern _nanocarrier nc;
  extern _receptor re;
  int i = re.gettype(antnum);
  int j = nc.get_ab_type(ncnum, abnum);
  return this->delg[i][j];
}
  
double _intparam :: get_delr(int i, int j){
  return this->delr[i][j];
}

double _intparam :: get_delrsq(int i, int j){
  return this->delrsq[i][j];
}

double _intparam :: get_delr(int ncnum, int abnum, int antnum){
  extern _nanocarrier nc;
  extern _receptor re;
  int i = re.gettype(antnum);
  int j = nc.get_ab_type(ncnum, abnum);
  return this->delr[i][j];
}

double _intparam :: get_delrsq(int ncnum, int abnum, int antnum){
  extern _nanocarrier nc;
  extern _receptor re;
  int i = re.gettype(antnum);
  int j = nc.get_ab_type(ncnum, abnum);
  return this->delrsq[i][j];
}

double* _intparam :: get_bellbond_parameters(int ncnum, int abnum, int antnum){
  extern _nanocarrier nc;
  extern _receptor re;
  int t1 = re.gettype(antnum);
  int t2 = nc.get_ab_type(ncnum, abnum);
  this->param[0] = this->delg[t1][t2];
  this->param[1] = this->kbond[t1][t2];
  this->param[2] = this->delr[t1][t2];
  return this->param;
}

double _intparam :: get_min_reaction_length(){
  return this->min_delr;
}

double _intparam :: get_max_reaction_length(){
  return this->max_delr;
}
