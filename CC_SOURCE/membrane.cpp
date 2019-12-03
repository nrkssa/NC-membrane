#include "_declarations.h"
#include "_constant.h"
#include "_membrane.h"
#include "_datainput.h"                            // definition for data reading from input.in
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_linklist.h"
#include "_selfavoidance.h"
#include "_debugger.h"

void _MEMBRANE ::init() {
  __module_makesurface_MOD_makesurface();                   // call to build the membrane
  this->nver  =  __module_datastruct_MOD_nver;              // store the parameters
  this->pbnum =  __module_datastruct_MOD_pbnum;
  this->ntr  =  __module_datastruct_MOD_ntr;
  this->nlink  =  __module_datastruct_MOD_tlink;
  this->verptr = &__module_datastruct_MOD_ver;
}

_MEMBRANE ::_MEMBRANE(int nver1, int pbnum1, int ntr1, int nlink1){
	this->nver = nver1;
	this->pbnum = pbnum1;
	this->ntr = ntr1;
	this->nlink = nlink1;
	this->verptr = &__module_datastruct_MOD_ver;
}

int _MEMBRANE :: getnum_vertices(){
   return this->nver;
}

int _MEMBRANE :: getnum_triangles(){
	return this->ntr;
}

int _MEMBRANE :: getnum_links(){
	return this->nlink;
}

int _MEMBRANE :: getnum_verticestot(){
	return this->pbnum;
}

double _MEMBRANE :: getxc(int i) {
  return (this->verptr+i)->vcoord[0][0];
 }

double _MEMBRANE :: getyc(int i) {
  return (this->verptr+i)->vcoord[1][0];
 }

double _MEMBRANE :: getzc(int i) {
  return (this->verptr+i)->vcoord[2][0];
 }

void _MEMBRANE :: setxc(int i,double rx) {
  (this->verptr+i)->vcoord[0][0] = rx;
 }

void _MEMBRANE :: setyc(int i,double ry) {
  (this->verptr+i)->vcoord[1][0] = ry; }

void _MEMBRANE :: setzc(int i,double rz) {
  (this->verptr+i)->vcoord[2][0] = rz;
 }

void _MEMBRANE :: setcellno(int i, int cellno){
  (this->verptr+i)->cellno = cellno;
}

void _membrane :: init() {
  extern _datainput data;
  flmemb.init();
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();
}

_membrane :: _membrane(double L1, double H1, int nver1, int pbnum1, int ntr1, int nlink1){
	this->L = L1;
	this->H = H1;
	_MEMBRANE(nver1, pbnum1, ntr1, nlink1);
}

void _membrane :: reset_box_dimensions(){
  extern _datainput data;
  this->L = data.getperiodic_box_length();
  this->H = data.getperiodic_box_height();
}

int _membrane :: getnum_vertices(){
  return flmemb.getnum_vertices();
}

int _membrane :: getnum_triangles(){
	return flmemb.getnum_triangles();
}

int _membrane :: getnum_links(){
	return flmemb.getnum_links();
}

int _membrane :: getnum_verticestot(){
	return flmemb.getnum_verticestot();
}

double _membrane :: getL(){
	return this->L;
}

double _membrane :: getH(){
	return this->H;
}

void _membrane :: setxc(int i, double x){
	flmemb.setxc(i,x);
}

double _membrane :: getxc(int i){
	return flmemb.getxc(i);
}

void _membrane :: setyc(int i, double y){
	flmemb.setyc(i,y);
}

double _membrane :: getyc(int i){
	return flmemb.getyc(i);
}

void _membrane :: setzc(int i, double z){
	flmemb.setzc(i,z);
}

double _membrane :: getzc(int i){
	return flmemb.getzc(i);
}

void _membrane :: setcellno(int i,int cellno){
	flmemb.setcellno(i,cellno);
}

void _membrane :: movevertex(bool eqmflag){
  if (eqmflag)
    __module_mcsmoves_MOD_link_flip();
  else
  __module_mcsmoves_MOD_vertex_move_biased();
}

void _membrane :: fliplink(bool eqmflag){
  if (eqmflag)
    __module_mcsmoves_MOD_vertex_move();
  else 
   __module_mcsmoves_MOD_link_flip_biased();
}

void _membrane :: membrane_montecarlo(bool eqmflag){
  if (eqmflag)
    __module_mcsmoves_MOD_membrane_montecarlo();
  else
    __module_mcsmoves_MOD_membrane_montecarlo_biased();
}

void _membrane :: dump_membrane_conf(int ensno){
	__module_writedata_MOD_membrane_dump1(&ensno);                                    
}

void _membrane :: dump_membrane_conf(int ensno, int window){
	__module_writedata_MOD_membrane_dump2(&ensno,&window);                                    
}


