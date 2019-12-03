#include "_declarations.h"
#include "_setup.h"
#include "_nanocarrier.h"
#include "_membrane.h"
#include "_receptor.h"

void _setup :: init() {
}

_setup :: _setup(){
}

int _setup :: getnum_members(char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;

	if (ch =='v')
		return nc.getnum_members();
	if (ch == 'c')
		return re.getnum_members();
	if (ch=='m')
		return me.getnum_vertices();
	return 0;
}

void _setup :: setxc(int i, double x, char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;
	if (ch == 'v')
		nc.setxc(i,x);
	if (ch == 'c')
		re.setxc(i,x);
	if (ch =='m')
		me.setxc(i,x);
}

void _setup :: setyc(int i, double x, char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;

	if (ch == 'v')
		nc.setyc(i,x);
	if (ch == 'c')
		re.setyc(i,x);
	if (ch =='m')
		me.setyc(i,x);
}


void _setup :: setzc(int i, double x, char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;

	if (ch == 'v')
		nc.setzc(i,x);
	if (ch == 'c')
		re.setzc(i,x);
	if (ch =='m')
		me.setzc(i,x);
}


double _setup :: getxc(int i, char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;

	if (ch == 'v')
		return nc.getxc(i);
	if (ch == 'c')
		return re.getxc(i);
	if (ch =='m')
		return me.getxc(i);
	return 0.0;
}

double _setup :: getyc(int i, char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;

	if (ch == 'v')
		return nc.getyc(i);
	if (ch == 'c')
		return re.getyc(i);
	if (ch =='m')
		return me.getyc(i);
	return 0.0;
}

double _setup :: getzc(int i, char ch){
	extern _membrane me;
	extern _nanocarrier nc;
	extern _receptor re;

	if (ch == 'v')
		return nc.getzc(i);
	if (ch == 'c')
		return re.getzc(i);
	if (ch =='m')
		return me.getzc(i);
	return 0.0;
}
