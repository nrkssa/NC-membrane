#ifndef __MEMBRANE_H__
#define __MEMBRANE_H__

class _MEMBRANE{
int nver,pbnum,ntr,nlink;
struct vertex *verptr;
struct antigen *antptr;
public: void init();
public: _MEMBRANE(){};
public: _MEMBRANE(int nver1, int pbnum1, int ntr1, int nlink1);
public: int getnum_vertices();
public: int getnum_triangles();
public: int getnum_links();
public: int getnum_verticestot();
public: double getxc(int i);
public: double getyc(int i);
public: double getzc(int i);
public: void setxc(int i,double rx);
public: void setyc(int i,double ry);
public: void setzc(int i,double rz);
public: void setcellno(int i,int cellno);
};

class _membrane {
struct vertex *verptr;
struct antigen *antptr;
int num_memb_vert;
double L, H;
int i;
 _MEMBRANE flmemb;
public: _membrane(){}
public: void init();
public: _membrane(double L1, double H1, int nver1, int pbnum1, int ntr1, int nlink1);
public: void reset_box_dimensions();
public: int getnum_vertices();
public: int getnum_triangles();
public: int getnum_links();
public: int getnum_verticestot();
public: double getL();
public: double getH();
public: void setxc(int i, double x);
public: double getxc(int i);
public: void setyc(int i, double y);
public: double getyc(int i);
public: void setzc(int i, double z);
public: double getzc(int i);
public: void setcellno(int i, int cellno);
public: void movevertex(bool eqmflag=false);
public: void fliplink(bool eqmflag=false);
public: void membrane_montecarlo(bool eqmflag=false);
public: void dump_membrane_conf(int ensno);
public: void dump_membrane_conf(int ensno, int window);
};
#endif
