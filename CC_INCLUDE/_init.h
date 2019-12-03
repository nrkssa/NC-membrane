#ifndef __INIT_H__
#define __INIT_H__
class _init{
int num_members, num_ab;
double radius, soft_radius, L, H, chem_cutoffu, ant_size,chem_cutoffl;
double overlap_ves,overlap_mem;
int k, j,gridsize;
double xt,yt;
struct vertex *verptr;
struct antigen *antptr;
char ch;
double z_start;

public: _init(){}
public: ~_init(){}
public: void init(int mode);
public: void reset_box_dimensions();
public: void random_distribution();
public: void place_at_nearest_location();
};
#endif
