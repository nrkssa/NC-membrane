#ifndef __RESTART_H__
#define __RESTART_H__
class _restart {
int num_ves, num_ant,  num_ab;
int i,j;
public: _restart(){}
public: void init(){}
public: void write_restart_conf(int ens, int window, int step);
public: void write_restart_conf(int ens, int window);
public: void write_restart_conf(int window);
public: void _write_restart_conf(ofstream &ofile);
public: void read_restart_conf(int procnum);
};
#endif
