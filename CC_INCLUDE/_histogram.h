#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__
class _histogram {
double delr, kz, delh, kh, rref,href;						// for umbrella sampling
double* kbias;
double temperature;
string approach_dir;

public: _histogram(){}
public: void init();
public: void set_kzkh_fortran(double kzdir,double kcurv);
public: void set_kz_fortran(double kzdir);
public: void set_kh_fortran(double kcurv);
public: void set_rref_fortran(double currrref);
public: void set_href_fortran(double currhref);
public: void set_zhref_fortran(double currrref,double currhref);
public: void advance_rref(string direction);
public: void recede_rref();
public: void advance_rref();
public: double getrref();
public: double gethref();
public: double getdelr();
public: double getdelh();
public: double getkz();
public: double getkh();
public: double* get_kbias();
public: double compute_arraysum(double **data,int datasize, int ncol);
public: double find_max_entry(double **array, int nrow, int ncol);
public: double find_max_entry(double **array1, double **array2, int nrow, int ncol);
public: double find_min_entry(double **array, int nrow, int ncol);
public: double find_min_entry(double **array1, double **array2, int nrow, int ncol);
};
#endif
