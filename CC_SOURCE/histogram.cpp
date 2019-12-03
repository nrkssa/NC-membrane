#include "_declarations.h"
#include "_histogram.h"

void _histogram :: init(){
}

double _histogram :: compute_arraysum(double **data, int datasize, int ncol){
 double sum=0.0;
 for (int i=0 ; i<datasize; i++){
  sum += data[i][ncol];
 }
 return (sum/datasize);
}

double _histogram:: find_max_entry(double **array, int nrow, int ncol){
  double maxdata=0.0;
  for (int j=0; j<ncol; j++){
   for  (int i=0 ; i<nrow ; i++){
    if (i==0 && j==0) 
     maxdata = array[i][j];
    else{
      if (array[i][j] > maxdata)
       maxdata=array[i][j];
    }
   }
  }
  return maxdata;
}

double _histogram:: find_max_entry(double **array1, double **array2, int nrow, int ncol){
  double maxdata=0.0;
  for (int j=0; j<ncol; j++){
    for  (int i=0 ; i<nrow ; i++){
      if (i==0 && j==0)
        maxdata = array1[i][j]-array2[i][j];
      else{
        if ((array1[i][j]-array2[i][j]) > maxdata)
          maxdata=array1[i][j]-array2[i][j];
      }
    }
  }
  return maxdata;
}

double _histogram:: find_min_entry(double **array, int nrow, int ncol){
  double mindata=0.0;
  for (int j=0; j<ncol; j++){
    for  (int i=0 ; i<nrow ; i++){
      if (i==0 && j==0) 
        mindata = array[i][j];
      else{
        if (array[i][j] < mindata)
          mindata=array[i][j];
      }
    }
  }
  return mindata;
}

double _histogram:: find_min_entry(double **array1, double **array2, int nrow, int ncol){
  double mindata=0.0;
  for (int j=0; j<ncol; j++){
    for  (int i=0 ; i<nrow ; i++){
      if (i==0 && j==0)
        mindata = array1[i][j]-array2[i][j];
      else{
        if ((array1[i][j]-array2[i][j]) < mindata)
          mindata=array1[i][j]-array2[i][j];
      }
    }
  }
  return mindata;
}