#include "_declarations.h"
#include "_distribution.h"
#include "_setup.h"
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_membrane.h"
#include "_linklist.h"
#include "_datainput.h"
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_selfavoidance.h"

void _distribution :: _distribution :: init () {
  extern _datainput data;
  ostringstream ss;                                                                                            //  passes the filename to print
  binsize_rdf=data.get_binsize_rdf();
  this->L=data.getperiodic_box_length();
  this->H=data.getperiodic_box_height();
  avg_antigen_density=(double)(data.getnum_members('c'))/(double)(pow(L,2));                                   // density of antigens
  verptr=&__module_datastruct_MOD_ver;
  antptr=&__module_datastruct_MOD_antig;

  Nbin_rdf=(int)((L)/binsize_rdf);                                                               // binsize to compute radial distribution function for the antigens
  rdf = new(nothrow) double[Nbin_rdf];
  assign_binvalue(rdf,Nbin_rdf,0.0);
  nsample_rdf=0.0;
}

void _distribution :: reset_box_dimensions(){
  extern _datainput data;
  this->L=data.getperiodic_box_length();
  this->H=data.getperiodic_box_height();
}


void _distribution :: assign_binvalue(double **dist_name, int nbin,  int column,double minval, double maxval, double binsize){    // initialize a column with values between minval and maxval with interval binsize
	int i;
	for (i=0 ;i<nbin;i++)
		dist_name[i][column]=minval+i*binsize;
}

void _distribution :: assign_binvalue(double **dist_name, int nbin,  int column, double value){                         // set all rows of a column to a predefined value zero
	int i;
	for (i=0 ;i<nbin;i++)
		dist_name[i][column]=value;
}

void _distribution :: assign_binvalue(double *dist_name, int nbin, double value){                         // set all rows of a column to a predefined value zero
	int i;
	for (i=0 ;i<nbin;i++)
		dist_name[i]=value;
}

void _distribution :: transfer_binvalues(double **target, double **source, int nbin,  int column){                      // append the value of a column in source array to the same column in target array.
	int i;
	for (i=0 ;i<nbin;i++)
		target[i][column]=target[i][column]+source[i][column];
}



void _distribution :: compute_antigen_rdf_distribution(int nbin, double binsize, double *hist_array,double *nsample){ 
  extern _setup a;
  double distance,rrx,rry,rrz;
  int i,j,binnum;
  int nantigen=a.getnum_members('c');

  *nsample += 1;
  for(i=0;i<nantigen-1;i++){                                                                    // Run over all antigens
    for(j=i+1;j<nantigen;j++){                                                                  // proceed with calculation of rdf with the rest if yes
      rrx=a.getxc(i,'c')-a.getxc(j,'c');
      rry=a.getyc(i,'c')-a.getyc(j,'c');
      rrz=a.getzc(i,'c')-a.getzc(j,'c');
      rrx=rrx-round(rrx/L)*L;                                                                 // periodic boundary along x
      rry=rry-round(rry/L)*L;                                               // periodic boundary along y
      distance=pow((rrx*rrx+rry*rry+rrz*rrz),0.5);
      if (distance < L/2){
       binnum=(int)(distance/binsize);
       if (binnum >= nbin) binnum = nbin-1;
       hist_array[binnum] += 2;
      }
    }
  }
  return;
}

void _distribution :: write_antigen_rdf_distribution(int nbin, double binsize, double *hist_array,double nsample,std::ostringstream &ss){ 
  extern _setup a;
  ofstream outfile;
  outfile.open(ss.str().c_str(), ios::out);                                                                                 // filename passed through string
  int i;
  int nantigen=a.getnum_members('c');
  double memarea = __module_datastruct_MOD_get_membrane_area();
  cout<<"nantigen,memarea "<<nantigen<<" "<<memarea<<endl;
  for(i=0;i<nbin;i++){
    double r = (i+0.5)*binsize;
    double da=(pow((i+1),2)-pow(i,2))*pow(binsize,2);
    double nid = PI*da*nantigen/memarea;
    outfile<<r<<"\t"<<hist_array[i]/(nsample*nantigen*nid)<<"\t"<<hist_array[i]<<endl;                                                    // the factor of 4 is included since we only have 
  }
  outfile.close();
  ss.str();
}

double _distribution :: array_sum(double **hist_array, int nbin, int column){
  int i;
  double sum=0.0;
  for (i=1;i<nbin;i++)
   sum=sum+hist_array[i][column];
  return(sum);
}


void _distribution :: compute_rdf_distribution(){                                                                                                    // distribution restricted to x and y of the NC
   compute_antigen_rdf_distribution(Nbin_rdf,binsize_rdf,rdf,&nsample_rdf);
   cout<<"nsample rdf "<<nsample_rdf<<endl;
}


void _distribution :: write_rdf_distribution(int rank,int frame){
  ostringstream ss;
  ss<<"rdf-ENS-"<<rank<<"_Frame-"<<frame<<".dat";
  write_antigen_rdf_distribution(Nbin_rdf, binsize_rdf, rdf,nsample_rdf,ss);
}


int _distribution :: compute_multivalency(int ves_num){
  extern _nanocarrier nc;
  int ab1;
  int multivalency=0;
  if (nc.getbonded(ves_num)){
   for(ab1=0; ab1<nc.getnum_ab(ves_num); ab1++)
    if (nc.getbonded(ves_num,ab1)) multivalency +=1;
   }
  return multivalency;
}

void _distribution :: check_antigen_vesicle_dist(int ves_num, int stepno, double type_of_move, double type_of_trans,double ratio_m,double flip_move_ratio){
  extern _receptor re;
  extern _datainput data;
  extern _nanocarrier nc;
  int  i;
  for (i=0 ; i< re.getnum_members();i++){
    double cutoff = pow(nc.getsoft_radius(ves_num)+re.getlength(i),2);
    double dist = this->dist_nc_receptor(ves_num,i);
    if(dist<cutoff){
      cout<<"Step "<<stepno<<endl;
      cout<<"Move "<<type_of_move<<" "<<type_of_trans<<" "<<ratio_m<<endl;
      cout<<"Ratio_m "<<ratio_m<<" flip_move_ratio "<<flip_move_ratio<<endl;
      cout<<"Cutoff distance "<<cutoff<<endl;
      cout<<str_green<<"Cutoff violation violated in Antigen "<<i<<" colliding with vesicle "<<ves_num<<" with distance "<<dist<<str_black<<endl;
      cin.get();
    }
  }
}

bool _distribution :: check_bondspring_dist(int mem_sel,char ch){
  extern _datainput data;
  extern _nanocarrier nc;
  extern _receptor re;
  int i,num_members,ant_no;
  int *temp;
  double dist;
  if (ch == 'v'){
    num_members = nc.getnum_ab(mem_sel);
    for (i=0; i<num_members; i++){
     if (nc.getbonded(mem_sel,i)==1){
      ant_no = nc.getbondedto(mem_sel,i);
      dist = this->dist_receptor_ligand(mem_sel,i,ant_no); 
      if ( dist > reaction_dist_sq){
       cout<<str_green<<"Bonded antigen-antibody violating max bond stretch conditions"<<endl;
       cout<<"NC  AB Ant "<<mem_sel<<" "<<i<<" "<<ant_no<<endl;
       cout<<"Computed distance for the bond "<<dist<<endl;
       cout<<"NC coords "<<nc.getxc(mem_sel)<<" "<<nc.getyc(mem_sel)<<" "<<nc.getzc(mem_sel)<<endl;
       cout<<"AB coords "<<nc.getx_antibody(mem_sel,i)<<" "<<nc.gety_antibody(mem_sel,i)<<" "<<nc.getz_antibody(mem_sel,i)<<endl;
       cout<<"AN coords "<<re.getxt(ant_no)<<" "<<re.getyt(ant_no)<<" "<<re.getzt(ant_no)<<str_black<<endl;
       return true;
      }
     }
    }
  }

  else if (ch=='c'){
   if (re.getbonded(mem_sel)==1){
    temp = re.getbondedto(mem_sel);
    dist = this->dist_receptor_ligand(temp[0],temp[1],mem_sel); 
    if ( dist >reaction_dist_sq){
      cout<<str_red<<"Bonded antigen-antibody violating max bond stretch conditions"<<endl;
      cout<<"Ant NC AB "<<mem_sel<<" "<<temp[0]<<" "<<temp[1]<<endl;
      cout<<"Computed distance for the bond "<<dist<<endl;
      cout<<"NC coords "<<nc.getxc(temp[0])<<" "<<nc.getyc(temp[0])<<" "<<nc.getzc(temp[0])<<endl;
      cout<<"AB coords "<<nc.getx_antibody(temp[0],temp[1])<<" "<<nc.gety_antibody(temp[0],temp[1])<<" "<<nc.getz_antibody(temp[0],temp[1])<<endl;
      cout<<"AN coords "<<re.getxt(mem_sel)<<" "<<re.getyt(mem_sel)<<" "<<re.getzt(mem_sel)<<str_black<<endl;
      return true;
    }
   }
  }
  return false;
}

double** _distribution:: compute_histogram(double **array1,double **array2, int nrow, int ncol, double maxval, double minval, double binsize, int nbin){
  double **histogram_data;
  int i,j,binnum;
  histogram_data = new(nothrow) double*[nbin];

  for (i=0; i<nbin; i++)
    histogram_data[i] = new(nothrow) double[ncol+1];                 // declare the array to store the histogram data

  for (i=0; i<ncol+1; i++)
    assign_binvalue(histogram_data,nbin,i,0.0);                     // set all elements to zero

  for (j=0;j<nbin;j++)
    histogram_data[j][0] = minval + j*binsize;                      // the first column contains the values of the x axis   

  for (i=0;i<ncol;i++){
    for(j=0;j<nrow;j++){
      binnum = int((array1[j][i]-array2[j][i]-minval)/binsize);
      if (binnum >=0 && binnum <nbin)
        histogram_data[binnum][i+1] += 1.0;                         // the first column is occupied by the bin values 
    }
  }
  return histogram_data;
}

double** _distribution:: compute_histogram(double **array1,int nrow, int ncol, double maxval, double minval, double binsize, int nbin){
  double **histogram_data;
  int i,j,binnum;
  histogram_data = new(nothrow) double*[nbin];


  for (i=0; i<nbin; i++)
    histogram_data[i] = new(nothrow) double[ncol+1];                 // declare the array to store the histogram data

  for (i=0; i<ncol+1; i++)
    assign_binvalue(histogram_data,nbin,i,0.0);                     // set all elements to zero

  for (j=0;j<nbin;j++)
    histogram_data[j][0] = minval + j*binsize;                      // the first column contains the values of the x axis   

  for (i=0;i<ncol;i++){
    for(j=0;j<nrow;j++){
      binnum = int((array1[j][i]-minval)/binsize);
      if (binnum >=0 && binnum <nbin)
       histogram_data[binnum][i+1] += 1.0;                         // the first column is occupied by the bin values 
    }
  }
  return histogram_data;
}

double _distribution :: dist_receptor_ligand(int ncnum, int abnum, int recep_num){
  extern _nanocarrier nc;
  extern _receptor re;
  //@
  double rrx, rry, rrz;
  rrx = nc.getx_antibody(ncnum, abnum) - re.getxt(recep_num);      // distance between antibody and antigen tip
  rry = nc.gety_antibody(ncnum, abnum) - re.getyt(recep_num);
  rrz = nc.getz_antibody(ncnum, abnum) - re.getzt(recep_num);
  rrx = rrx - this->L*round(rrx/this->L);                          // implements the PBC along x and y
  rry = rry - this->L*round(rry/this->L);
  return (rrx*rrx + rry*rry + rrz*rrz);                            // real distance
}

double _distribution :: dist_nc_receptor(int ncnum, int recep_num){
  extern _nanocarrier nc;
  extern _receptor re;
  //@
  double rrx, rry, rrz;
  rrx = nc.getxc(ncnum) - re.getxt(recep_num);             // distance between antibody and antigen tip
  rry = nc.getyc(ncnum) - re.getyt(recep_num);
  rrz = nc.getzc(ncnum) - re.getzt(recep_num);
  rrx = rrx - this->L*round(rrx/this->L);                 // implements the PBC along x and y
  rry = rry - this->L*round(rry/this->L);
  return (rrx*rrx + rry*rry + rrz*rrz);                   // real distance
}

double _distribution :: dist_nc_nc(int ncnum1, int ncnum2){
  extern _nanocarrier nc;
  //@
  double rrx, rry, rrz;
  rrx = nc.getxc(ncnum1) - nc.getxc(ncnum2);               // distance between antibody and antigen tip
  rry = nc.getyc(ncnum1) - nc.getyc(ncnum2);
  rrz = nc.getzc(ncnum1) - nc.getzc(ncnum2);
  rrx = rrx - this->L*round(rrx/this->L);                  // implements the PBC along x and y
  rry = rry - this->L*round(rry/this->L);
  return (rrx*rrx + rry*rry + rrz*rrz);                    // real distance
}