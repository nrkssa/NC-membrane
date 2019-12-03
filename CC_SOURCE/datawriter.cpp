#include "_declarations.h"
#include <iomanip>
#include "_datawriter.h"
#include "_datainput.h"
#include "_setup.h"
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_membrane.h"
#include "_fortran_structures.h"
#include "_fortran_modules.h"
#include "_histogram.h"
#include "_distribution.h"

void _datawriter :: init(){
    extern _datainput data;
    extern _nanocarrier nc;
    extern _receptor re;
    extern _membrane me;
    //@
    this->num_ves = nc.getnum_members();
    this->num_ant = re.getnum_members();
    this->nver = me.getnum_vertices();
    this->pbnum = me.getnum_verticestot();
    this->ntr = me.getnum_triangles();
    this->nlink = me.getnum_links();
    //@

    this->num_ab = new (nothrow) int[this->num_ves];

    for (int i=0; i<this->num_ves; i++){
      this->num_ab[i] = nc.getnum_ab(i);
    }

    this->boxl = data.getperiodic_box_length();
    this->boxh = data.getperiodic_box_height();
//    this->eqm_bond_dist=data.geteqm_bond_dist();
    if (boxl<boxh) this->scale = boxl/10;
    if (boxl>=boxh) this->scale = boxh/10;
    this->verptr=&__module_datastruct_MOD_ver;                // Initialize pointers to vertex,triangle and antigen 
    this->triptr=&__module_datastruct_MOD_tri;
    this->antptr=&__module_datastruct_MOD_antig;
    this->ncpointer=&__module_datastruct_MOD_nc_f;
}

_datawriter :: ~_datawriter(){
	outfile2.close();
}

void _datawriter:: reset_box_dimensions(){
    extern _datainput data;
    this->boxl = data.getperiodic_box_length();
    this->boxh = data.getperiodic_box_height();
}

bool _datawriter :: file_exists(stringstream filename){
	ifstream infile(filename.str().c_str());
	return infile.good();
}

string _datawriter :: modify_nc_header(string inp, int j){
	stringstream splitter(inp),ss;
	istream_iterator<string> begin(splitter);
	istream_iterator<string> end;
	vector<string> vstrings(begin, end);
	ss<<vstrings[0]<<" Name= \"NC"<<j<<"\"  "<<vstrings[1]<<"  "<<vstrings[2];
	return ss.str().c_str();
}


void _datawriter:: write_antigen_xyz(int rank,int window,int step){
	extern _nanocarrier nc;
        extern _receptor re;
	extern _histogram hist;
	double *meanr;
	std::ostringstream ss;
	ss<<"../SYS_STATE/DUMP/antigen-"<<rank<<"-"<<window<<"-"<<step<<".dat";
	std::ofstream fs;
	fs.open(ss.str().c_str(),ios::out);

	fs<<"#boxlen:  "<<(float)this->boxl<<"  area:  "<<(float)__module_datastruct_MOD_get_membrane_area()<<endl;
	for (int i=0;i<num_ves;i++){
	  fs<<"#ncpos: "<<i<<" "<<nc.getxc(i)<<" "<<nc.getyc(i)<<" "<<nc.getzc(i)<<endl;
	  meanr = nc.get_shadowr(i);
    	  fs<<"#shadowpos: "<<meanr[0]<<"\t"<<meanr[1]<<"\t"<<meanr[2]<<" zref: "<<nc.get_biasref(i)<<"\n";
	}

	for (int j=0; j<this->num_ant; j++){
	int vshadow = (antptr+j)->vertex-1;
	fs<<re.getxc(j)<<" "<<re.getyc(j)<<" "<<re.getzc(j)<<" "<<re.getbonded(j)<<" "<<(this->verptr+vshadow)->shadownc<<endl;
	}
	fs.close();
	ss.str("");
}

void _datawriter :: VTK_system(int ensemb,int window,int state_num)
{
	extern _datainput data;
	stringstream ss1;
	ofstream of1;
	int j=0;

	if (state_num >=0)
		ss1<<"../SYS_STATE/SYSTEM/system-ENS-"<<ensemb<<"-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
		ss1<<"../SYS_STATE/SYSTEM/system-ENS-"<<ensemb<<"-window-"<<window<<".vtu";

	of1.open(ss1.str().c_str(),ios::out);
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;

	for (j=0; j<this->num_ves; j++)
		 _VTK_NanoCarrier(of1,j);
	 _VTK_membrane(of1);
	 _VTK_Antigens(of1);
	 of1 <<"</UnstructuredGrid>"<<endl;
	 of1 <<"</VTKFile>"<<endl;
	 of1.close();
}

void _datawriter :: VTK_system(int window,int state_num)
{
	extern _datainput data;
	stringstream ss1;
	ofstream of1;
	int j=0;
	
	if (state_num >=0)
		ss1<<"../SYS_STATE/SYSTEM/system-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
		ss1<<"../SYS_STATE/SYSTEM/system-window-"<<window<<".vtu";
	
	of1.open(ss1.str().c_str(),ios::out);
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	
	for (j=0; j<this->num_ves; j++)
		_VTK_NanoCarrier(of1,j);
	_VTK_membrane(of1);
	_VTK_Antigens(of1);
	of1 <<"</UnstructuredGrid>"<<endl;
	of1 <<"</VTKFile>"<<endl;
	of1.close();
}

void _datawriter :: VTK_membrane(int ensemb, int window,int state_num)
{
	stringstream ss1;
	ofstream of1;
	if (state_num >=0)
	 ss1<<"../SYS_STATE/MEMB/membrane-ENS-"<<ensemb<<"-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
	 ss1<<"../SYS_STATE/MEMB/membrane-ENS-"<<ensemb<<"-window-"<<window<<".vtu";

	of1.open(ss1.str().c_str(),ios::out);
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	_VTK_membrane(of1);
	of1 <<"</UnstructuredGrid>"<<endl;
	of1 <<"</VTKFile>"<<endl;
	of1.close();
}

void _datawriter :: VTK_membrane(int window,int state_num)
{
	stringstream ss1;
	ofstream of1;
	if (state_num >=0)
		ss1<<"../SYS_STATE/MEMB/membrane-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
		ss1<<"../SYS_STATE/MEMB/membrane-window-"<<window<<".vtu";
	
	of1.open(ss1.str().c_str(),ios::out);
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	_VTK_membrane(of1);
	of1 <<"</UnstructuredGrid>"<<endl;
	of1 <<"</VTKFile>"<<endl;
	of1.close();
}


void _datawriter :: VTK_NanoCarrier(int ensemb,int window,int state_num)
{
	extern _datainput data;
	stringstream ss1;
	ofstream of1;

	if (state_num >=0)
	 ss1<<"../SYS_STATE/NC-DATA/NC-ENS-"<<ensemb<<"-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
	 ss1<<"../SYS_STATE/NC-DATA/NC-ENS-"<<ensemb<<"-window-"<<window<<".vtu";

	of1.open(ss1.str().c_str(),ios::out);
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	for (int j=0; j<this->num_ves; j++)
	 	_VTK_NanoCarrier(of1,j);

	 of1 <<"</UnstructuredGrid>"<<endl;
	 of1 <<"</VTKFile>"<<endl;
	 of1.close();
}


void _datawriter :: VTK_Antigens(int ensemb,int window,int state_num)
{
	stringstream ss1;
	ofstream of1;

	if (state_num >=0)
		ss1<<"../SYS_STATE/ANTIGEN/Antigen-ENS-"<<ensemb<<"-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
		ss1<<"../SYS_STATE/ANTIGEN/Antigen-ENS-"<<ensemb<<"-window-"<<window<<".vtu";

	of1.open(ss1.str().c_str(),ios::out); 
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	_VTK_Antigens(of1);
	of1 <<"</UnstructuredGrid>"<<endl;
	of1 <<"</VTKFile>"<<endl;

}

void _datawriter :: VTK_NanoCarrier(int window,int state_num)
{
	extern _datainput data;
	stringstream ss1;
	ofstream of1;
	
	if (state_num >=0)
		ss1<<"../SYS_STATE/NC-DATA/NC-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
		ss1<<"../SYS_STATE/NC-DATA/NC-window-"<<window<<".vtu";
	
	of1.open(ss1.str().c_str(),ios::out);
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	for (int j=0; j<this->num_ves; j++)
		_VTK_NanoCarrier(of1,j);
	
	of1 <<"</UnstructuredGrid>"<<endl;
	of1 <<"</VTKFile>"<<endl;
	of1.close();
}

void _datawriter :: VTK_Antigens(int window,int state_num)
{
	stringstream ss1;
	ofstream of1;
	
	if (state_num >=0)
		ss1<<"../SYS_STATE/ANTIGEN/Antigen-window-"<<window<<"-state-"<<state_num<<".vtu";
	else
		ss1<<"../SYS_STATE/ANTIGEN/Antigen-window-"<<window<<".vtu";
	
	of1.open(ss1.str().c_str(),ios::out); 
	of1<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">"<<endl;
	of1<<"<UnstructuredGrid>"<<endl;
	_VTK_Antigens(of1);
	of1 <<"</UnstructuredGrid>"<<endl;
	of1 <<"</VTKFile>"<<endl;
	
}


void _datawriter :: _VTK_Antigens(ofstream &of1)
{
	int i,N_tri;
	extern _receptor re;
	extern _datainput data;
	extern _distribution distrib;
	int* temp;
	double ant1,ant2,ant3;
	N_tri=0;
	
	of1<<"<Piece Name= \"receptor\"  NumberOfPoints= \""<<this->num_ant<<"\"  NumberOfCells=\""<< N_tri <<"\">"<<endl;
	of1<<"<PointData Scalars=\"scalars\">"<<endl;
        of1<<"<DataArray type=\"Int32\" Name=\"Antibody_number\" Format=\"ascii\">" <<endl;
	
	  for(i=0; i< this->num_ant; i++){
	  if (re.getbonded(i) == 1){
	   temp = re.getbondedto(i);
	   of1<<temp[1]<<endl;
	  }
	  
	   else{
	    of1<<"-1"<<endl;
	   }
	  }
	  of1<<"</DataArray>" <<endl;
	  
	  of1<<"<DataArray type=\"Int32\" Name=\"bonded\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
	   of1<<re.getbonded(i)<<endl;
	  }
	  of1<<"</DataArray>" <<endl;
	  
	  of1<<"<DataArray type=\"Float32\" Name=\"theta\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
		of1<<(antptr+i)->theta<<endl;
  	  }
	  of1<<"</DataArray>" <<endl;
	  
          of1<<"<DataArray type=\"Float32\" Name=\"phi\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
	   of1<<(antptr+i)->phi<<endl;
  	  }
	  of1<<"</DataArray>" <<endl;
	  
	  of1<<"<DataArray type=\"Int32\" Name=\"type\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
	   of1<<(antptr+i)->ant_type<<endl;
	  }
	  of1<<"</DataArray>" <<endl;
	  
	  of1<<"<DataArray type=\"Float32\" Name=\"length\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
	   of1<<(antptr+i)->length<<endl;
	  }
	  of1<<"</DataArray>" <<endl;
	  
	  of1<<"<DataArray type=\"Float32\" Name=\"radius\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
	   of1<<(antptr+i)->radius<<endl;
	  }
	  of1<<"</DataArray>" <<endl;
	 
	  of1<<"<DataArray type=\"Int32\" Name=\"vertex\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
	   of1<<(antptr+i)->vertex<<endl;
	  }
	  of1<<"</DataArray>" <<endl;
	  
 	  of1<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Antigen_vec\" Format=\"ascii\">" <<endl;
	  for(i=0; i<this->num_ant; i++){
  	  ant1 = re.getxt(i)-re.getxc(i); 
  	  ant2 = re.getyt(i)-re.getyc(i); 
  	  ant3 = re.getzt(i)-re.getzc(i); 
	  ant1 = ant1-round(ant1/boxl)*boxl;                                                                                // Periodic Box conditions for the antigens
	  ant2 = ant2-round(ant2/boxl)*boxl;
	  of1<<ant1<<"\t"<<ant2<<"\t"<<ant3<<endl;
	  }
	  of1<<"</DataArray>" <<endl;
	  
        of1<<"</PointData>"<<endl;
	
        of1<<"<Points>"<<endl;
	of1<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;
	for(i=0; i<this->num_ant; i++){
	of1<<re.getxc(i)<<'\t'<<re.getyc(i)<< '\t'<<re.getzc(i)<<endl;
	}
	of1<<"</DataArray>"<<endl;
        of1<<"</Points>"<<endl;
	
	of1<<"<Cells>"<<endl;
	of1<<"<DataArray type=\"Int32\"  Name=\"connectivity\" Format=\"ascii\">"<<endl;
	of1<<"</DataArray>"<<endl;
	
	of1<<"<DataArray type=\"Int32\" Name=\"offsets\"  Format=\"ascii\">"<<endl;
	of1<<"</DataArray>"<<endl;
	
	of1<<"<DataArray type=\"Int32\" Name=\"types\"  Format=\"ascii\">"<<endl;
	of1<<"</DataArray>"<<endl;
	of1<<"</Cells>"<<endl; 
	of1<<"</Piece>"<<endl;
}

void _datawriter :: _VTK_membrane(ofstream &of1)
{
	int i;
	extern _histogram hist;
	extern _datainput data;
	stringstream ss1;


	for(i=0; i<this->num_ves; i++){                                             //get the center of mass position of the membrane under the shadow of nanocarrier i
		int j = i+1;
		__module_datastruct_MOD_ncshadow_r(&j);
	}

	of1<<"<Piece Name= \"membrane\" NumberOfPoints= \""<<this->pbnum+this->num_ves<<"\"  NumberOfCells=\""<< this->ntr <<"\">"<<endl;
	of1<<"<PointData Scalars=\"scalars\">"<<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"Boundary\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->nver;i++)
	{
	of1<<"0"<<endl;
	}
	for(i=this->nver;i<this->pbnum;i++)
	{
	of1<<"1"<<endl;
	}
	for(i=0;i<this->num_ves;i++)
		of1<<100+i<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Float32\" Name=\"H\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->mcur<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0.0"<<endl;
	of1<<"</DataArray>" <<endl;


	of1<<"<DataArray type=\"Int32\" Name=\"nantigen\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->nantigen<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0"<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"boxvert\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->boxvert<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0"<<endl;
	of1<<"</DataArray>" <<endl;

	for(j=0; j<this->num_ves; j++){
	of1<<"<DataArray type=\"Int32\" Name=\"ncshadow"<<j+1<<"\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->shadownc[j]<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<-1<<endl;
	of1<<"</DataArray>" <<endl;
	}

	of1<<"<DataArray type=\"Float32\" Name=\"delhref\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
          of1<<(verptr+i)->mcur<<endl;
	for(i=0;i<this->num_ves;i++)
	  of1<<"0.0"<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Float32\" Name=\"c1\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->cur1<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0.0"<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Float32\" Name=\"c2\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->cur2<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0.0"<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"czero-flag\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->czero_flag<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0"<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Float32\" Name=\"czero\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->czero<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0.0"<<endl;
	of1<<"</DataArray>" <<endl;


	of1<<"<DataArray type=\"Int32\" Name=\"cellno\" Format=\"ascii\">" <<endl;
	for(i=0;i<this->pbnum;i++)
		of1<<(verptr+i)->cellno<<endl;
	for(i=0;i<this->num_ves;i++)
		of1<<"0"<<endl;
	of1<<"</DataArray>" <<endl;

	of1<<"</PointData>"<<endl;


	of1<<"<CellData Scalars=\"scalars\">"<<endl;
	of1<<"<DataArray type=\"Int32\" Name=\"face-pbflag\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
		 of1<<(triptr+i)->pbflag<<endl;
	of1<<"</DataArray>"<<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"face-nantigen\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
		 of1<<(triptr+i)->nantigen<<endl;
	of1<<"</DataArray>"<<endl;

	of1<<"<DataArray type=\"Float32\" Name=\"face-vol\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
		 of1<<(triptr+i)->vol<<endl;
	of1<<"</DataArray>"<<endl;

	of1<<"<DataArray type=\"Float32\" Name=\"face-area\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
		 of1<<(triptr+i)->ar<<endl;
	of1<<"</DataArray>"<<endl;


	of1<<"<DataArray type=\"Int32\" Name=\"boxtri\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
		 of1<<(triptr+i)->boxtriangle<<endl;
	of1<<"</DataArray>"<<endl;

	
	of1<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\"  Name=\"FaceNormal\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
		 of1<<(triptr+i)->fnor[0][0]<<'\t'<<(triptr+i)->fnor[1][0]<<'\t'<<(triptr+i)->fnor[2][0]<<endl;
	of1<<"</DataArray>"<<endl;
	of1<<"</CellData>"<<endl;
		
	
	of1<<"<Points>"<<endl;
	of1<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;
	for(i=0;i<this->pbnum;i++)
	{of1<<(verptr+i)->vcoord[0][0]<<'\t'<<(verptr+i)->vcoord[1][0] << '\t'<<(verptr+i)->vcoord[2][0]<<'\n';
	}
	for(i=0;i<this->num_ves;i++)
		of1<<(ncpointer+i)->meanr[0][0]<<"\t"<<(ncpointer+i)->meanr[1][0]<<"\t"<<(ncpointer+i)->meanr[2][0]<<"\n";
	of1<<"</DataArray>"<<endl;
	of1<<"</Points>"<<endl;

	of1<<"<Cells>"<<endl;
	of1<<"<DataArray type=\"Int32\"  Name=\"connectivity\" Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
	{
	 of1<<(triptr+i)->vert[0]-1<<'\t'<<(triptr+i)->vert[1]-1<<'\t'<<(triptr+i)->vert[2]-1<<endl;
	}
	of1<<"</DataArray>"<<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"offsets\"  Format=\"ascii\">"<<endl;
	for (i=1;i<=this->ntr;i++)
	{
	of1<<i*3<<'\n';
	}
	of1<<"</DataArray>"<<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"types\"  Format=\"ascii\">"<<endl;
	for (i=0;i<this->ntr;i++)
	{
	of1<<"5"<<'\n';
	}
			of1<<"</DataArray>"<<endl;
	of1<<"</Cells>"<<endl; 
	of1<<"</Piece>"<<endl;
}

void _datawriter :: _VTK_NanoCarrier(ofstream &of1, int j)
{
	int i;
	extern _datainput data;
	extern _nanocarrier nc;

	string field,splitword,empty;
	stringstream trfile,headerfile;
				
	nc.setxy_npb_antibody(j);
	trfile <<"../PARAMETERS/SPH_RAD1_"<<this->num_ab[j]<<"/polygon-connectivity-offset.vtu";
	headerfile <<"../PARAMETERS/SPH_RAD1_"<<this->num_ab[j]<<"/vtu-header.vtu";
	ifstream trdata (trfile.str().c_str(), ios::in);
	ifstream headdata (headerfile.str().c_str(), ios::in);

	if (trdata.good() && headdata.good()){

	while (!headdata.eof()){                                             		// Read triangulation,connectivity and offset data from a file
	 getline(headdata,field,'\0');
	 of1<<modify_nc_header(field,j)<<endl;
	}
	headdata.close();

	of1<<"<PointData Scalars=\"scalars\">"<<endl;
	

	of1<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\"  Name=\"Antibody\" Format=\"ascii\">" <<endl;               //Antibody Data
	for(i=0;i<this->num_ab[j];i++)
	{
	of1<<nc.getx_npb_antibody(j,i)-nc.getxf_npb_antibody(j,i)<<" "<<nc.gety_npb_antibody(j,i)-
	     nc.getyf_npb_antibody(j,i)<<" "<<nc.getz_antibody(j,i)-nc.getzf_antibody(j,i)<<endl;
	}
	of1<<"</DataArray>" <<endl;

	of1<<"<DataArray type=\"Int32\" Name=\"Bonded-State\" Format=\"ascii\">" <<endl;                                       // Bonded State
	for(i=0;i<this->num_ab[j];i++)
	{
	of1<<nc.getbonded(j,i)<<endl;
	}
	of1<<"</DataArray>" <<endl;
	
	
	int *ab_type=nc.get_ab_types(j);
	of1<<"<DataArray type=\"Int32\" Name=\"AB-TYPE\" Format=\"ascii\">" <<endl;                                       // Bonded State
	for(i=0;i<this->num_ab[j];i++)
		of1<<ab_type[i]<<endl;
	of1<<"</DataArray>" <<endl;
	

	of1<<"<DataArray type=\"Int32\" Name=\"Antibody-number\" Format=\"ascii\">" <<endl;                                       // Bonded State
	for(i=0;i<this->num_ab[j];i++)
		of1<<i<<endl;
	of1<<"</DataArray>" <<endl;
	
	

	of1<<"<DataArray type=\"Int32\" Name=\"Bonded-Antigen\" Format=\"ascii\">" <<endl;                                     // Antigen Bonded to
	for(i=0;i<this->num_ab[j];i++)
	{
	of1<<nc.getbondedto(j,i)<<endl;
	}
	of1<<"</DataArray>" <<endl;
	of1<<"</PointData>"<<endl;
	
	of1<<"<Points>"<<endl;
	of1<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;
	for(i=0;i<this->num_ab[j];i++){
	of1<<nc.getxf_npb_antibody(j,i)<<'\t'<<nc.getyf_npb_antibody(j,i) << '\t'<<nc.getzf_antibody(j,i)<<endl;
	}
	of1<<"</DataArray>"<<endl;
	of1<<"</Points>"<<endl;

	while (!trdata.eof()){                                                          // Read triangulation,connectivity and offset data from a file
	 getline(trdata,field,'\0');
	 of1<<field<<endl;
	}
	trdata.close();
	of1<<"</Piece>"<<endl;
	headerfile.str("");
	trfile.str("");
	}
}
