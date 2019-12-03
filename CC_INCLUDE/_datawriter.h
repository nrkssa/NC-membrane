class _datawriter {
struct vertex *verptr;
struct triangle *triptr;
struct antigen *antptr;
struct nanocarrier *ncpointer;
int num_ves, num_ant, *num_ab, num_memb, nver, ntr, nlink, pbnum;
int i,j,k;
double * abcd;
double scale; 				                        	//to scale the system to proper units. If box is in mm, and our units are nm, we get huge numbers and vmd cannot handle it
double boxl, boxh,eqm_bond_dist;
ofstream outfile1, outfile2;
public: _datawriter(){}
public: void init();
public: ~_datawriter();
public: void reset_box_dimensions();
public: void write_antigen_xyz(int rank,int window,int step);
public:void VTK_system(int window, int state_num=-1);
public:void VTK_system(int ensemb,int window, int state_num);
public:void VTK_membrane(int window, int state_num=-1);
public:void VTK_membrane(int ensemb,int window, int state_num);
public:void _VTK_membrane(ofstream &of1);
public:void VTK_NanoCarrier(int window, int state_num=-1);
public:void VTK_NanoCarrier(int ensemb,int window, int state_num);
public:void _VTK_NanoCarrier(ofstream &of1, int j);
public:void VTK_Antigens(int window, int state_num=-1);
public:void VTK_Antigens(int ensemb,int window, int state_num);
public:void _VTK_Antigens(ofstream &of1);
public:bool file_exists(stringstream filename);
public:string modify_nc_header(string inp, int j);
};
