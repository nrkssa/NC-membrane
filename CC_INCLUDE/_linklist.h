#ifndef __LINKLIST_H__
#define __LINKLIST_H__

class linkedlist{
int *lscl, *head,*cells_zdir;	          // array of lists and array of head cells
double L, H, cellIxy, cellIz;                     // length of periodic box, cellIxy is inverse of length per linking cell in x and y dirn, cellIz is in z dirn
int icell,i, k;
int ncell,num_members, Mxy, Mz;                   // Mxy is number of cell in x or y dimension, Mz is number of link cell in z-dirn; ncell is total # of linking cells,
int cell[27],cell2ring[125],cell1ring[27];        // nearest cells
char ch;

public: linkedlist(){}
public: void init(char ch1, double cellsize);
public: void reset_box_dimensions();
public: void calc();
public: void calc(int sel_mem, double x, double y, double z);
public: void locate_allcells(int obj_no);
public: int* getlscl();
public: int* gethead();
public: void printheadlist();
public: void printlscllist();
public: double getcellIxy ();
public: double getcellIz ();
public: int getMxy();
public: int getMz();
public: int* getcells(double x, double y, double z);
public: int getcell(double x, double y, double z);                      //Returns the cell to which the vesicle belongs
public: int getcellnumber(double x, double y, double z);                      //Returns the cell to which the vesicle belongs
public: int getcellnumber(int i, char ch);
public: int* get2ringcells(int i,char ch);
public: int* get1ringcells(int i,char ch);
public: int* getcells_zdirection(int i,char ch);
public: void is_member_present_in_cell(int sel_mem);
public: void print_celllist(int sel_cell);
};

class _linklist {
linkedlist veslink,antigenlink,membvertlink;
int errorflag[1],derrorflag[1];
public: _linklist() {}
public: void init();
public: void reset_box_dimensions();
public: void init(char ch);
public: void calc(char ch);
public: void calc(int sel_mem, double x, double y, double z, char ch); 
public: void locate_allcells(int obj_no,char ch);
public: int* getlscl(char ch); 
public: int* gethead(char ch); 
public: void printheadlist(char ch);
public: void printlscllist(char ch);
public: double getcellIxy(char ch);
public: double getcellIz(char ch);
public: int getMxy(char ch);
public: int getMz(char ch);
public: int* getcells(double x, double y, double z, char ch); 
public: int getcell(double x, double y, double z,char ch);              //Returns the cell to which the vesicle belongs
public: int getcellnumber(int i, char ch);
public: int* get2ringcells(int i,char ch);
public: int* get1ringcells(int i,char ch);
public: int* getcells_zdirection(int i,char ch);
public: void is_member_present_in_cell(int sel_mem,char ch);
public: void print_celllist(int sel_cell,char ch);
public: void throw_errormessage(string input,char ch);
};
#endif
