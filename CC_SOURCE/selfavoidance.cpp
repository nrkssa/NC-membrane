#include "_declarations.h"
#include "_selfavoidance.h"
#include "_datainput.h"               // definition for data reading from input.in
#include "_setup.h"                   // Lays down the datastructure for system variables
#include "_linklist.h"                // Construct a linklist connecting space to nano carrier and memb.
#include "_selfavoidance.h"
#include "_fortran_structures.h"
#include "_nanocarrier.h"
#include "_receptor.h"
#include "_membrane.h"

void _selfavoidance :: init(){
	extern _datainput data;
	extern _nanocarrier nc;
	extern _receptor re;
	int i,j;
	
	this->L = data.getperiodic_box_length();
	this->H = data.getperiodic_box_height();
	this->num_nc = nc.getnum_members();
    this->soft_radius = new (nothrow) double[num_nc];  //soft radius for each vesicle
    this->radius = new (nothrow) double[num_nc];       //soft radius for each vesicle
	
	for (i=0; i<this->num_nc; i++){                    // cutoff distance for each NC-NC pair
	  this->soft_radius[i] = nc.getsoft_radius(i);
	  this->radius[i] = nc.getradius(i);
	}
	
	this->vesicle_vesicle_cutoff = new(nothrow) double*[num_nc];
	for (i=0; i<this->num_nc; i++){
	this->vesicle_vesicle_cutoff[i] = new(nothrow) double[num_nc];
	}
	
	// Distance between various NCs
	for (i=0;i<this->num_nc; i++){
	  for (j=0;j<this->num_nc; j++){
	    vesicle_vesicle_cutoff[i][j] = pow(this->soft_radius[i]+this->soft_radius[j],2);  //1 (v-v)
	  }
	}
	
	// Distance between NC and _membrane
	this->vesicle_membrane_cutoff = new(nothrow) double[num_nc];
	for (i=0; i<num_nc; i++){
	  this->vesicle_membrane_cutoff[i] = pow(this->soft_radius[i],2);                     //2 (v-m)
	}
	
	// Distance between various num_antigens
	this->num_antigens = re.getnum_members();
	this->ant_radius = new (nothrow) double[this->num_antigens];                        //soft radius for each vesicle
	for (i=0; i<this->num_antigens; i++){                                               // cutoff distance for each NC-NC pair
	  this->ant_radius[i] = re.getradius(i);
	}
	
	// Distance between antigens and vesicles
	this->antigen_vesicle_cutoff_dist = new(nothrow) double*[this->num_nc];
	for (i=0; i<this->num_nc; i++){
	this->antigen_vesicle_cutoff_dist[i] = new(nothrow) double[this->num_antigens];
	}
	
	for (i=0;i<this->num_nc; i++){
	  for (j=0;j<this->num_antigens; j++){
	    antigen_vesicle_cutoff_dist[i][j] = pow(this->radius[i],2);                      //3 (v-c)
	  }
	}
	
	this->antigen_antigen_cutoff = new(nothrow) double*[this->num_antigens];
	for (i=0; i<this->num_antigens; i++){
	this->antigen_antigen_cutoff[i] = new(nothrow) double[this->num_antigens];
	}
	
	for (i=0;i<this->num_antigens; i++){
	  for (j=0;j<this->num_antigens; j++){
	    this->antigen_antigen_cutoff[i][j] = pow(this->ant_radius[i]+this->ant_radius[j],2);  //4 (c-c)
	  }
	}
	
	this->verptr=&__module_datastruct_MOD_ver;
}

void _selfavoidance :: reset_box_dimensions(){
	extern _datainput data;
	this->L = data.getperiodic_box_length();
	this->H = data.getperiodic_box_height();
}

double _selfavoidance :: sa_distance(char ch1, char ch2, int m1, int m2){
      double cutoff=0.0;
      if ((ch1 == 'c') && (ch2 == 'c')){
         cutoff = this->antigen_antigen_cutoff[m1][m2];
      }
      else if (ch1 == 'v' && ch2 == 'v'){
         cutoff = this->vesicle_vesicle_cutoff[m1][m2];
      }
      else if((ch1 == 'c') && (ch2 == 'v')){
         cutoff = this->antigen_vesicle_cutoff_dist[m2][m1];
      }
      else if((ch1 == 'v') && (ch2 == 'c')){
         cutoff = this->antigen_vesicle_cutoff_dist[m1][m2];
      }
      else if((ch1 == 'v') && (ch2 == 'm')){
         cutoff= this->vesicle_membrane_cutoff[m1];
      }
      else if ((ch1=='m') && (ch2=='v')){
         cutoff = this->vesicle_membrane_cutoff[m2];
      }
      else{
       cout<<"Invalid arguments "<<ch1<<" and "<<ch2<< "to function "<<__FUNCTION__<<endl;
       exit(1);
      }
      return cutoff;
}

//@  vesicle-vesicle or receptor-receptor  
bool _selfavoidance ::  overlap (int mem_sel, char ch){
    extern _setup a;
    extern _linklist link1;
    unsigned int i;
    int mem, *lscl, *head;
    double cutoff, rrx, rry, rrz, dist;
	
    tempx = a.getxc(mem_sel, ch);
    tempy = a.getyc(mem_sel, ch);
    tempz = a.getzc(mem_sel, ch);
	
	if (ch=='v'){
  	  for(mem=0; mem<this->num_nc; mem++){
		if (mem != mem_sel){
		 rrx = tempx - a.getxc(mem, ch);
		 rry = tempy - a.getyc(mem, ch);
		 rrz = tempz - a.getzc(mem, ch);
		 rrx = rrx - this->L*round(rrx/this->L);             //calculate minimum image distance
		 rry = rry - this->L*round(rry/this->L);
		 dist = (rrx*rrx + rry*rry + rrz*rrz);              // distance is square of actual distance
		 cutoff = this->sa_distance(ch,ch,mem_sel,mem);
		 if (dist < cutoff) return true;
		}
	  }
	}
    
    else{
     lscl = link1.getlscl(ch);                                                                 // pointers to the linklist lscl   
     head = link1.gethead(ch);                                                                 // pointers to the header cell linklist
     cell = link1.getcells(tempx, tempy, tempz,ch);                                        // All 26 neighbours are returned here
     for (i=0;i<27;i++){                                                                  //To check for the overlap: vesicle-vesicle or antigen-antigen
      mem = head[cell[i]];                                                               // get the number of nano carrier in cell i
	  if (mem == -1) continue;                                                           // Check for occupied cells
	  do {
	   if (mem == mem_sel) {mem = lscl[mem]; continue;}                           // skip calculation with itself       
	    rrx = tempx - a.getxc(mem, ch);
	    rry = tempy - a.getyc(mem, ch);
	    rrz = tempz - a.getzc(mem, ch);
	    rrx = rrx - this->L*round(rrx/this->L);             //calculate minimum image distance
	    rry = rry - this->L*round(rry/this->L);
	    dist = (rrx*rrx + rry*rry + rrz*rrz);              // distance is square of actual distance
        cutoff = this->sa_distance(ch,ch,mem_sel,mem);
	    if (dist < cutoff) return true;                                  // two objects cannot come closer than their diamter
	    mem = lscl[mem];
	  } while (mem != -1);
	 }
	}
	return false;
}

//@ Any kind
bool _selfavoidance ::  overlap (int mem_sel, char ch, char ch1){         // Checking overlap of one given object with another class of objects
    extern _nanocarrier nc;
    extern _membrane me;
    extern _receptor re;
    extern _datainput data;
    extern _linklist link1;
    
    int i, mem, *lscl, *head, *cell;
    double cutoff, rrx, rry, rrz, dist;
        
    if (ch == ch1){
      return this->overlap(mem_sel, ch);                                                     // same as the above function
    }
    
    else if(ch=='v'){
       tempx = nc.getxc(mem_sel);     
       tempy = nc.getyc(mem_sel);     
       tempz = nc.getzc(mem_sel);     

	  if( ch1 == 'c' ){
	   for (i=0 ; i< num_antigens; i++){
	   rrx = tempx - re.getxt(i);
	   rry = tempy - re.getyt(i);
	   rrz = tempz - re.getzt(i);
	   rrx = rrx - this->L*round(rrx/this->L);                   //calculate minimum image distance
	   rry = rry - this->L*round(rry/this->L);
	   dist = (rrx*rrx + rry*rry + rrz*rrz);                     // distance is square of actual distance
	   cutoff = this->sa_distance(ch,ch1,mem_sel,i);
	   if (dist < cutoff) return true;
	   }
       return false;
      }

	  if( ch1 == 'm'){
	   lscl = link1.getlscl(ch1);                                                         // pointers to the linklist lscl   
	   head = link1.gethead(ch1);                                                         // pointers to the header cell linklist
	   cell = link1.getcells(tempx,tempy,tempz, ch);                                      // All 26 neighbours are returned here
	   for (i=0;i<27;i++){                                                                // To check for the overlap: vesicle-vesicle or antigen-antigen
	    mem = head[cell[i]];                                                            // get the number of nano carrier in cell i
	    if (mem == -1) continue;                                                        // Check for occupied cells
	    do {
	     rrx=tempx - me.getxc(mem);
	     rry=tempy - me.getyc(mem);
	     rrz=tempz - me.getzc(mem);
	     rrx = rrx - this->L*round(rrx/this->L);                                    //calculate minimum image distance
	     rry = rry - this->L*round(rry/this->L);
	     dist = (rrx*rrx + rry*rry + rrz*rrz);                                      // distance is square of actual distance
         cutoff = this->sa_distance(ch,ch1,mem_sel,mem);
	     if (dist < cutoff) return true;                                            // two objects cannot come closer than their diamter
	     mem = lscl[mem];
	    } while (mem !=-1);
	   }  
	   return false;
      }
   }

   else if (ch == 'c'){
    lscl = link1.getlscl(ch1);                                                                         // pointers to the linklist lscl   
    head = link1.gethead(ch1);                                                                         // pointers to the header cell linklist
    tempx = re.getxt(mem_sel);
    tempy = re.getyt(mem_sel);
    tempz = re.getzt(mem_sel);	
    cell = link1.getcells(re.getxc(mem_sel), re.getyc(mem_sel), re.getzc(mem_sel), ch);        // All 26 neighbours are returned here

    if( ch1 == 'v' ){
	for (i=0;i<27;i++){                                                                 //To check for the overlap: vesicle-vesicle or antigen-antigen
	  mem = head[cell[i]];                                                               // get the number of nano carrier in cell i
	  if (mem == -1) continue;                                                           // Check for occupied cells
	  do {
	   rrx = tempx - nc.getxc(mem);
	   rry = tempy - nc.getyc(mem);
	   rrz = tempz - nc.getzc(mem);
	   rrx = rrx - this->L*round(rrx/this->L);                   //calculate minimum image distance
	   rry = rry - this->L*round(rry/this->L);
	   dist = (rrx*rrx + rry*rry + rrz*rrz);                    // distance is square of actual distance
	   cutoff = this->sa_distance(ch,ch1,mem_sel,mem);
	   if (dist < cutoff) return true;                          // two objects cannot come closer than their diamter
	   mem = lscl[mem];
	  } while (mem !=-1);
	}
       return false;
    }

 }
 return false;                                                  // this return statement is to suppress the warning: control-reaches end of non-void function.
}

bool _selfavoidance ::  overlap (int mem_sel, char ch, double dx, double dy, double dz){         //Checking overlap of one vesicle with all other vesicles
    extern _setup a;
    extern _linklist link1;
    int mem,i,*lscl, *head, *cell;
    double cutoff, rrx, rry, rrz, dist;

    lscl = link1.getlscl(ch);                                                       // pointers to the linklist lscl   
    head = link1.gethead(ch);                                                       // pointers to the header cell linklist
	
    tempx = a.getxc(mem_sel, ch)+dx;                                                // new displaced position
    tempy = a.getyc(mem_sel, ch)+dy;
    tempz = a.getzc(mem_sel, ch)+dz;
	
	if (ch=='v'){
	 for(mem=0; mem<this->num_nc; mem++){
	  if (mem != mem_sel){
		rrx = tempx - a.getxc(mem, ch);
		rry = tempy - a.getyc(mem, ch);
		rrz = tempz - a.getzc(mem, ch);
		rrx = rrx - this->L*round(rrx/this->L);             //calculate minimum image distance
		rry = rry - this->L*round(rry/this->L);
		dist = (rrx*rrx + rry*rry + rrz*rrz);              // distance is square of actual distance
		cutoff = this->sa_distance(ch,ch,mem_sel,mem);
		if (dist < cutoff) return true;
		}
	 }
	}

	else {
     cell = link1.getcells(tempx,tempy,tempz, ch);                                   // All 26 neighbours are returned here
     for (i=0;i<27;i++){                                                             // To check for the overlap: vesicle-vesicle or antigen-antigen
	  mem = head[cell[i]];                                                        // get the number of nano carrier in cell i
	  if (mem == -1) continue;                                                    // Check for occupied cells
	  do {
	   if (mem == mem_sel) {mem = lscl[mem]; continue;}                         // skip calculation with itself       	   
	   rrx = tempx - a.getxc(mem, ch);
	   rrx = rrx - this->L*round(rrx/this->L);                    //calculate minimum image distance
	   rry = tempy - a.getyc(mem, ch);
	   rry = rry - this->L*round(rry/this->L);
	   rrz = tempz - a.getzc(mem, ch);
	   dist = (rrx*rrx + rry*rry + rrz*rrz);                     // distance is square of actual distance
	   cutoff = this->sa_distance(ch,ch,mem_sel,mem);
	   if (dist < cutoff) return true;                           // two objects cannot come closer than their diamter
	   mem = lscl[mem];
	  }while (mem !=-1);
     }
	}
    return false;
}

//@ checking overlap of one given object displaced by dx,dy,dz with another class of objects
bool _selfavoidance ::  overlap (int mem_sel, char ch, char ch1, double dx, double dy, double dz){     
    extern _nanocarrier nc;
    extern _membrane me;
    extern _receptor re;
    extern _linklist link1;
    extern _datainput data;
    int mem,i;
    double cutoff, rrx, rry, rrz, dist;
    int *lscl, *head, *cell;

    if (ch == ch1){
      return this->overlap(mem_sel,ch,dx,dy,dz);                                      // same as the above function
    }

   else if((ch == 'v')){ 
    tempx = nc.getxc(mem_sel)+dx;
    tempy = nc.getyc(mem_sel)+dy;
    tempz = nc.getzc(mem_sel)+dz;
    
    if( ch1=='c' ){
     for (i=0 ; i< num_antigens; i++){
      rrx = tempx - re.getxt(i);
      rry = tempy - re.getyt(i);
      rrz = tempz - re.getzt(i);
      rrx = rrx - this->L*round(rrx/this->L);                   //calculate minimum image distance
      rry = rry - this->L*round(rry/this->L);
      dist = (rrx*rrx + rry*rry + rrz*rrz);                     // distance is square of actual distance
      cutoff = this->sa_distance(ch,ch1,mem_sel,i);
		//cout<<str_green<<rrx<<" "<<rry<<" "<<rrz<<" "<<dist<<" "<<cutoff<<endl;
      if (dist < cutoff) return true;
     }
     return false;
    }

    if( ch1 == 'm'){
	lscl = link1.getlscl(ch1);                                                  // pointers to the linklist lscl   
	head = link1.gethead(ch1);                                                  // pointers to the header cell linklist
	cell = link1.getcells(tempx,tempy,tempz,ch);                               // All 26 neighbours are returned here
	for (i=0; i<27; i++){                                                       //To check for the overlap: vesicle-vesicle or antigen-antigen
	   mem = head[cell[i]];                                                     // get the number of nano carrier in cell i
	   if (mem == -1) continue;                                                 // Check for occupied cells
	   do {
	    rrx = tempx - me.getxc(mem);
	    rry = tempy - me.getyc(mem);
	    rrz = tempz - me.getzc(mem);
	    rrx = rrx - this->L*round(rrx/this->L);                   //calculate minimum image distance
	    rry = rry - this->L*round(rry/this->L);
	    dist = (rrx*rrx + rry*rry + rrz*rrz);                     // distance is square of actual distance
	    cutoff = this->sa_distance(ch,ch1,mem_sel,mem);
	    if (dist < cutoff) return true;                           // two objects cannot come closer than their diamter
	    mem = lscl[mem];
	   } while (mem !=-1);
	}
	return false;
	}
   }

    else if (ch=='c'){
        lscl = link1.getlscl(ch1);                                                               // pointers to the linklist lscl   
        head = link1.gethead(ch1);                                                               // pointers to the header cell linklist
        tempx = re.getxt(mem_sel)+dx;
        tempy = re.getyt(mem_sel)+dy;
        tempz = re.getzt(mem_sel)+dz;

        cell = link1.getcells(re.getxc(mem_sel)+dx, re.getyc(mem_sel)+dy, re.getzc(mem_sel)+dz, ch);  // All 26 neighbours are returned here

        if( ch1 == 'v' ){
         for (i=0; i<27; i++){                                                       //To check for the overlap: vesicle-vesicle or antigen-antigen
          mem = head[cell[i]];                                                       // get the number of nano carrier in cell i
          if (mem == -1) continue;                                                   // Check for occupied cells
          do {
           rrx = tempx - nc.getxc(mem);
           rry = tempy - nc.getyc(mem);
           rrz = tempz - nc.getzc(mem);
           rrx = rrx - this->L*round(rrx/this->L);                //calculate minimum image distance
           rry = rry - this->L*round(rry/this->L);
           dist = (rrx*rrx + rry*rry + rrz*rrz);                  // distance is square of actual distance
	   cutoff = this->sa_distance(ch,ch1,mem_sel,mem);
           if (dist < cutoff) return true;                        // two objects cannot come closer than their diamter
            mem = lscl[mem];
          } while (mem !=-1);
        }
        return false;
       }

    }
    return false;                                           // this return statement is to suppress the warning: control-reaches end of non-void function.

}

bool _selfavoidance ::  overlap_all (int mem_sel, char ch, char ch1){           //Checking overlap of one vesicle with all other vesicles
    extern _setup a;
    int mem,num_members;
    double cutoff, rrx, rry, rrz, dist;      
    tempx = a.getxc(mem_sel, ch);
    tempy = a.getyc(mem_sel, ch);
    tempz = a.getzc(mem_sel, ch);
    num_members = a.getnum_members(ch1);

    for (mem=0; mem<num_members; mem++){                         // To check for the overlap: vesicle-vesicle or antigen-antigen
      if ((ch == ch1) && (mem == mem_sel))  continue;              // skip calculation with itself       
       rrx = tempx - a.getxc(mem, ch1);
       rry = tempy - a.getyc(mem, ch1);
       rrz = tempz - a.getzc(mem, ch1);
       rrx = rrx - this->L*round(rrx/this->L);                   // calculate minimum image distance
       rry = rry - this->L*round(rry/this->L);
       dist = (rrx*rrx + rry*rry + rrz*rrz);                    // distance is square of actual distance
       cutoff = this->sa_distance(ch,ch1,mem_sel,mem);
       if (dist < cutoff) return true;                           // two objects cannot come closer than their diamter
    }
    return false;
}

bool _selfavoidance ::  overlap_all (char ch, char ch1){                      //Checking overlap of one vesicle with all other vesicles
    extern _setup a;
    extern _datainput data;
    int m,num_members1;
    bool flag = false;
    
    num_members1 = a.getnum_members(ch);
    for (m=0 ; m<num_members1; m++){
      flag = this->overlap_all(m,ch, ch1);
      if (flag) break;
    }
     return false;
}

bool _selfavoidance ::  overlap_antigens (int mem1, int mem2){         
    extern _receptor re;
    double rrx, rry, rrz, dist;
    double tempx1,tempy1,tempz1;
    bool res = false;
    
    tempx = re.getxc(mem1);
    tempy = re.getyc(mem1);
    tempz = re.getzc(mem1);
    tempx1 = re.getxc(mem2);
    tempy1 = re.getyc(mem2);
    tempz1 = re.getzc(mem2);

    rrx = tempx - tempx1;
    rry = tempy - tempy1;
    rrz = tempz - tempz1;
    rrx = rrx - this->L*round(rrx/this->L);                    //calculate minimum image distance
    rry = rry - this->L*round(rry/this->L);
    dist = (rrx*rrx + rry*rry + rrz*rrz);                      // distance is square of actual distance
    if (dist < this->sa_distance('c','c',mem1,mem2)){
      res=true;
    }
    return res;
}
