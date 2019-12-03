// Detailed notes can be found in ``readme.md`` found one level higher in the directory tree
#include <mpi.h>                                                                             // mpi.h should come first first for compilation with openmpi
#include "_declarations.h"
#include "_datainput.h"                          // definition for data reading from input.in
#include "_nanocarrier.h"                       // Lays down the datastructure for system variables
#include "_membrane.h"
#include "_receptor.h"
#include "_setup.h"
#include "_linklist.h"               // Construct a linklist connecting space to nano carrier and memb.
#include "_potential.h"              // Computes various potential energies	
#include "_init.h"                   // Initialize the system with parameters read 
#include "_distribution.h"           // Compute various distributions
#include "_datawriter.h"             // Write down the data in known formats for visualization
#include "_histogram.h"
#include "_movement.h"              // moves that allow the system to explore various states	
#include "_restart.h"
#include "_fortran_modules.h"
#include "_fortran_structures.h"
#include "_selfavoidance.h"
#include "_debugger.h"
#include "_interaction_parameters.h"

MTRand mtrand1;
_datainput data;                                 // data is an object of type _datainput
_nanocarrier nc;                                 // nanocarrier
_receptor re;                                    // receptor
_membrane me;                                    // membrane
_setup a;
_linklist link1;
_potential p;
_init init;                                     // Calls to initialize the Nanocarrier, Antibodies and Antigens	
_datawriter w;
_histogram hist;
_movement mv;
_restart RES;
_distribution distrib;
_selfavoidance sa;
_debugger debug;
_intparam intparam;


bool debug_mode, frequent_dump, continue_run;         // if set to true, the code prints to screen and dumps files at regular interval  
int processor, nprocessors, window;

int main(int argc, char* argv[]) {

  //@ MPI intialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocessors);
  MPI_Comm_rank(MPI_COMM_WORLD, &processor);
  __module_datastruct_MOD_processor_number = processor;
  window=processor;
  
  //@ random number initialization (one for each MPI core)
  mtrand1.seed();
  
  
  //@ To restart or not
  continue_run = false;
  if ((argc > 1) && (strcmp(argv[1],"true")==0)) {
    continue_run = true;
  }
  
    debug_mode = false;
    __module_datastruct_MOD_debug_mode = false;
	
   //@ To print debug info or not
  if ((argc > 2) && (strcmp(argv[2],"debug")==0)) {
    debug_mode = true;
    __module_datastruct_MOD_debug_mode = true;
  }

  //@ To dump configurations frequently
  frequent_dump = false;  
  if ((argc > 3) && (strcmp(argv[3],"dump")==0))  
	  frequent_dump = true;
  
  
  if (debug_mode){
	  cout<<str_red<<"\t<MPI Info> "<<str_black<<endl;
	  cout<<"\t --> MPI Initialization begins"<<endl;
	  cout<<"\t --> Initialization of MPI over"<<endl;
	  cout<<"\t --> Job is running on "<<nprocessors<<" processors"<<endl;
	  cout<<"\t --> Processor is "<<processor<<endl;
	  cout<<"\t\t\t ============================ \n"<<endl;
	}
  
  //@ Initialize the restart module
  RES.init();
  
  //@ System Initialization
  if (!continue_run){           // start fresh
	cout<<str_green <<"\t Running a newly initialized system" <<endl<<str_black;
  	data.init();                  // Read the input parameters and initialize class _datainput
  	nc.init();                    // Initialize the Nanocarrier and Antigens
  	re.init();
    cout<<"Coming here 1"<<endl;
  	intparam.init();
    cout<<"Coming here 2"<<endl;
    me.init();
  }
  else{                                  //restart
	cout<<str_blue <<"\t \t <--> Restarting from a previously stored state <-->" <<endl<<str_black;
	RES.read_restart_conf(processor);
  }
  
  runparam RP = data.get_runparameters();
  
  
  //@ set the linklist structure for each component
  link1.init();                                            // Initialize the linklist 
  link1.calc('v');
  link1.calc('c');
  link1.calc('m');
  
  //@ interaction parameters, potentials and self avoidance
  w.init();
  hist.init();                                           // Histogram class is initialized
  mv.init();                                             // movement class is intialized 
  distrib.init();
  p.init();
  sa.init();
  debug.init();

  stringstream s1,s2,s3,s4,s5,s6,s7,s8,s9,s10;
  ofstream nc_angle,nc_ab_data,antigen_data,multivalency_file,zhist_file,Hhist_file,ti_ener; 

  int  antigen, datastep;
  struct nanocarrier *ncpointer;
  ncpointer =&__module_datastruct_MOD_nc_f;

  int type_of_move;
  double type_of_trans, ratio_m;

  link1.calc('m');
  link1.calc('c');                                                                   // Link list for antigen (antigenlink.calc)

  //@ Membrane geometry parameters
  int nver = me.getnum_vertices();
  int nlink = me.getnum_links();
  double flip_move_ratio = (double)(nver)/(double)(nlink);
  
  
  //@ Equilibrate the membrane if the starting configuration is a flat membrane
  if (! continue_run){
  double halfDTMC = (int) RP.DTMCsteps/4;
  int twindow = 0;
  int ccstep = 0;
  
  cout<<RP.system<<endl;
  
  if (__module_datastruct_MOD_memb_restart_flag == false && RP.DTMCflag){      
   if (debug_mode) cout << str_blue<<"Equilibrating membrane "<<endl;
   
    for (int dtmcs = 0 ; dtmcs < RP.DTMCsteps; dtmcs++){
      cout <<"Coming here in step1";
      if (((dtmcs%1000)==0) && debug_mode) cout<<"Completed membrane equilibration step "<<dtmcs<<endl;
      for (int iloop=0; iloop<nlink; iloop++){
       type_of_trans = mtrand1.rand();
       if (type_of_trans > flip_move_ratio) 
        me.movevertex(true);                                                   
       else
        me.fliplink(true);
      }
      if((dtmcs >= halfDTMC) && (dtmcs%RP.conf_freq == 0)){
         ccstep += 1;
         __module_writedata_MOD_write_membrane_xyz(&RP.ensno,&twindow,&ccstep);
         __module_writedata_MOD_write_membrane_area(&RP.ensno,&twindow,&ccstep);       // print membrane area at this interval
         if (ccstep%100 == 0){
            int dstepp = 1000+(int) ccstep/100;
            me.dump_membrane_conf(RP.ensno,dstepp);
            w.VTK_membrane(RP.ensno,twindow,dstepp);
         }
       }
    }
       
    if (debug_mode) {
        cout << str_blue<<" \t\t --> Membrane Equilibration complete"<<endl;
        cout<<endmarker<<endl;
     }
     
    if (RP.system.compare("membrane_DTMC") == 0) 
        exit(1);
    }
 
   //@ reinitialize the link list after membrane equilibration
   link1.calc('m');
   link1.calc('c');

	//@ Put the NC in the required position and update the linklist
	init.init(2);
	link1.calc('v');                                             // compute the linklist for vesicle (veslink.calc)
	mv.reset_bond_and_equilibrate(20000, RP.DTMCflag);
	}

  w.VTK_membrane(RP.ensno,0);
  w.VTK_Antigens(RP.ensno,0);
  w.VTK_NanoCarrier(RP.ensno,0);

  double total_energy = p.total_energy();

  if (debug_mode){
    cout<<str_black<<" \t --> Total energy after equilibriation "<<total_energy<<endl;
    cout<<"\t --> Window conf interval: "<<RP.conf_freq<<endl;
    cout<<"\t --> Window antab interval: "<<RP.antab_freq<<endl;
    cout<<"\t --> On screen print interval: "<<RP.osdp_freq<<endl;
    cout<<"\t --> Antigen membrane configurations interval: "<<RP.ant_write_freq<<endl;
  }

  if (debug_mode){
    cout<<str_red<<"\n \t  <-------- Running in debug mode : \
    Results are printed to screen (recommended only for test runs) --------------> \n"<<str_black<<endl;
  }
  else{
    cout<<str_red<<"\n \t <---------- Running in silent mode : \
    Only essential files are written (recommended for production runs) ------------>\n "<<str_black<<endl;
  }

  // @ store the number of nanocarriers and number of antibodies in each NC
  int num_nanocarrier = nc.getnum_members();
  int num_antigens = re.getnum_members();
  int num_ab[num_nanocarrier];
  for (int i=0; i<num_nanocarrier; i++)
    num_ab[i] = nc.getnum_ab(i);
  
  int curr_multivalency;                                                                       // Store multivalency as a array
  ratio_m = (num_nanocarrier*num_ab[0]*5.0)/(num_nanocarrier*num_ab[0]*5.0+num_antigens);      // ratio of antibody movement to antigen movement
  
  if (debug_mode){
    for(int i=0;i<num_nanocarrier;i++){
    cout<<"\t\t NC# "<<i<<",  Mode = "<<nc.get_bias_mode(i);
    cout<<", k = "<<nc.get_kbias(i);
    cout<<", r_0 = "<<nc.get_biasref(i);
    cout<<", dr = "<<nc.get_biasbinsize(i);
    cout<<", approach = "<<nc.get_biasdir(i);
    cout<<", lambda = "<<nc.get_lambda(i)<<endl;
    }
  }
  
  /*set up the arrays to store the data */
  int hist_array_size = (int)((RP.step_end-RP.sampling_cutoff)/(RP.samp_freq));                              // frequency at which the WHAM info is sampled
  int syst_array_size = (int)((RP.step_end-RP.sampling_cutoff)/(RP.antab_freq));                             // frequency at which the various geometric quatities are sampled
  
  for (window = RP.win_start; window < RP.win_end; window++){
	  
	cout<<str_green<<endl<<"\t --> Running window: "<<window<<"/"<<RP.win_end<<endl<<endl<<endl<<str_black;
	cout<<str_blue<<" \t Bias values: ";
	for (int i=0;i<num_nanocarrier;i++)
		cout<<" NC#"<<i<<" : "<<nc.get_biasref(i)<<", ";
	cout<<endl;
    
    //@ Reinitialization of arrays to store the various measures required in our calculations
    int hist_array_counter = 0, syst_array_counter=0,conf_no=0;
    double **histogram_array,**memb_meanr, **localcurv_array = NULL;                 // set to NULL because this array is used only when the bias_mode is 'H'
    double **multivalency_array = NULL;                                                                
    double **bond_energy,*flexure_energy,*membrane_energy;
	double **orientationx=NULL,**orientationy=NULL,**orientationz=NULL;
	vector < string > ncabdata,antdata;  
    
    //@ nanocarrier separation, curvature below, and multivalency
    histogram_array = new(nothrow) double*[hist_array_size];                          // Store all the histogram data in the array and finally write it to a file
    memb_meanr = new(nothrow) double*[hist_array_size];                                                                      
    localcurv_array = new (nothrow) double*[hist_array_size];                      // array to store the curvature of the membrane in the shadow of the nanocarrier  
    multivalency_array = new (nothrow) double*[hist_array_size];                   //array to store the number of multivalent interactions

    for (int i=0;i<hist_array_size;i++){
      histogram_array[i]=new(nothrow) double[num_nanocarrier];
      memb_meanr[i]=new(nothrow) double[num_nanocarrier];
      localcurv_array[i] = new(nothrow) double[num_nanocarrier];  
	  multivalency_array[i] = new(nothrow) double[num_nanocarrier];
    }
  
    //@ system measures
    orientationx = new (nothrow) double*[syst_array_size];                          // array to store the curvature of the membrane in the shadow of the nanocarrier
    orientationy = new (nothrow) double*[syst_array_size];                          // array to store the curvature of the membrane in the shadow of the nanocarrier
    orientationz = new (nothrow) double*[syst_array_size];                          // array to store the curvature of the membrane in the shadow of the nanocarrier
    bond_energy = new(nothrow) double*[syst_array_size];                           // array to store the energy due to formation of bonds
    for (int i=0;i<syst_array_size;i++){
      orientationx[i] = new(nothrow) double[num_nanocarrier];  
      orientationy[i] = new(nothrow) double[num_nanocarrier];  
      orientationz[i] = new(nothrow) double[num_nanocarrier];  
      bond_energy[i] = new(nothrow) double[num_nanocarrier];
    }
    flexure_energy = new(nothrow) double[syst_array_size];                         // array to store the energy cost due to flexure of antigens
    membrane_energy = new(nothrow) double[syst_array_size];                        // array to store the bending energy of the membrane
    //@ end of initializations 


   //@ Monte Carlo dynamics 
    for (int curr_step = RP.step_start; curr_step < RP.step_end; curr_step++){        // Production run ( a pre-equilibrated membrane has been used here)
      __module_datastruct_MOD_curr_step= curr_step;                                    // if in case you like to access current step number from fortran

      if ((curr_step%RP.osdp_freq == 0) && debug_mode){                                // print data to file and screen at intervals of sampling_interval
        cout<<str_blue<<"\n \t Window --> "<<window+1<<"/"<<RP.win_end<<"\t: Step --> "<<curr_step<<"/"<<RP.step_end<<endl;
        cout<<"\t Total Energy = "<<p.total_energy()<<endl<<endl;
        for(int i=0; i<num_nanocarrier; i++){
         double* stepsize = nc.get_mc_stepsizes(i);
         int* bond_counter = nc.get_bond_counters(i);
         cout<<str_green<<"\t NC# "<<i<<", <m> = "<<bond_counter[0]<<str_black;
         cout<<str_blue<<", N_sha = "<<(ncpointer+i)->nshadow_vert<<
	      ", <R> = "<<setprecision(5)<<nc.get_shadowr_dist(i)<<
	      ", <H> = "<<setprecision(5)<<(ncpointer+i)->meancurv<<endl;
         cout<<"\t Step Size Statistics: "<<stepsize[0]<<" "<<stepsize[1]<<" "<<stepsize[2]<<" "<<stepsize[3]<<str_black<<endl;
         cout<<str_red<<"\t Number of bond formed and broken : "<<bond_counter[1]<<" "<<bond_counter[2]<<endl<<str_black;
         cout<<"\t|"<<endl;
         cout<<"\tV"<<endl;
         nc.reset_bond_counters(i);
        }
        cout<<"\t===================================>"<<endl;
      }

      type_of_move = mtrand1.rand(RP.nmcss, RP.nmcse);                                        // random number generator itself returns 0 or 36
      if (type_of_move <= 6.0){
        type_of_trans = mtrand1.rand(); 
        if (type_of_trans <= ratio_m/2) nc.translation();                                     // translate the nanocarrier
        if (type_of_trans > ratio_m/2 && type_of_trans <= ratio_m) nc.rotation();             // rotate the nanocarrier
        if (type_of_trans > ratio_m) re.receptor_hopping();                                   // receptors hopping from one vertex to another
      }

      else if ((type_of_move > 6.0) && (type_of_move < 24.0)){
        type_of_trans = mtrand1.rand();
        if (type_of_trans <= 0.6) {
          nc.bond_formation();
		}
        else{
          re.receptor_diffusion();
		}
      }

      else if ((type_of_move >= 24.0) && (type_of_move <=36.0)){                     // thermalize individual vertex and links in the membrane
        type_of_trans = mtrand1.rand();
        if (type_of_trans > flip_move_ratio )
          me.fliplink();
        else
          me.movevertex();
      }
  
      if (((curr_step%5000) == 0) && RP.DTMCflag) me.membrane_montecarlo();          //@ Additional thermalization of the membrane

  
      //@ Compute various measures in the system
      if (curr_step >= RP.sampling_cutoff){
        
        if (curr_step%RP.samp_freq == 0){
          for(int i=0;i<num_nanocarrier;i++){
            histogram_array[hist_array_counter][i] = nc.get_shadowr_dist(i);
            localcurv_array[hist_array_counter][i]= nc.get_shadowcurv(i);
			multivalency_array[hist_array_counter][i] = nc.get_multivalency(i);
        }
        hist_array_counter++;
        }
  
       //@ write the datafiles for antigen and membrane configurations
       //@ control the frequency of this module by setting the value of antigen_memb_write_interval
       
       if((curr_step%RP.ant_write_freq==0)){
		cout<<"writing this file at "<<curr_step<<" "<<RP.ant_write_freq<<endl;
         datastep = (int) curr_step/RP.ant_write_freq;
         w.write_antigen_xyz(processor,window,datastep);
         __module_writedata_MOD_write_membrane_xyz(&RP.ensno,&window,&datastep);
         __module_writedata_MOD_write_membrane_area(&RP.ensno,&window,&datastep);       // print membrane area at this interval
         distrib.compute_rdf_distribution();                                          // compute gofr and also dump the datafile for further post processing
         if (datastep%10 == 0)
           distrib.write_rdf_distribution(processor,window);
       }
       
       if (curr_step%RP.ant_traj_freq==0){
		re.store_trajectory(RP.ensno);   
		}
  
       bool bonded_flag = false;
	   
       if((curr_step%RP.antab_freq==0)){
		flexure_energy[syst_array_counter] = p.total_Benergy();
		membrane_energy[syst_array_counter] = __module_datastruct_MOD_total_membrane_energy();  
		stringstream ncabdatat;
  	  	stringstream antdatat;
		 
		for (int i=0;i<num_nanocarrier;i++){                                                             //output the antigen bond length and angle
			curr_multivalency = nc.get_multivalency(i);
			orientationx[syst_array_counter][i] = nc.getanglex(i);
			orientationy[syst_array_counter][i] = nc.getangley(i);
			orientationz[syst_array_counter][i] = nc.getanglez(i);
			bond_energy[syst_array_counter][i] = p.calc(i,'v');
			
			if(curr_multivalency != 0) {
				ncabdatat<<curr_multivalency<<" "<<nc.getxc(i)<<" "<<nc.getyc(i)<<" "<<nc.getzc(i)<<" ";
				for (int j=0; j<num_ab[i]; j++) {
					if (nc.getbonded(i,j)){
						ncabdatat<<nc.getx_antibody(i,j)<<" "<<nc.gety_antibody(i,j)<<" "<<nc.getz_antibody(i,j)<<" ";
						antigen = nc.getbondedto(i,j);
						ncabdatat<<re.getxt(antigen)<<" "<<re.getyt(antigen)<<" "<<re.getzt(antigen)<<" ";
 	            		antigen = nc.getbondedto(i,j);
 	            		antdatat<<antigen<<" "<<nc.compute_bond_length(i,j,antigen)<<" "<<re.get_theta(antigen)<<" "<<re.get_phi(antigen)<<" ";
						bonded_flag = true;	      	 						
					}
				}
			}
		 }
		 
		 if (bonded_flag){
 		 ncabdata.push_back(ncabdatat.str());
		 antdata.push_back(antdatat.str());
	     }
		 syst_array_counter++;
		 ncabdatat.str("");
		 antdatat.str("");
	    }
     }
      
      //@ dump the configuration files
      if ( frequent_dump && (curr_step%RP.conf_freq == 0)){
        conf_no +=1;
        w.VTK_membrane(RP.ensno,window,conf_no);
		w.VTK_NanoCarrier(RP.ensno,window,conf_no);
		w.VTK_Antigens(RP.ensno,window,conf_no);
	  }
      else if((curr_step>RP.sampling_cutoff) && (curr_step%RP.conf_freq==0)){
        conf_no +=1;
		w.VTK_membrane(RP.ensno,window,conf_no);
		w.VTK_NanoCarrier(RP.ensno,window,conf_no);
		w.VTK_Antigens(RP.ensno,window,conf_no);
     }
    }                                                                                         //@end of Monte Carlo Sampling for each window

	RES.write_restart_conf(RP.ensno,window);
    me.dump_membrane_conf(RP.ensno);
	re.dump_trajectory(RP.ensno);

    // Write statistical measures such as multivalency, only when equilibrium calculations/WHAM calculations are performed 
    if ( RP.system.compare("DTMC_membrane") != 0) {
		//@ Write the data for the NC angle
		s4<<"NC-angle-ENS-"<<RP.ensno<<"_Frame-"<<window<<".dat";
		nc_angle.open(s4.str().c_str(), ios::out);
		for (int j=0; j<num_nanocarrier; j++)
			nc_angle<<"NC#"<<j<<"(x) "<<"NC#"<<j<<"(y) "<<"NC#"<<j<<"(z) "<<" ";
		nc_angle<<endl;
			
		for (int i=0;i<syst_array_size;i++){
			for (int j=0; j<num_nanocarrier; j++)
				nc_angle<<orientationx[i][j]<<" "<<orientationy[i][j]<<" "<<orientationz[i][j]<<" ";
			nc_angle<<endl;
		}
		nc_angle.close(); s4.str("");                                                        // close the file
		//--------------------------------->
		
		//@ Write bound antibody-antigen data
		s5<<"NC-Antibody-Data-ENS-"<<RP.ensno<<"_Frame-"<<window<<".dat";                  // x,y, phi,theta,psi for the nanocarrier
		nc_ab_data.open(s5.str().c_str(), ios::out | ios::binary); 
		for (size_t ii=0; ii<ncabdata.size(); ii++)
				nc_ab_data<<ncabdata[ii]<<endl;			
		nc_ab_data.close(); s5.str("");
		//--------------------------------->
		
		//@ Write Antigen flexure data 
		s6<<"Antigen-Data-ENS-"<<RP.ensno<<"_Frame-"<<window<<".dat";                    // bonddist,theta,phi for the antigen
		antigen_data.open(s6.str().c_str(), ios::out);
		for (size_t ii=0; ii<antdata.size(); ii++)
				antigen_data<<antdata[ii]<<endl;
		antigen_data.close(); s6.str("");		
		//--------------------------------->

		
		//@ Write Multivalency data
		s7<<"Multivalency-ENS-"<<RP.ensno<<"_Frame-"<<window<<".dat";         // bonddist,theta,phi for the antigen
		multivalency_file.open(s7.str().c_str(), ios::out);
		//@ header for multivalency file
		for (int i=0; i<num_nanocarrier; i++)
			multivalency_file<<"NC#"<<i<<" ";
		multivalency_file<<endl;
		//@ multivalency data
		for(int i1=0; i1<hist_array_size;i1++){
			for(int i=0;i<num_nanocarrier;i++){
				multivalency_file<<multivalency_array[i1][i]<<" ";
			} 
			multivalency_file<<endl;
		}
		multivalency_file.close(); s7.str(" ");
		//--------------------------------->
		
		//@ Write the energy data and also use it for TI calculation 
		s10<<"Energy-ENS-"<<RP.ensno<<"-Frame-"<<window<<".dat";
		ti_ener.open(s10.str().c_str(), ios::out);
		for (int i=0; i<num_nanocarrier; i++){
			ti_ener<<"#NC: "<<i<<"  #mode: "<<nc.get_bias_mode(i)<<"  #lambda: "<<nc.get_lambda(i)<<endl;
		}
		for(int i=0; i<syst_array_size; i++){
			for(int j=0; j<num_nanocarrier; j++)
				ti_ener<<bond_energy[i][j]<<"  ";                            //@ 2 to 2+num_nc column are the bonded energy which will be used in TI calculations				
			ti_ener<<flexure_energy[i]<<"  "<<membrane_energy[i]<<endl;      //@ last two column are the flexure and membrane energy
		}
		ti_ener.close();s10.str("");	
		//--------------------------------->
		
		
		s8<<"Histogram_ENS-"<<RP.ensno<<"_Frame-"<<window<<".rawdata";                  // write the histogram data
      	zhist_file.open(s8.str().c_str(),ios::out|ios::binary);
	  	s9<<"Curvature_Histogram_ENS-"<<RP.ensno<<"_Frame-"<<window<<".rawdata";        // write the histogram data
      	Hhist_file.open(s9.str().c_str(),ios::out|ios::binary);
      
      	//@ header for R histogram data
      	zhist_file<<"#Frame: "<<window<<endl;
      	for (int i=0; i<num_nanocarrier; i++)
       		zhist_file<<"mode: "<<nc.get_bias_mode(i)<<"  kbias: "<<nc.get_kbias(i)<<"   rref: "<<nc.get_biasref(i)<<endl;
      
		zhist_file<<"# freq: "<<RP.samp_freq<<" ";
      	for (int i=0; i<num_nanocarrier; i++)
       		zhist_file<<", NC#"<<i<<" ";
      	zhist_file<<endl;
      	zhist_file<<"####### end of header #############"<<endl;
      
      	//@ header for H histogram data
      	Hhist_file<<"#Frame: "<<window<<endl;
      	for(int i=0;i<num_nanocarrier;i++)
        	Hhist_file<<"mode: "<<nc.get_bias_mode(i)<<"  kbias: "<<nc.get_kbias(i)<<"   rref: "<<nc.get_biasref(i)<<endl;
      	Hhist_file<<"# freq: "<<RP.samp_freq<<" ";
      	for (int i=0; i<num_nanocarrier; i++)
       		Hhist_file<<", NC#"<<i<<" ";
      	Hhist_file<<endl;
      	Hhist_file<<"####### end of header #############"<<endl;
      
      	for(int i1=0; i1<hist_array_size;i1++){
		  zhist_file<<setprecision(10);
		  Hhist_file<<setprecision(10);
       	for(int i=0;i<num_nanocarrier;i++){
         zhist_file<<histogram_array[i1][i]<<" ";
         Hhist_file<<localcurv_array[i1][i]<<" ";
       	} 
      	zhist_file<<endl;
      	Hhist_file<<endl;
      	}
	    zhist_file.close(); s8.str(" ");
    	Hhist_file.close(); s9.str(" ");
    }	

    delete [] histogram_array;                                                                                 // delete histogram_array
    delete [] memb_meanr;                                                                                      // delete histogram_array
    delete [] localcurv_array;                                                                                 // delete the local curvature array
    delete [] multivalency_array;
    delete [] bond_energy;
    delete [] flexure_energy;
    delete [] membrane_energy;
	delete [] orientationx;
	delete [] orientationy;
	delete [] orientationz;
	//antdata.clear();
	//ncabdata.clear();
	

    nc.shift_biaswindow();                                                     // shift the biaswindow to the next position if a PMF calculation is requested
  }  //@ looping over windows end here
 
  cout<<"\nProgram completed on processor "<<processor<<endl;
  cout<<"Waiting for the rest to complete \n"<<endl;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(1);
}
