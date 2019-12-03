def cpp_parameters(parameter):
	fp=open('./CC_INPUT_FILES/input.in','w')
	for f in parameter:		
		fp.write(str(f[0])+'\t'+str(f[1])+'\n')
	fp.close()

def fort_parameters(parameter):
	fp=open('./FC_INPUT_FILES/parameters.in','w')
	for f in parameter:
		fp.write(str(f[0])+'\t'+str(f[1])+'\n')
	fp.close()
	
#!/bin/env python
import numpy,sys,os,time,shutil,stat,subprocess,glob



# fortran variables
antigen_density  = 0.002                                                      # 0.002 -> 2000 antigens per micron^2
kappa,pressure=20,0.0                                                         # kappa, pressure, and the excess surface area paramter
blen_fac=[1.25,1.5]

czero,czero_pattern,czero_conc = 0.0, 'CIRCULAR', 0.0
if(czero_conc == 0.0):
	czero=0.0
	
annulus_start,annulus_end = 3,10
antigen_restart_flag = '.False.'
fixed_frame_flag = '.False.'
restart_memb=0                                                               # set this flag to 1 if the membrane is to be restarted




#C++ variables
num_moves=100000000
start_ens,end_ens = 1,5                                                      # number of ensemble
N,init_bond_length= 50,10                                                     # grid size and initial grid spacing     
geom='PLANAR'
periodicity=2
height=20
flexural_rigidity=7000
vermov_step_adjust=10000
bias_mode ='Z'
bias_pot_fluctz,z_start,zbiasstart_window,numz_windows,zwindow_step_size=1.0,84.45,0,50,0.05
bias_pot_fluctH,H_start,Hbiasstart_window,numH_windows,Hwindow_step_size=1.0,0.5,0,20,0.05

window_conf_interval=num_moves/100
window_antab_data_interval=num_moves/10000
osd_print_interval=num_moves/100000
dcd_frequency=100


source_dir='../../../../../../../../..'
codehome=os.getcwd()
exec_name='NC_MEMB_fixed'

if (bias_mode == 'H'):
	dirname1='Fixed_Hbias_N-'+str(N)+"_bl-"+str(init_bond_length)
	dirname1_rs='Fixed_Hbias_N-'+str(N)+"_bl-"+str(init_bond_length)+'-rs'
elif (bias_mode == 'Z'):
	dirname1='Fixed_Zbias_N-'+str(N)+"_bl-"+str(init_bond_length)
	dirname1_rs='Fixed_Zbias_N-'+str(N)+"_bl-"+str(init_bond_length)+'-rs'
else:
	dirname1='Fixed-N'+str(N)+"_bl-"+str(init_bond_length)
	dirname1_rs='Fixed-N'+str(N)+"_bl-"+str(init_bond_length)+'-rs'
	
subprocess.call(["mkdir","-pv",dirname1])
os.chdir(dirname1)

curdir_out=os.getcwd()
for blen_scale_factor in blen_fac:

	dirname10='blen-scale-factor-'+str(blen_scale_factor)
	subprocess.call(["mkdir","-pv",dirname10])
	os.chdir(dirname10)
	
	dirname11='Flexural_Rigidity-'+str(flexural_rigidity)
	subprocess.call(["mkdir","-pv",dirname11])
	os.chdir(dirname11)
	
	if geom=='SINUSOIDAL':
		dirname2='SINUSOIDAL_Period'+str(periodicity)+'_depth'+str(height)
	if geom=='PLANAR':
		dirname2=geom
	subprocess.call(["mkdir","-pv",dirname2])
	os.chdir(dirname2)
	dirname3='zstart-'+str(z_start)
	subprocess.call(["mkdir","-pv",dirname3])
	os.chdir(dirname3)
	
	dirname40='kappa-'+str(kappa)+'-pressure-'+str(pressure)+'-czero-'+str(czero)+'-conc-'+str(czero_conc)
	subprocess.call(["mkdir","-pv",dirname40])
	os.chdir(dirname40)
	num_antigens=int((N*init_bond_length)**2*antigen_density)                            # computed from the density predefined in the system (2000 antigens per micron^2)
	
	
		
	num_antibody=[75,100,162]
	for nab in num_antibody:
		curdir=os.getcwd()
		dirname4='AB-'+str(nab)
		subprocess.call(["mkdir","-pv",dirname4])
		os.chdir(dirname4)
		subdir=os.getcwd()


		for ens in range(start_ens,end_ens):
			run_exe='equm-'+str(nab)+'-N'+str(N)+'E-'+str(ens)
			dirname5='ENS-'+str(ens)
			os.mkdir(dirname5)
			os.chdir(dirname5)
			subdir1=os.getcwd()

			if bias_mode == 'H':
				firstval,binsize,nwindow = H_start,Hwindow_step_size,numH_windows
			elif bias_mode == 'Z':
				firstval,binsize,nwindow = z_start,zwindow_step_size,numz_windows
			else:
				print 'Invalid option for bias_mode : ',bias_mode
				print 'This script is only meant to split each bias window and run them individually'
				print 'Check other scripts in the same direcetory'
				sys.exit('Exiting')

			for nwin in range(0,nwindow):
				dirname51 = 'window-'+str(nwin)
				os.mkdir(dirname51)
				os.chdir(dirname51)

				os.mkdir('RUNDIR')
				memb_restart,res_file_num=".False.",0
				if restart_memb ==1:
					restart_dir=codehome+'/'+dirname1_rs+'/'+dirname10+'/'+dirname11+'/'+dirname2+'/'+dirname3+'/'+dirname40+'/'+dirname4+'/'+dirname5+'/SYS_STATE/DUMP'
					print restart_dir
					nfiles=glob.glob(restart_dir+'/*.dump').__len__()-1
					memb_restart,res_file_num=".True.",nfiles
					memb_res_file=restart_dir+'/membrane_state-'+str(nfiles)+'.dump'
				print memb_restart,res_file_num

				if bias_mode == 'Z':
					zref,zstart,zend = z_start-(nwin*zwindow_step_size),nwin,nwin+1
					Href,Hstart,Hend = 0.0,0,0
				elif bias_mode == 'H':
					Href,Hstart,Hend = H_start-(nwin*Hwindow_step_size),nwin,nwin+1
					zref,zstart,zend = 0.0,0,0
	
	
				cpp_param=[['periodic_box_length_height_ratio',2.0],
				['link_cell_length',50],
				['num_vesicles',1],
				['radius',50.0],
				['antibody_size',15],
				['num_antigens',num_antigens],
				['antigen_radius',1.5],
				['num_antibody',nab],
				['flexural_rigidity',flexural_rigidity],
				['eqm_bond_dist',19.0],
				['eqm_bond_energy',-7.98e-20],
				['spring_constant',1000.0],
				['temperature',300],
				['dcd_freq',dcd_frequency],
				['N',N],
				['memb_init_bond_length',init_bond_length],
				['num_moves',num_moves],
				['init_transstep_size',0.5], 
				['init_anglstep_size',0.2],
				['bias_mode',bias_mode],
				['bias_pot_fluctz',bias_pot_fluctz],
				['z_start',zref],
				['zstart_window',zstart],
				['numz_windows',zend],
				['zwindow_step_size',zwindow_step_size],
				['bias_pot_fluctH',bias_pot_fluctH],
				['H_start',Href],
				['Hstart_window',Hstart],
				['numH_windows',Hend],
				['Hwindow_step_size',Hwindow_step_size],
				['binsize_dist',0.1],
				['binsize_H',0.01],
				['hist_samp_freq',10000],
				['binsize_rdf',1.0],
				['approach_dir','advance'],
				['window_conf_interval',window_conf_interval],
				['window_antab_data_interval',window_antab_data_interval],
				['osd_print_interval',osd_print_interval]]                    #options are advance or recede
			
				fort_param=[[kappa,"!Bending Rigidity"],
					[czero ,"!czero"],
					[czero_pattern ,"!czero_pattern"],
					[czero_conc ,"!czero_conc"],
					[annulus_start ,"!annulus_start_ring"],
					[annulus_end ,"!annulus_end_ring"],
					[pressure,"!pressure"],
					[geom,"!System_Geom"],
					[height,"!Sinusoidal height"],
					[periodicity,"!Periodicity"],
					[blen_scale_factor,"!blen_scale_factor"],
					[memb_restart,res_file_num],
					[antigen_restart_flag,"!antigen_restart_flag"],
					[vermov_step_adjust,"!datainterval"],
					[fixed_frame_flag,"! is the frame fixed or not ?"]]
	
	
	
				os.mkdir('CC_INPUT_FILES')
				cpp_parameters(cpp_param)
				os.mkdir('FC_INPUT_FILES')
				fort_parameters(fort_param)
				os.mkdir('SYS_STATE')
				os.mkdir('SYS_STATE/MEMB')
				os.mkdir('SYS_STATE/NC-DATA')
				os.mkdir('SYS_STATE/ANTIGEN')		
				os.mkdir('SYS_STATE/DUMP')
				os.mkdir('SYS_STATE/RESTART')
	
	
				shutil.copyfile(codehome+'/'+exec_name,'RUNDIR/'+run_exe)
				shutil.copytree(codehome+'/CC_INPUT_FILES/SPH_RAD1_'+str(nab),'./CC_INPUT_FILES/SPH_RAD1_'+str(nab))
				shutil.copyfile(codehome+'/CC_INPUT_FILES/sphere'+str(nab)+'.in','./CC_INPUT_FILES/sphere'+str(nab)+'.in')
				if restart_memb==1:
					shutil.copyfile(memb_res_file,'./SYS_STATE/RESTART/membrane_state-'+str(nfiles)+'.dump')
		
				os.chdir('./RUNDIR')
				os.chmod(run_exe,stat.S_IRWXU)
				runcmd="nohup ./"+str(run_exe)+"&"
#				os.system(runcmd)
				os.chdir(subdir1)

			os.chdir(subdir)
		os.chdir(curdir)
	os.chdir(curdir_out)
