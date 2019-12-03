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
import numpy,sys,os,time,shutil,stat,subprocess,glob,random
start_ens,end_ens = int(sys.argv[1]),int(sys.argv[2])


if sys.argv.__len__() <3:
	print 'Enter the number for starting and ending ensembles'
	sys.exit('Exiting')
	
# fortran variables
antigen_density  = 0.001                  # 0.002 -> 2000 antigens per micron^2
kappa,pressure = 20,0.0                   # kappa, pressure, and the excess surface area parameter
blen_fac=[0.25]
dcd_frequency=100

n_processor=16

czero,czero_pattern,czero_conc = 0.0, 'CIRCULAR', 0.0
if(czero_conc == 0.0):
	czero=0.0
	
annulus_start,annulus_end = 3,10
antigen_restart_flag = '.False.'
fixed_frame_flag = '.False.'
restart_memb=0                                # set this flag to 1 if the membrane is to be restarted


nc_size,num_antibody = 50,[162]
antigen_size,antibody_size = 19.0,15.0
#z_start_min = nc_size + antigen_size + antibody_size +0.5
z_start_min = 84.5


#C++ variables
num_moves=600000000
N,init_bond_length= 50,10                     # grid size and initial grid spacing  
geom='PLANAR'
periodicity=2
height=20
flexural_rigidity=7000
vermov_step_adjust=10000
bias_mode ='Z'
bias_pot_fluctz,z_start,zstart_window,numz_windows,zwindow_step_size=1.0,z_start_min,0,1,0.1
bias_pot_fluctH,H_start,Hstart_window,numH_windows,Hwindow_step_size=0.0,0.024,0,48,0.0015
window_conf_interval=num_moves/100
window_antab_data_interval=num_moves/60000
osd_print_interval=num_moves/100000

source_dir=os.getcwd()
startconf_dir=source_dir+'/../..'
codehome=os.getcwd()
exec_name='NC_MEMB_fixed'

if (bias_mode == 'H'):
	dirname1='Fixed-Equm_Hbias_N-'+str(N)+"_bl-"+str(init_bond_length)
elif (bias_mode == 'Z'):
	dirname1='Fixed-Equm_Rbias_N-'+str(N)+"_bl-"+str(init_bond_length)
else:
	dirname1='Fixed-Equm-N'+str(N)+"_bl-"+str(init_bond_length)

dirname1_rs='Memb-N'+str(N)+'_bl-'+str(init_bond_length)
	
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
	num_antigens=int((N*init_bond_length)**2*antigen_density)           # computed from the density predefined in the system (2000 antigens per micron^2)
	
		
	for nab in num_antibody:
		curdir=os.getcwd()
		dirname4='AB-'+str(nab)+'_nantigen-'+str(num_antigens)
		subprocess.call(["mkdir","-pv",dirname4])
		os.chdir(dirname4)
		subdir=os.getcwd()


		for ens in range(start_ens,end_ens):
			run_exe='equm-'+str(nab)+'-N'+str(N)+'E-'+str(ens)
			dirname5='ENS-'+str(ens)
			subprocess.call(["mkdir","-pv",dirname5])
			os.chdir(dirname5) 
		
			dirname50='ENS-'+str(ens)
			os.mkdir('RUNDIR')
			memb_restart,res_file_num=".False.",0
			if restart_memb ==1:
				restart_dir=startconf_dir+'/'+dirname1_rs+'/'+dirname10+'/'+dirname40+'/'+dirname5+'/SYS_STATE/DUMP'
				print restart_dir
				nfiles=glob.glob(restart_dir+'/*.dump').__len__()-random.randint(1,6)
				memb_restart,res_file_num=".True.",nfiles
				memb_res_file=restart_dir+'/membrane_state-'+str(nfiles)+'.dump'
			print memb_restart,res_file_num

			if bias_mode == 'Z':
				zref,zstart,zend = z_start,zstart_window,numz_windows
				Href,Hstart,Hend = 0.0,0,0
			elif bias_mode == 'H':
				Href,Hstart,Hend = H_start,Hstart_window,numH_windows
				zref,zstart,zend = z_start,0,0
			elif bias_mode == 'N':
				zref,zstart,zend = z_start,0,0
				Href,Hstart,Hend = 0.0,0,0


			cpp_param=[['periodic_box_length_height_ratio',2.0],
			['link_cell_length',nc_size],
			['num_vesicles',1],
			['radius',nc_size],
			['antibody_size',antibody_size],
			['num_antigens',num_antigens],
			['antigen_radius',1.5],
			['num_antibody',nab],
			['flexural_rigidity',flexural_rigidity],
			['eqm_bond_dist',antigen_size],
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
				[fixed_frame_flag,"! is the frame fixed or not ?"],
				[2.0,"! shadow_size"]]

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
			os.chdir('../')

                        fp=open('subscript.pbs','w')
                        fp.write('#!/bin/bash \n')
                        if bias_mode == 'Z':
                                fp.write('#PBS -N '+str(kappa)+'-'+str(blen_scale_factor)+'-'+str(zref)+'-'+str(nab)+'-'+str(ens)+'\n')
                        elif bias_mode == 'H':
                                fp.write('#PBS -N '+str(kappa)+'-'+str(blen_scale_factor)+'-'+str(Href)+'-'+str(nab)+'-'+str(ens)+'\n')
                        else:
                                fp.write('#PBS -N '+str(kappa)+'-'+str(blen_scale_factor)+'-equm-'+str(nab)+'-'+str(ens)+'\n')
                        fp.write('#PBS -M ram.n.krishnan@gmail.com \n')
                        fp.write('#PBS -o output.dat\n')
                        fp.write('#PBS -l nodes=1:ppn='+str(n_processor)+',pmem=1500mb \n')
			fp.write('#PBS -q long\n')
                        fp.write('#PBS -l walltime=36:00:00\n')
                        fp.write('module load openmpi-1.6.4-gcc \n')
                        fp.write('cd $PBS_O_WORKDIR\n')
                        fp.write('jobid=`echo $PBS_JOBID|cut -f1 -d .`\n')
                        fp.write('rundir=/scratch-local/job-$jobid\n')
                        fp.write('mkdir $rundir\n')
                        fp.write('mv * $rundir\n')
                        fp.write('cd $rundir\n')
			fp.write('cd RUNDIR\n')
                        fp.write("qstat -f >>$PBS_O_WORKDIR/jobdetails.dat\n")
                        fp.write('mpirun -np ' + str(n_processor) + ' ' + str(run_exe) + '\n')
                        fp.write('cd ../\n')
                        fp.write('mv * $PBS_O_WORKDIR\n')
                        fp.write('rmdir $rundir\n')
                        fp.close()
                        runcmd="qsub subscript.pbs"
                        os.system(runcmd)

			os.chdir(subdir)
		os.chdir(curdir)
	os.chdir(curdir_out)
