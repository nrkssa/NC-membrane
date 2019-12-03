#!/bin/env python                                          

# code to convert an intel code compiled with -assume nounderscore option to a code running on 
# gnu compilers and compiled with -fno-underscoring option

import glob,os
intel_module="module_"
gnu_module="__module_"
intel_mod="_mp_"
gnu_mod="_MOD_"

printflag=0
curdir=os.getcwd()

os.chdir('CC_SOURCE/')                                        #convert all CC_SOURCE files to Intel version
ccsource=glob.glob('*.cpp')
print ('\n*************************************************')
for ccode in ccsource:
	fp=open(ccode,'r')
	newname='../../../Parallel-MPI-Intel/multi-particles-revamp/CC_SOURCE/'+ccode
	fp1=open(newname,'w')
	for f in fp:
		fnew=f
		if f.find(gnu_module)>=0:
			fnew=fnew.replace(gnu_module,intel_module);
		if f.find(gnu_mod)>=0:
			fnew=fnew.replace(gnu_mod,intel_mod);
		fp1.write(fnew)
	fp.close()
	fp1.close()
	print ('completed (GNU->INTEL) : ',ccode)
os.chdir(curdir)
curdir=os.getcwd()

print ('\n*************************************************')
os.chdir('CC_INCLUDE/')
cinclude=glob.glob('*.h')
for inccode in cinclude:
	fp=open(inccode,'r')
	newname='../../../Parallel-MPI-Intel/multi-particles-revamp/CC_INCLUDE/'+inccode
	fp1=open(newname,'w')
	for f in fp:
		fnew=f
		if f.find(gnu_module) >= 0:
			fnew=fnew.replace(gnu_module,intel_module);
		if f.find(gnu_mod) >= 0:
			fnew=fnew.replace(gnu_mod,intel_mod);
		fp1.write(fnew)
	fp.close()
	fp1.close()
	print ('completed (GNU->INTEL) : ',inccode)
os.chdir(curdir)
curdir=os.getcwd()

print ('\n*************************************************')
os.chdir('FC_SOURCE/')                                        #convert all CC_SOURCE files to Intel version
fcsource=glob.glob('*.f')
for fcode in fcsource:
	fp=open(fcode,'r')
	newname='../../../Parallel-MPI-Intel/multi-particles-revamp/FC_SOURCE/'+fcode
	fp1=open(newname,'w')
	for f in fp:
		fp1.write(f)
	fp.close()
	fp1.close()
	print ('completed (GNU->INTEL) : ',fcode)
os.chdir(curdir)
curdir=os.getcwd()

print ('\n*************************************************')
os.chdir('FC_INCLUDE/')
finclude=glob.glob('*.h')
for infcode in finclude:
	fp=open(infcode,'r')
	newname='../../../Parallel-MPI-Intel/multi-particles-revamp/FC_INCLUDE/'+infcode
	fp1=open(newname,'w')
	for f in fp:
		fp1.write(f)
	fp.close()
	fp1.close()
	print ('completed (GNU->INTEL) : ',infcode)
os.chdir(curdir)
