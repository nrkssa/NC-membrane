def compute_variance(dataset):
	ncol,nrow=dataset[0].__len__(),dataset.__len__()
	variance=numpy.zeros(ncol,'d')
	dataset2=dataset**2
	for i in range(0,ncol):
		avg=sum(dataset[:,i])/nrow
		avg2=sum(dataset2[:,i])/nrow
		variance[i]=avg2-avg**2
	return variance


#!/bin/env python
import os,sys,numpy,glob
boxl=int(sys.argv[1])
nc_rad=100
binsize=0.5
nbin=int(2*nc_rad/binsize)
bindata=numpy.zeros([nbin,nbin],'d')
ensdir=glob.glob('ENS*')
filename='NC-Antibody-Data-ENS-0_Frame-1.dat'
for ens in ensdir:
	outerdir=os.getcwd()
	os.chdir(ens)
	input_file='./RUNDIR/'+filename
	print os.getcwd()+input_file
	fp=open(input_file,'r')
	for f in fp:
		g=f.split()
		l=(g.__len__()-4)/2
		if l>0:
			g1=[float(x) for x in g]
			for i in range(0,l):
				delx,dely=g1[2*i+4]-g1[2],g1[2*i+5]-g1[3]
				delx -= round(delx/boxl)*boxl
				dely -= round(dely/boxl)*boxl
				binx = int((delx+nc_rad)/binsize)
				biny = int((dely+nc_rad)/binsize)
				bindata[binx][biny] += 1

	fp.close()
	os.chdir(outerdir)

nsum=numpy.sum(bindata)
print nsum,nbin
avgfile='./NC-Bond-Distribution.dat'
fp=open(avgfile,'w')
for i in range(0,nbin):
	for j in range(0,nbin):
		fp.write(str((i*binsize-nc_rad))+'\t'+str((j*binsize-nc_rad))+'\t'+str(float(bindata[i][j]/nsum))+'\n')
	fp.write('\n')
fp.close()
print '----------------------------------------------------------------------------------------------------------------------------------------------------'

