#!/bin/env python
import os,sys,numpy,glob
ensdir=glob.glob('ENS*')
filename=sys.argv[1]
binsize=1.0
zmin,zmax=1.0,20.0
ensavg=numpy.zeros([zmax-zmin-1,3],'d')
ensavg[:,0],nzero_ens=numpy.arange(zmin,zmax-1,binsize),0
nbin=(zmax-zmin)/binsize
for ens in ensdir:
	print ens
	outerdir=os.getcwd()
	os.chdir(ens)
	try:
		data=numpy.loadtxt('./RUNDIR/'+filename)
		print os.getcwd()
		if (data.__len__()>0):
			fdata=data[:,1]
			hist,edges=numpy.histogram(fdata,numpy.linspace(zmin,zmax,num=nbin))
			binvalue=[0.5*(edges[i]+edges[i+1]) for i in range(0,edges.__len__()-1)]
			histarray=numpy.zeros([binvalue.__len__(),3],'d')
			if sum(hist)>0:
				histarray[:,0],histarray[:,1],histarray[:,2]=binvalue,hist,1.0*hist/sum(hist)
			else:
				histarray[:,0],histarray[:,1],histarray[:,2]=binvalue,hist,1.0*hist
	
			nzero_ens +=1
			histfile='hist-'+filename
			numpy.savetxt(histfile,histarray,delimiter='\t')
			ensavg[:,1]=ensavg[:,1]+histarray[:,1]
			ensavg[:,2]=ensavg[:,2]+histarray[:,2]
	except:
		print './RUNDIR/',filename,'does not contains any data - Not included into average'
	os.chdir(outerdir)

ensavg[:,1]=ensavg[:,1]/nzero_ens
ensavg[:,2]=ensavg[:,2]/nzero_ens
avgfile='Avg-'+filename
numpy.savetxt(avgfile,ensavg,delimiter='\t')
print '-------------------------------------------------------------------------------------'

