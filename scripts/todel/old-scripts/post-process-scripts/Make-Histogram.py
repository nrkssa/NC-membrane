#!/bin/env python
import os,sys,numpy,glob
ensdir=glob.glob('ENS*')
filename=sys.argv[1]
binsize=1.0
zmin,zmax=1.0,20.0
ensavg=numpy.zeros([zmax-zmin-1,3],'d')
ensavg[:,0]=numpy.arange(zmin,zmax-1,binsize)
for ens in ensdir:
	outerdir=os.getcwd()
	os.chdir(ens)
	nbin=(zmax-zmin)/binsize
	fp=open('./RUNDIR/'+filename,'r')
	print 'Making histogram for ',os.getcwd()+'/RUNDIR/'+filename
	data=[line.strip().split() for line in fp]
	fp.close()
	fdata=numpy.array(data[2:-1],dtype='d')
	hist,edges=numpy.histogram(fdata,numpy.linspace(zmin,zmax,num=nbin))
	binvalue=[0.5*(edges[i]+edges[i+1]) for i in range(0,edges.__len__()-1)]
	histarray=numpy.zeros([binvalue.__len__(),3],'d')
	if sum(hist)>0:
		histarray[:,0],histarray[:,1],histarray[:,2]=binvalue,hist,1.0*hist/sum(hist)
	else:
		histarray[:,0],histarray[:,1],histarray[:,2]=binvalue,hist,1.0*hist

	histfile='hist-'+filename
	numpy.savetxt(histfile,histarray,delimiter='\t')
	ensavg[:,1]=ensavg[:,1]+histarray[:,1]
	ensavg[:,2]=ensavg[:,2]+histarray[:,2]
	os.chdir(outerdir)

ensavg[:,1]=ensavg[:,1]/ensdir.__len__()
ensavg[:,2]=ensavg[:,2]/ensdir.__len__()
avgfile='Avg-'+filename
numpy.savetxt(avgfile,ensavg,delimiter='\t')
print '-------------------------------------------------------------------------------------'

