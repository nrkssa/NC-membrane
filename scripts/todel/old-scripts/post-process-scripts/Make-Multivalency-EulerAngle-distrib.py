#!/bin/env python
import os,sys,numpy,glob
ensdir=glob.glob('ENS*')
filename='NC-multivalency-angle-ENS-0_Frame-1.dat'
column=int(sys.argv[1])
if column==2:
	avgfile='Averaged-Multivalency.dat'
	binsize=1.0
	zmin,zmax=1.0,20.0
if column==3:
	avgfile='Averaged-EulerAngle1.dat'
	binsize=0.1
	zmin,zmax=0.0,2*numpy.pi
if column==4:
	avgfile='Averaged-EulerAngle2.dat'
	binsize=0.1
	zmin,zmax=0.0,2*numpy.pi
if column==5:
	avgfile='Averaged-EulerAngle3.dat'
	binsize=0.1
	zmin,zmax=0.0,2*numpy.pi

nbin=int(round((zmax-zmin)/binsize))
ensavg=numpy.zeros([nbin-1,3],'d')
ensavg[:,0],nzero_ens=numpy.arange(zmin,zmax-binsize,binsize),0
for ens in ensdir:
	outerdir=os.getcwd()
	os.chdir(ens)
	try:
		data=numpy.loadtxt('./RUNDIR/'+filename)
		print os.getcwd()+ '/RUNDIR/'+filename
		if (data.__len__()>0):
			fdata=data[:,column]
			hist,edges=numpy.histogram(fdata,numpy.linspace(zmin,zmax,num=nbin))
			binvalue=[0.5*(edges[i]+edges[i+1]) for i in range(0,edges.__len__()-1)]
			histarray=numpy.zeros([binvalue.__len__(),3],'d')

			if sum(hist)>0:
				histarray[:,0],histarray[:,1],histarray[:,2]=binvalue,hist,1.0*hist/sum(hist)
			else:
				histarray[:,0],histarray[:,1],histarray[:,2]=binvalue,hist,1.0*hist
	
			nzero_ens +=1
			ensavg[:,1]=ensavg[:,1]+histarray[:,1]
			ensavg[:,2]=ensavg[:,2]+histarray[:,2]
	except:
		print './RUNDIR/',filename,'does not contains any data - Not included into average'
	os.chdir(outerdir)

ensavg[:,1]=ensavg[:,1]/nzero_ens
ensavg[:,2]=ensavg[:,2]/nzero_ens


numpy.savetxt(avgfile,ensavg,delimiter='\t')
print '-------------------------------------------------------------------------------------'

