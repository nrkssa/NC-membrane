#!/bin/env python
import os,sys,numpy,glob
ensdir=glob.glob('ENS*')
ensavg,nens_zero,nens=[],0,0
for ens in ensdir:
	nens+=1
	outerdir=os.getcwd()
	os.chdir(ens)
	rdffiles=glob.glob('./RUNDIR/rdf-Distribution-XY_block*')
	for rdf in rdffiles:
		print 'averaging ',ens,rdf
		try:
			data=numpy.loadtxt(rdf)
			if (data.__len__()>0):
				print nens_zero
				if(nens_zero==0):
					data_sum=data
					nens_zero +=1
				else:
					data_sum+=data
					nens_zero +=1
		except:
			print rdf,'does not contains any data - Not included into average'
	os.chdir(outerdir)
	fdata=data_sum
avgfile='Avg-rdf-distribution.dat'
numpy.savetxt(avgfile,fdata/nens_zero,delimiter='\t')
print '-------------------------------------------------------------------------------------'

