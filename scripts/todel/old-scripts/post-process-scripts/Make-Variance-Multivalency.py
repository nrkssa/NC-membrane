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
global data
ensdir=glob.glob('ENS*')
filename='NC-Antibody-Data-ENS-0_Frame-1.dat'
avgdata=numpy.zeros([20,4],'d')
for ens in ensdir:
	outerdir=os.getcwd()
	os.chdir(ens)
	input_file='./RUNDIR/'+filename
	print os.getcwd()+input_file
	linecount_cmd='wc -l '+input_file
	nlines=int(os.popen(linecount_cmd).read().split()[0])
	data=numpy.zeros([nlines,2],'d')
	n,curr_mval,start_row=0,0,0
	fp=open(input_file,'r')
	for f in fp:
		g=f.split()
		data[n][0],data[n][1]=float(g[2]),float(g[3])
		l=(g.__len__()-4)/2
		if n==0:
			curr_mval=l
			start_row=n
		else:
			if (l != curr_mval):
				end_row=n
				dataset=data[start_row:end_row]
				variance=compute_variance(dataset)
				avgdata[curr_mval][0] += variance[0]
				avgdata[curr_mval][1] += variance[1]
				avgdata[curr_mval][2] += sum(variance)
				avgdata[curr_mval][3] += 1
				curr_mval=l
				start_row=n
		n+=1
	fp.close()
	os.chdir(outerdir)

avgfile='./NC-Variance-Multivalency.dat'
fp=open(avgfile,'w')
for i in range(0,20):
	if (avgdata[i][3]>0):
		fp.write(str(i)+'\t'+str(avgdata[i][0]/avgdata[i][3])+'\t'+str(avgdata[i][1]/avgdata[i][3])+'\t'+str(avgdata[i][2]/avgdata[i][3])+'\n')
	else:
		fp.write(str(i)+'\t'+str(avgdata[i][0])+'\t'+str(avgdata[i][1])+'\t'+str(avgdata[i][2])+'\n')
fp.close()
print '----------------------------------------------------------------------------------------------------------------------------------------------------'

