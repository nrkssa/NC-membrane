def square_lattice(L,N):
	"make a NXN square lattice"
	data=numpy.zeros([(N+1)*(N+1),3],'d')
	binx,biny=L/N,L/N
	n=0
	for i in range(0,N+1):
		for j in range(0,N+1):
			data[n,0],data[n,1] = i*binx,j*biny
			n += 1
	
	return [data,min(binx,biny)]

def random_lattice(L,N,partsize):
	"make a NXN square lattice"
	tlen = (N+1)**2
	data=numpy.zeros([tlen,3],'d')

	if tlen*numpy.pi*partsize**2>L**2:
		cutoff = numpy.sqrt(L**2/(numpy.pi*tlen))
		print 'resizing overlap cutoff to ',cutoff
	else:
		cutoff=partsize

	n=0
	data[n,0],data[n,1]=random.randrange(0,L),random.randrange(0,L)
	n += 1
	while n<tlen:
		data[n,0],data[n,1]=random.randrange(0,L),random.randrange(0,L)
		if selfavoidance(data[n-1,:],data[n,0],cutoff):
			n += 1
	
	return [data,partsize]

def selfavoidance(data,target,cutoff):
	for d in data:
		dis = numpy.sum((d-target)**2)
		if dis<cutoff:
			return False
	return True

def compute_gofr(data,L,binsize,nsample,gofr):
	"compute gofr"
	for i in range(0,len(data)-1):
		for j in range(i+1,len(data)):
			rvec = data[i,:]-data[j,:]
			rvec[0] -= round(rvec[0]/L)*L
			rvec[1] -= round(rvec[1]/L)*L
			dis=numpy.sqrt(numpy.sum(rvec**2))
		 	if dis<0.5*L:
				bnum=int(dis/binsize)
#				print i,j,dis,bnum
#				print data[i,:]
#				print data[j,:]
#				print '------------------------------>\n'
				gofr[bnum][1] += 2
	
	nsample = nsample +1
	return nsample
	

def normalize_gofr(nparticles,nsample,gofr,binsize,rho):
	print nparticles,nsample
	nbins = len(gofr)
	gofr1 = numpy.zeros([nbins,2],'d')
	for i in range(0,len(gofr)):
		gofr1[i,0] = (i+0.5)*binsize
		da = ((i+1)**2-i**2)*binsize**2
		nid = numpy.pi*da*rho
		print i,da,nid,gofr[i,1]
		gofr1[i,1] = gofr[i,1]/(nid*nsample*nparticles)
	return gofr1

	


#!/bin/env python
import numpy,os,glob,sys
import matplotlib.pyplot as plt
import brewer2mpl as mpl
import random
import matplotlib.gridspec as gs

fig = plt.figure(figsize=(3.0,2.0))
grid = gs.GridSpec(2,1,wspace=0.2,hspace=0.2)
ax1 = fig.add_subplot(grid[0:1,0:1])
ax2 = fig.add_subplot(grid[1:2,0:1])
colors = mpl.get_map('Set3','qualitative',5).mpl_colors
binsize=1.0
filestring=sys.argv[1]
ofile='gofr-'+str(filestring.split('/')[1])+'.eps'

filenames=glob.glob(filestring)


i,nsample,leng=0,0,[]
data,nant,areas=[],[],[]

for filename in filenames:
	metadata=open(filename,'r').readline().split()
	L,area=float(metadata[1]),float(metadata[3])
	data.append(numpy.loadtxt(filename,skiprows=1))
	leng.append(L)
	areas.append(area)
	nant.append(len(data[-1]))
	ax1.scatter(data[-1][:,0],data[-1][:,1],color=colors[i%len(colors)])
	i += 1

L=max(leng)
print L
nbins=int(L/binsize)
gofr = numpy.zeros([nbins,2],'d')
for d in data:
	print 'generating gofr for ',filename
	nsample=compute_gofr(d,L,binsize,nsample,gofr)



	
rho=min(nant)/area
gofrn = normalize_gofr(min(nant),nsample,gofr,binsize,rho)
print 'plotting gofr'
ax2.plot(gofrn[:,0],gofrn[:,1],lw=1,ms=0,color='b')
print '-------------------------------->\n'
plt.savefig(ofile,dpi=300)
plt.show()

