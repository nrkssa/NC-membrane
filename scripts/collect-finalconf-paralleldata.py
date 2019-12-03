#!/bin/env python
import numpy 
import os,glob,sys
import copy
from operator import itemgetter

try:
	mode = sys.argv[1]

except:
	mode = 'last';


nclist,antlist,memblist=[],[],[]
winnumnc, winnumant,winnummemb = [],[],[]

#----------------------- NC data    --------------------------------------------

winlist=glob.glob('./SYS_STATE/NC-DATA/Nanocarrier-1-window-*-state-1.vtu')
winlsplit = [ w.split('-') for w in winlist]
windownum = [ int(w[4]) for w in winlsplit]
windows = sorted(windownum)

for w in windows:
	filelist = glob.glob('./SYS_STATE/NC-DATA/Nanocarrier-1-window-'+str(w)+'-state-*.vtu')
	filelsplit = [f.split('-') for f in filelist]

	for i in range(0,filelsplit.__len__()):
		a = filelsplit[i][-1]
		filelsplit[i].remove(a)
		filelsplit[i].append(a.split('.')[0])
		filelsplit[i].append(a.split('.')[1])
		filelsplit[i][-2] = int(filelsplit[i][-2])
	
	ncfilesort = sorted(filelsplit,key=itemgetter(-2))
	for i in range(0,ncfilesort.__len__()):
		a,b = ncfilesort[i][-2],ncfilesort[i][-1]
		ncfilesort[i].remove(a)
		ncfilesort[i].remove(b)
		ncfilesort[i].append(str(a)+'.'+str(b))
	
	ncfilelist = ['-'.join(ncfilename) for ncfilename in ncfilesort]

	if mode == 'all':
		for i in range(0,len(ncfilelist)):
			nclist.append(ncfilelist[i])
			winnumnc.append(w)
	else:
		nclist.append(ncfilelist[-1])
		winnumnc.append(w)


#-------------------------------------------------------------------

#----------------------- Memebrane data   --------------------------------------------
winlist=glob.glob('./SYS_STATE/MEMB/membrane-window-*-state-1.vtu')
winlsplit = [ w.split('-') for w in winlist]
windownum = [ int(w[2]) for w in winlsplit]
windows = sorted(windownum)
for w in windows:
	filelist = glob.glob('./SYS_STATE/MEMB/membrane-window-'+str(w)+'-state-*.vtu')
	filelsplit = [f.split('-') for f in filelist]

	for i in range(0,filelsplit.__len__()):
		a = filelsplit[i][-1]
		filelsplit[i].remove(a)
		filelsplit[i].append(a.split('.')[0])
		filelsplit[i].append(a.split('.')[1])
		filelsplit[i][-2] = int(filelsplit[i][-2])
	
	ncfilesort = sorted(filelsplit,key=itemgetter(-2))
	for i in range(0,ncfilesort.__len__()):
		a,b = ncfilesort[i][-2],ncfilesort[i][-1]
		ncfilesort[i].remove(a)
		ncfilesort[i].remove(b)
		ncfilesort[i].append(str(a)+'.'+str(b))
	
	ncfilelist = ['-'.join(ncfilename) for ncfilename in ncfilesort]

	if mode == 'all':
		for i in range(0,len(ncfilelist)):
			memblist.append(ncfilelist[i])
			winnummemb.append(w)
	else:
		memblist.append(ncfilelist[-1])
		winnummemb.append(w)

print memblist
#-------------------------------------------------------------------

#----------------------- Antigen data   --------------------------------------------
winlist=glob.glob('./SYS_STATE/ANTIGEN/Antigen-window-*-state-1.vtu')
winlsplit = [ w.split('-') for w in winlist]
windownum = [ int(w[2]) for w in winlsplit]
windows = sorted(windownum)
for w in windows:
	filelist = glob.glob('./SYS_STATE/ANTIGEN/Antigen-window-'+str(w)+'-state-*.vtu')
	filelsplit = [f.split('-') for f in filelist]

	for i in range(0,filelsplit.__len__()):
		a = filelsplit[i][-1]
		filelsplit[i].remove(a)
		filelsplit[i].append(a.split('.')[0])
		filelsplit[i].append(a.split('.')[1])
		filelsplit[i][-2] = int(filelsplit[i][-2])
	
	ncfilesort = sorted(filelsplit,key=itemgetter(-2))
	for i in range(0,ncfilesort.__len__()):
		a,b = ncfilesort[i][-2],ncfilesort[i][-1]
		ncfilesort[i].remove(a)
		ncfilesort[i].remove(b)
		ncfilesort[i].append(str(a)+'.'+str(b))
	
	ncfilelist = ['-'.join(ncfilename) for ncfilename in ncfilesort]


	if mode == 'all':
		for i in range(0,len(ncfilelist)):
			antlist.append(ncfilelist[i])
			winnumant.append(w)
	else:
		antlist.append(ncfilelist[-1])
		winnumant.append(w)

print antlist

#---------------------------------------- Write the data now ----------------------------------------
fp1=open('./loadfiles-memb.pvd','w')
fp1.write('<?xml version="1.0"?>\n')
fp1.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"> \n')
fp1.write('<Collection>\n')
n=0
for i in range(0,len(memblist)):
		print i,len(memblist),winnummemb[i],memblist[i]
		fp1.write('<DataSet timestep="'+str(i)+'" group="'+str(int(i/(winnummemb[i]+1)))+'"  part="'+str(winnummemb[i])+'" file="'+memblist[i]+'" />\n')
fp1.write('</Collection>\n')
fp1.write('</VTKFile>\n')
fp1.close()

fp1=open('./loadfiles-nc.pvd','w')
fp1.write('<?xml version="1.0"?>\n')
fp1.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"> \n')
fp1.write('<Collection>\n')
n=0
for i in range(0,len(nclist)):
		fp1.write('<DataSet timestep="'+str(i)+'" group="'+str(int(i/(winnumnc[i]+1)))+'"  part="'+str(winnumnc[i])+'" file="'+nclist[i]+'" />\n')
fp1.write('</Collection>\n')
fp1.write('</VTKFile>\n')
fp1.close()


fp1=open('./loadfiles-antigen.pvd','w')
fp1.write('<?xml version="1.0"?>\n')
fp1.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"> \n')
fp1.write('<Collection>\n')
n=0
for i in range(0,len(antlist)):
		fp1.write('<DataSet timestep="'+str(i)+'" group="'+str(int(i/(winnumant[i]+1)))+'"  part="'+str(winnumant[i])+'" file="'+antlist[i]+'" />\n')
fp1.write('</Collection>\n')
fp1.write('</VTKFile>\n')
fp1.close()
