class antibody:
	def __init__(self,ntype):
		
	

class nanocarrier():
	def __init__(self,index):
		self.ncno = index
		self.ab_type =1
		self.ab_radius=0.0
		self.ab_length=0.0
		
		self.nc_shape='none'
		self.nc_radius=0.0
		self.nc_num_ab=0
		self.nc_rstart=0
		self.nc_hstart=0
		self.ab_type_conc=0
		self.ab_mode='none'
		
		self.bias_mode='N'
		self.biasstrength=1.0
		self.biasref=0.0
		self.binsize=1.0
		self.biasdir='approach'
		self.TI_lambda='1.0'
		
		self.ab={'ab_type':self.ab_type,'ab_radius':self.ab_radius,
			      'ab_length':self.ab_length}
		self.nc={'nc_shape':self.nc_shape,'nc_radius':self.nc_radius,
				  'nc_num_ab':self.nc_num_ab,'nc_rstart':self.nc_rstart,
				  'nc_hstart':self.nc_hstart,'ab_type_conc':self.ab_type_conc,
				  'ab_mode':self.ab_mode
				  }
		self.bias={'bias_mode':self.bias_mode,'biasstrength':self.biasstrength,
			        'biasref':self.biasref,'binsize':self.binsize,
			        'biasdir':self.biasdir,'TI_lambda':self.TI_lambda
			        }
	
	def write_ncfile(self):
		filename='PARAMETERS/nc-'+str(self.ncno)+'.ncin'
		fp=open(filename,'w')
		
		fp.write('<antibody> \n')
		for a in self.ab:
			fp.write(a+', '+str(self.ab[a])+'\n')
		fp.write('</antibody> \n\n')
		
		fp.write('<nanocarrier> \n')
		for a in self.nc:
			fp.write(a+', '+str(self.nc[a])+'\n')
		fp.write('</nanocarrier> \n\n')
		
		fp.write('<bias> \n')
		for a in self.bias:
			fp.write(a+', '+str(self.bias[a])+'\n')
		fp.write('</bias> \n\n')
		fp.close()

#!/bin/env python
import numpy
nc=nanocarrier(0)
nc.write_ncfile()

			