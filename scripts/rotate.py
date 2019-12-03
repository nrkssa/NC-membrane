#!/bin/env python
from numpy import *
rot_mat = zeros([3,3],'d')
vect = zeros([3,1],'d')+1.0
print (vect)

phi,psi,theta = 1.0,0.5,0.5
cosphi = cos(phi);
cospsi = cos(psi);	
costheta = cos(theta);
sinphi = sin(phi);
sinpsi = sin(psi);
sintheta = sin(theta);
	
rot_mat[0][0] = cosphi*cospsi - sinphi*costheta*sinpsi;
rot_mat[1][0] = sinphi*cospsi + cosphi*costheta*sinpsi;
rot_mat[2][0] = sintheta*sinpsi;		                    
rot_mat[0][1] = -cosphi*sinpsi - sinphi*costheta*cospsi;
rot_mat[1][1] = cosphi*costheta*cospsi -sinphi*sinpsi;
rot_mat[2][1] = sintheta*cospsi;
rot_mat[0][2] = sinphi*sintheta;
rot_mat[1][2] = -cosphi*sintheta;
rot_mat[2][2] = costheta;

nvect = dot(rot_mat,vect)
print (nvect,sum(nvect**2))

phi,psi,theta = 1.0, 0.5,0.5
cosphi = cos(phi);
cospsi = cos(psi);							
costheta = cos(theta);
sinphi = sin(phi);
sinpsi = sin(psi);
sintheta = sin(theta);
	
rot_mat[0][0] = cosphi*cospsi - sinphi*costheta*sinpsi;
rot_mat[1][0] = sinphi*cospsi + cosphi*costheta*sinpsi;
rot_mat[2][0] = sintheta*sinpsi;		                    
rot_mat[0][1] = -cosphi*sinpsi - sinphi*costheta*cospsi;
rot_mat[1][1] = cosphi*costheta*cospsi -sinphi*sinpsi;
rot_mat[2][1] = sintheta*cospsi;
rot_mat[0][2] = sinphi*sintheta;
rot_mat[1][2] = -cosphi*sintheta;
rot_mat[2][2] = costheta;

nvect1 = dot(rot_mat.transpose(),nvect)
print (nvect1,sum(nvect1**2))
