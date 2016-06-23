# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 19:27:03 2015

@author: Inom Mirzaev

Wraps the C-code for cell division process in 3D
"""

from __future__ import division

import ctypes
import mayavi.mlab as mlab

mydll = ctypes.CDLL('/home/inom/Dropbox/Research/floc_growth/latest_codes/cell_division.so')
    
import numpy as np
import numpy.ctypeslib as npct
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

start = time.time()

array_2d_double = npct.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
array_1d_double = npct.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')

cell_division= mydll.cell_division

#Argument types for c-code
cell_division.argtypes = ( ctypes.c_int, 	#length of allowed radius vector
                           ctypes.c_int, 	#length of maximum location matrix
						  ctypes.c_int, 	#number of cells at the beginning 
						   ctypes.c_int, 	#number of steps
						   ctypes.c_int, 	#number of phi discretization points 
						   ctypes.c_int, 	#number of theta discretization points 
                           array_2d_double, #location matrix
                           array_1d_double, #allowed radius vector
                           array_1d_double  #total volume at each step
                           )
                           
#Return type for c-code                           
cell_division.restype  = ctypes.c_int

def mycell_division(loc_mat , allowed_rad, numsteps, volume, num_phi, num_theta ):


    numcells = cell_division( len(allowed_rad) , len(loc_mat) , 3, numsteps, num_phi, num_theta, loc_mat, allowed_rad, volume )
    
    return (loc_mat, volume, numcells)



numsteps = 10

mat_dim = numsteps**3

vol = np.zeros( numsteps, np.float64 )
num_phi   = 20
num_theta = 40

loc_mat = np.zeros((mat_dim,5), dtype=np.float64)

loc_mat[0] = [0,0,0, 0.5,1]
loc_mat[1] = [0,1,0, 0.5,1]
loc_mat[2] = [0,0,1, 0.5,1]


allowed_rad = np.arange(0.5, 0.40, -0.01)


(loc_mat, vol, numcells)=mycell_division(loc_mat , allowed_rad,  numsteps, vol, num_phi, num_theta)

loc_mat = loc_mat[0:numcells,:]

np.savetxt('loc_mat_40.out', loc_mat)
np.savetxt('volumebysize_40.out' , vol)

end = time.time()

print 'Number of cells at the end ' + str(numcells)

print 'Time elapsed ' + str( round( ( end - start ) , 2 ) ) + ' second'


plt.close('all')    
 
growth = vol[range( 1 , numsteps ) ] - vol[ range( numsteps - 1 ) ]   



plt.figure(0)

sfac = 1 / vol[-2] 
plt.plot( sfac*vol[ range( numsteps -1 ) ] , sfac*growth , linewidth=2  )

(x1, x2, y1, y2) = plt.axis()
x2 = np.max(sfac * vol[-2])
plt.axis( (x1, x2, y1, y2) )


plt.xlabel('$x$')
plt.ylabel('$g(x)$')   
plt.show()
plt.savefig('growth_rate_40steps.png', dpi=400)

"""
mlab.close(all=True)


mlab.figure(  bgcolor=(1,1,1) )

mlab.points3d( loc_mat[:, 0], loc_mat[:, 1], loc_mat[:, 2] , 
               0.5*np.ones( len( loc_mat ) ), scale_factor=2.0 , 
               resolution=20 )
               
mlab.view(distance = 75 )"""







