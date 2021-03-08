from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math


#
PI = 3.14
 
Sigma_h0 = 100.0  # [m]

D_h = 1.0  # [m2/s]

C0 = 1.0  # [kg/m2]

x = [0.0]
C = [2.5]

N_time = 5000
N_ring = 100

for i in range(N_ring-1):
	x.append( (i+1)*10.0 )
	C.append( 0.0 )


Sigma_h = Sigma_h0

print(x)

for i in range(N_ring):
	C[i] = C0 / (math.sqrt(2.0*PI)*Sigma_h) * math.exp(-0.5*( x[i]**2/Sigma_h**2))
print(C[0:50])

for t in range(1,N_time,1):
	
	print(t)
	Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * t)

	for i in range(N_ring):
		C[i] = C0 / (math.sqrt(2.0*PI)*Sigma_h) * math.exp( -0.5* ( x[i]**2/Sigma_h**2 ) ) 
	print(C[0:30])



