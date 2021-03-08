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
Sigma_v0 = 100.0

D_h = 1.0  # [m2/s]
D_v = 1.0

C0 = 100.0* PI*Sigma_h0*Sigma_v0 # [kg/m3 * m2 = kg/m]

x = [1.0]
z = [1.0]
C = [0.0]

N_time = 5000
N_ring = 500

D_r = 2.0

for i in range(N_ring-1):
	x.append( (i+1)*D_r + 1.0 )
	z.append( (i+1)*D_r + 1.0 )
	C.append( 0.0 )


Sigma_h = Sigma_h0
Sigma_v = Sigma_v0


file1 = open("Data_gaussian_2m.txt","w")

for i in range(N_ring):
	C[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))

print(C[0:250])

file1.write(str(C[0:N_ring]))
file1.write('\n')

for t in range(1,N_time,1):
	
	#print(t)
	Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * t)
	Sigma_v = math.sqrt(Sigma_v0**2 + 2.0 * D_v * t)
	
	for i in range(N_ring):
		C[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))
	#print(C)
	file1.write(str(C))
	file1.write('\n')



