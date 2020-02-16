#
# Create the initial concentration for Plume model
#
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
 
Sigma_h0 = 1000.0  # [m]
Sigma_v0 = 1000.0
Dr       = 100.0

D_h = 1.0  # [m2/s]
D_v = 1.0

C0 = 100.0*PI*Sigma_h0*Sigma_v0  # [kg/m]

print(C0)

x = [0.5*Dr]
z = [0.0]
C = [0.0]

N_time = 3600*24	# For 10 days
N_ring = 150

for i in range(N_ring-1):
	x.append( (i+1.5)*Dr )
	z.append( 0.0 )
	C.append( 0.0 )


Sigma_h = Sigma_h0
Sigma_v = Sigma_v0
print(x)

for i in range(N_ring):
	C[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))
print(C[0:N_ring])

# quit()

###
#print(' ')
#print(N_time)
#Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * N_time)
#Sigma_v = math.sqrt(Sigma_v0**2 + 2.0 * D_v * N_time)
#
#for i in range(N_ring):
#        C[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))
#print(C[0:30])
#quit()
###


for t in range(1,N_time,600):
	
	print(' ')
	print(t)
	Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * t)
	Sigma_v = math.sqrt(Sigma_v0**2 + 2.0 * D_v * t)

	for i in range(N_ring):
		C[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))
	print(C[0:30])


