from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math
#pandas

#FILEDIR = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_standard_tracer/'
FILEDIR  = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_gc_timing_plume_100m/'
FILEDIR2 = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_gc_timing_plume_100m/'
FILEDIR3 = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_gc_timing_plume_100m/'

model_txt  = np.loadtxt(FILEDIR+'Plume_theta_max_min_radius_50000.txt')
model2_txt = np.loadtxt(FILEDIR2+'Plume_theta_max_min_radius_50000.txt')
model3_txt = np.loadtxt(FILEDIR3+'Plume_theta_max_min_radius_50000.txt')

#concnt_model = model_txt[t,:]

#print(len(model_txt[t,:])) # 


#
# guassian ------------------------------
PI = 3.14

Sigma_h0 = 100.0  # [m]
Sigma_v0 = 100.0

D_h = 1.0  # [m2/s]
D_v = 1.0

C0 = 100.0* PI*Sigma_h0*Sigma_v0  # [kg/m3 * m2 = kg/m]

x = [50.0]
z = [50.0]
concnt_gaus = [0.0]

N_time = 5000
N_ring = 500
D_ring = 100.0 # [m]

for i in range(N_ring-1):
        x.append( (i+1)*D_ring + x[0] )
        z.append( (i+1)*D_ring + x[0] )
        concnt_gaus.append( 0.0 )


Sigma_h = Sigma_h0
Sigma_v = Sigma_v0

#concnt_gaus = [0.0]


for j in range(42):
	t = j*100	
	Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * t)
	Sigma_v = math.sqrt(Sigma_v0**2 + 2.0 * D_v * t)
	
	print(Sigma_h)

	for i in range(N_ring):
		concnt_gaus[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))
	
	
	
	# plot -----------------------------------
	
	plt.figure(figsize=(7,8))
		
	plt.plot(x[0:50], concnt_gaus[0:50], 'b-', x[0:50], model_txt[t,0:50], 'r--')
	
	plt.title( 'time = 100s * '+str(j) )

	plt.legend(('Gaussian Analytic', 'GC'),loc='upper right')

	#plt.ylim((0.0, 55.0))
	#plt.yscale('log')
	
	#plt.yticks(y)

	plt.savefig(str(j)+'_xy.png')
	plt.clf()
	plt.cla()



