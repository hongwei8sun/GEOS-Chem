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
FILEDIR = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_gc_timing_Dr200m/'

model_txt  = np.loadtxt(FILEDIR+'Plume_theta_max_min_radius.txt')

#concnt_model = model_txt[t,:]

#print(len(model_txt[t,:])) # 

# guassian ------------------------------
PI = 3.14

Sigma_h0 = 1000.0  # [m]
Sigma_v0 = 1000.0
Dr = 200.0

D_h = 1.0  # [m2/s]
D_v = 1.0

C0 = 100.0* PI*Sigma_h0*Sigma_v0 # [kg/m3 * m2 = kg/m]

x = [0.5*Dr]
z = [0.5*Dr]
concnt_gaus = [0.0]

N_ring = 25

for i in range(N_ring-1):
        x.append( (i+1.5)*Dr + 1.0 )
        z.append( (i+1.5)*Dr + 1.0 )
        concnt_gaus.append( 0.0 )


Sigma_h = Sigma_h0
Sigma_v = Sigma_v0


for j in range(100):
	t = j*3600	
	Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * t)
	Sigma_v = math.sqrt(Sigma_v0**2 + 2.0 * D_v * t)
	
	for i in range(N_ring):
		concnt_gaus[i] = C0 / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i]**2/Sigma_h**2 + z[i]**2/Sigma_v**2 ))
	
	
	
	# plot -----------------------------------
	
	plt.figure(figsize=(7,8))
		
	plt.plot(x, model_txt[t,0:25], 'b-', x, concnt_gaus, 'r--')
	
	plt.title( 'time = '+str(j)+' hour' )

	#plt.ylim((0.0, 55.0))
	#plt.yscale('log')
	
	#plt.yticks(y)

	plt.savefig(str(j)+'_xy.png')
	plt.clf()
	plt.cla()



