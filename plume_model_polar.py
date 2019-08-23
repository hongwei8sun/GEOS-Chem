from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math
#pandas


# guassian ------------------------------
PI = 3.14

Sigma_h0 = 1000.0  			# [m]
Sigma_v0 = 1000.0
Dr = 100.0  				#[m]

D_h = 1.0  				# diffusion coefficient, [m2/s]
D_v = 1.0

Q = 100.0* PI*Sigma_h0*Sigma_v0 	# [kg/m3 * m2 = kg/m]

x = [0.5*Dr]
z = [0.5*Dr]
concnt_gaussian = [0.0]			# concentration from gaussian analytical solution

N_ring = 50

for i_ring in range(N_ring-1):
        x.append( (i_ring+1.5)*Dr )
        z.append( (i_ring+1.5)*Dr )
        concnt_gaussian.append( 0.0 )


Sigma_h = Sigma_h0
Sigma_v = Sigma_v0


###
# let z=0.0 when calculate concentration distribution !!!
###

# initial concentration distribution from gaussion analytical results---------
for i_ring in range(N_ring-1):
        concnt_gaussian[i_ring] = Q / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i_ring]**2/Sigma_h**2 + 0.0**2/Sigma_v**2 ))


# plume modela & xy model -------------------------------
Dif          = 1.0     				# equal to D_h and D_v, [m2/s]
r            = [0.5*Dr]				# equal to x and z
concnt_model = [0.0]				# concentration from model results
concnt_old   = [0.0]

for i_ring in range(N_ring-1):
    r.append( (i_ring+1.5)*Dr )
    concnt_model.append( 0.0 )
    concnt_old.append( 0.0 )

# use gaussian distribution as model initial concentration:
for i_ring in range(0,N_ring-1,1):
    concnt_model[i_ring] = concnt_gaussian[i_ring]

Dt = 1.0 					#[s]
for i_time in range(1,99999,1): 		# [s]
    for i_ring in range(0,N_ring-1,1):
        concnt_old[i_ring] = concnt_model[i_ring]

    # calculate new concentration after dilusion
    concnt_model[0] = concnt_old[0] + Dt*Dif*2.0*(concnt_old[1]-concnt_old[0])/(Dr*Dr)
    for i_ring in range(1,N_ring-1,1):
        concnt_model[i_ring] = concnt_old[i_ring] + Dt*Dif*( (r[i_ring]+0.5*Dr)*(concnt_old[i_ring+1]-concnt_old[i_ring]) + (r[i_ring]-0.5*Dr)*(concnt_old[i_ring-1]-concnt_old[i_ring]) )/(r[i_ring]*Dr*Dr)
    
    #compare model results with gaussian analytical results every 3600s (1 hour)
    if i_time%3600==0 :
        t = i_time
        Sigma_h = math.sqrt(Sigma_h0**2 + 2.0 * D_h * t)
        Sigma_v = math.sqrt(Sigma_v0**2 + 2.0 * D_v * t)
        
        for i_ring in range(0,N_ring-1,1):
            concnt_gaussian[i_ring] = Q / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i_ring]**2/Sigma_h**2 + 0.0**2/Sigma_v**2 ))
        
        
        # plot every 3600s (1hour) -----------------------------------
        plt.figure(figsize=(7,8))
        
        plt.plot(x, concnt_model[0:50], 'b-', x, concnt_gaussian[0:50], 'r--')
        
        plt.title( 'time = '+str(i_time/3600)+' hour' )
        
        plt.xlabel('distance (m)')
        plt.ylabel('concentration (kg/m3)')
        
        plt.legend(('model results', 'gaussian analytical results'),
           loc='upper right')
        
        plt.ylim((0.0, 55.0))
        #plt.yscale('log')
        
        #plt.yticks(y)
        
        plt.savefig(str(i_time/3600)+'_xy.png')
        plt.clf()
        plt.cla()



