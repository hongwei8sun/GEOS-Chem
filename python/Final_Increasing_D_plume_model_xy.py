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

Sigma_h0 = 50.0  			# [m]
Sigma_v0 = 50.0
Dr = 5.0  				#[m]

D_h = 1.0  				# diffusion coefficient, [m2/s]
D_v = 1.0

Q = 100.0* PI*Sigma_h0*Sigma_v0 	# [kg/m3 * m2 = kg/m]

x = [0.5*Dr]
z = [0.0]
concnt_gaussian = [0.0]			# concentration from gaussian analytical solution

N_ring = 50

for i_ring in range(N_ring-1):
    x.append( (i_ring+1.5)*Dr )
    z.append( 0.0 )
    concnt_gaussian.append( 0.0 )
print(x)

R = [Dr]
S = [PI*R[0]**2]
for i_ring in range(N_ring-1):
    R.append( (i_ring+2)*Dr )
    S.append( PI*R[i_ring+1]**2 - PI*R[i_ring]**2 )


Sigma_h = Sigma_h0
Sigma_v = Sigma_v0

# initial concentration distribution from gaussion analytical results---------
for i_ring in range(N_ring-1):
        concnt_gaussian[i_ring] = Q / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i_ring]**2/Sigma_h**2 + z[i_ring]**2/Sigma_v**2 ))


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

Dt = 1 					#[s]
for i_time in range(1,200,Dt): 		# [s]
    Dif = Dif + 5
    D_h = D_h + 5
    D_v = D_v + 5

    CFL = Dr - 2.0*Dif/Dr*Dt
    print('CFL',CFL)


    # Combine 2 adjacent rings into 1 ring:
    if (CFL<0) :

        #update the concentration for new rings
        for i_ring in range(int(N_ring/2)-1):
            concnt_model[i_ring] = (concnt_model[2*i_ring]*S[2*i_ring]+concnt_model[2*i_ring+1]*S[2*i_ring+1]) / (S[2*i_ring]+S[2*i_ring+1])
            print('S1+S2',S[2*i_ring]+S[2*i_ring+1])
        for i_ring in range(int(N_ring/2)-1,N_ring-1,1):
            concnt_model[i_ring] = 0.0

        # update parameters:
        del x, R, S, r
        Dr = Dr*2

        x = [0.5*Dr]
        r = [0.5*Dr]
        R = [Dr]
        S = [PI*R[0]**2]
        for i_ring in range(N_ring-1):
            x.append( (i_ring+1.5)*Dr )
            r.append( (i_ring+1.5)*Dr )
            R.append( (i_ring+2)*Dr )
            S.append( PI*R[i_ring+1]**2 - PI*R[i_ring]**2 )

        total = []
        for i in range(0, len(S)):
            total.append(concnt_model[i] * S[i])
        print('amount 000: ', sum(total))
        print('S',S)


    # calculate new concentration after dilusion
    for i_ring in range(0,N_ring-1,1):
        concnt_old[i_ring] = concnt_model[i_ring]

    concnt_model[0] = concnt_old[0] + Dt*Dif*2.0*(concnt_old[1]-concnt_old[0])/(Dr*Dr)
    for i_ring in range(1,N_ring-1,1):
        concnt_model[i_ring] = concnt_old[i_ring] + Dt*Dif*( (r[i_ring]+0.5*Dr)*(concnt_old[i_ring+1]-concnt_old[i_ring]) + (r[i_ring]-0.5*Dr)*(concnt_old[i_ring-1]-concnt_old[i_ring]) )/(r[i_ring]*Dr*Dr)

    total = [] 
    for i in range(0, len(S)): 
        total.append(concnt_model[i] * S[i])
    print('amount1: ', sum(total))

    #compare model results with gaussian analytical results every 3600s (1 hour)
    Sigma_h = math.sqrt(Sigma_h**2 + 2.0 * D_h * Dt)
    Sigma_v = math.sqrt(Sigma_v**2 + 2.0 * D_v * Dt)
    if (i_time-1)%1==0 :
        t = i_time
        print('t', t, ';     Dr', Dr)
        
        for i_ring in range(0,N_ring-1,1):
            concnt_gaussian[i_ring] = Q / (2.0*PI*Sigma_h*Sigma_v) * math.exp(-0.5*( x[i_ring]**2/Sigma_h**2 + z[i_ring]**2/Sigma_v**2 ))
        
        # plot every 3600s (1hour) -----------------------------------
        plt.figure(figsize=(5,8))
        
        plt.plot(x, concnt_model, 'b-', x, concnt_gaussian, 'r--')
        
        plt.title( 'time = '+str(int(i_time/1))+' second' )
        
        plt.xlabel('distance (m)')
        plt.ylabel('concentration (kg/m3)')
        
        plt.legend(('model results', 'gaussian analytical results'),
           loc='upper right')
        
        plt.ylim((0.0, 55.0))
        plt.xlim((0.0, Dr*N_ring))
        #plt.yscale('log')
        
        #plt.yticks(y)
        
        plt.savefig(str(int(i_time/1))+'_xy.png')
        plt.clf()
        plt.cla()



