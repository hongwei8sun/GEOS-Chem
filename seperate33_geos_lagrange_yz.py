from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math
#pandas

NA = 6.022e+23
#------------------------------------------------
FILEDIR2 = '/n/home12/hongwei/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES2x25_Feb/'
NcFile   = Dataset(FILEDIR2+'GEOSChem.SpeciesConc_inst.20150201_0000z.nc4','r',format='NETCDF4_CLASSIC')

lat             = NcFile.variables['lat']
lon             = NcFile.variables['lon']
geos            = NcFile.variables['SpeciesConc_PASV1']
lagrange        = NcFile.variables['LAGRG']

#------------------------------------------------
FILEDIR = '/n/home12/hongwei/hongwei/merra2_2x25_standard_Feb/'

# grid height
H = np.loadtxt("Height_73_72_levels_km.txt")	# [km]

H2 = H[0:145:2]
Dh2 = np.zeros(len(H2)-1)

for i in range(len(H2)-1):
    Dh2[i] = H2[i]-H2[i+1]

Dh = Dh2[::-1]*1000*100         # km -> cm

Height = H[1:145:2]
z = Height[::-1]

###

geos_Xmean = np.mean(geos[:,:,:,:],axis=3)

lagrange_Xmean = np.mean(lagrange[:,:,:,:],axis=3)


# plot  -----------------------------------------
iz = 50
y  = lat
z  = z[0:iz]

Y, Z = np.meshgrid(y, z)

i = 26

plt.figure(figsize=(16,8))

# lagrange
plt.subplot(121)
levels = np.linspace(1e-11,120e-11,101)
plt.contourf(Y, Z, lagrange_Xmean[i,0:iz,:], levels, cmap='Reds')
plt.colorbar();


plt.xlabel('Lat [deg]')
plt.ylabel('Height [km]')
plt.title('Lagrange: Tracer concentration in Day='+str(i+1)+' [mol/mol]', fontsize=10)

# geos
plt.subplot(122)
plt.contourf(Y, Z, geos_Xmean[i,0:iz,:], levels, cmap='Reds')
plt.colorbar();


plt.xlabel('Lat [deg]')
plt.ylabel('Height [km]')
plt.title('Euler: Tracer concentration in Day='+str(i+1)+' [mol/mol]', fontsize=10)

plt.savefig(str(i+1)+'_mol_yz.png')
plt.clf()
plt.cla()

