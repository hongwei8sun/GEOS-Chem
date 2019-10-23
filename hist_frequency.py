import seaborn as sns
sns.set_style("white")
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math
#pandas
#------------------------------------------------
#------------------------------------------------
FILEDIR = '/n/home12/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES4x5/'

#NcFile = Dataset(FILEDIR+'Lagrange_Geos_Concentration_12deg.nc','r',format='NETCDF4_CLASSIC')
NcFile   = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4','r',format='NETCDF4_CLASSIC')

geos      = NcFile.variables['SpeciesConc_PASV1']
geos_mean = np.sum(geos[:,:,:,:], axis=1)

lagrange        = NcFile.variables['LAGRG']
lagrange_mean   = np.sum(lagrange[:,:,:,:], axis=1)

Nt   = len(geos_mean[:,0,0])
Nlat = len(geos_mean[0,:,0])
Nlon = len(geos_mean[0,0,:])

print(geos_mean.shape)

geos_mean_Nt     = geos_mean[Nt-2,:,:]
lagrange_mean_Nt = lagrange_mean[Nt-2,:,:]

geos_mean_1D     = geos_mean_Nt.reshape(Nlat*Nlon)
lagrange_mean_1D = geos_mean_Nt.reshape(Nlat*Nlon)


# Import data
x1 = geos_mean_1D
x2 = lagrange_mean_1D
print(x1.shape)

# Plot
plt.figure(figsize=(10,7), dpi= 80)

plt.subplot(2, 1, 1)
plt.hist(x1, bins=50,range=(1e-9,0.8e-7))
plt.gca().set(title='Frequency Histogram', ylabel='Frequency');
#plt.xlim(50,75)
#plt.legend();

plt.subplot(2, 1, 2)
plt.hist(x1, bins=50,range=(1e-9,0.8e-7))

plt.savefig('Frequency_Distribution.png')
plt.clf()
plt.cla()
