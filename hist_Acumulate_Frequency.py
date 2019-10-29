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
FILEDIR = '/n/home12/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES4x5_1mon/'

#NcFile = Dataset(FILEDIR+'Lagrange_Geos_Concentration_12deg.nc','r',format='NETCDF4_CLASSIC')
NcFile   = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4','r',format='NETCDF4_CLASSIC')

lat             = NcFile.variables['lat']
print(lat)
geos            = NcFile.variables['SpeciesConc_PASV1']
geos_Zsum       = np.sum(geos[:,:,:,:], axis=1)
geos_Zsum_Xmean = np.sum(geos_Zsum[:,:,:], axis=2)

lagrange            = NcFile.variables['LAGRG']
lagrange_Zsum       = np.sum(lagrange[:,:,:,:], axis=1)


print(geos_Zsum.shape)
print(lagrange_Zsum.shape)

Nt   = len(geos_Zsum[:,0,0])
Nlat = len(geos_Zsum[0,:,0])
Nlon = len(geos_Zsum[0,0,:])

geos_Zsum_Nt     = geos_Zsum[Nt-3,:,:]
lagrange_Zsum_Nt = lagrange_Zsum[Nt-3,:,:]
print(Nt)

geos_Zsum_1D     = geos_Zsum_Nt.reshape(Nlat*Nlon)
lagrange_Zsum_1D = lagrange_Zsum_Nt.reshape(Nlat*Nlon)


# Import data
x1 = geos_Zsum_1D
x2 = lagrange_Zsum_1D
print(x1.shape)

# Plot
plt.figure(figsize=(10,7), dpi= 80)

#plt.hist(x1, bins=50, range=(1e-10,7e-9))
plt.hist(x1, bins=50, range=(0,7e-9), density=True, histtype='step', cumulative=True,
        linewidth=2.5, label='Euler')
#plt.gca().set(title='Frequency Histogram', ylabel='Frequency')
#plt.ylim(0,500)
#plt.legend();

#plt.hist(x2, bins=50, range=(1e-10,7e-9))
plt.hist(x2, bins=50, range=(0,7e-9), density=True, histtype='step', cumulative=True, linewidth=1.0, label='Lagrange')
#plt.ylim(0,500)

plt.legend(loc='right')
plt.title('Cumulative step histograms')
plt.xlabel('Tracer Concentration (mol/mol)')
plt.ylabel('Frequency Number')

plt.savefig('Accumulate_Frequency_Distribution.png')
plt.clf()
plt.cla()
