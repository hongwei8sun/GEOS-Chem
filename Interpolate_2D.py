from   mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
from   netCDF4 import Dataset as open_ncfile
import numpy as np

from scipy import interpolate

#-- open netcdf file
nc = open_ncfile('/n/home12/hongwei/hongwei/data/goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T3NVASM.5.12.4/2015/01/MERRA2_400.tavg3_3d_asm_Nv.20150101.nc4')

#-- read variable
var = nc.variables['U'][0,71,:,:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

x = np.arange(-5.01, 5.01, 0.25)
y = np.arange(-5.01, 5.01, 0.25)
xx, yy = np.meshgrid(x, y)
z = np.sin(xx**2+yy**2)
f = interpolate.interp2d(x, y, z, kind='cubic')

xnew = np.arange(-5.01, 5.01, 1e-2)
ynew = np.arange(-5.01, 5.01, 1e-2)
znew = f(xnew, ynew)
plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
plt.show()
