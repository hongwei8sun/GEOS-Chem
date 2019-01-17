import numpy as np
#from scipy.io import netcdf
#import netCDF4
from netCDF4 import Dataset


filename = 'GEOSChem.Restart.20160701_0000z.nc4'
f = Dataset(filename,'r+',format='NETCDF4_CLASSIC')

#lat = f.variables['lat']
# lon = f.variables['lon']
# lev = f.variables['lev']
OH = f.variables['SpeciesRst_OH']

print(OH)
#print(lat[23])
# print(lon[36])

# create dimensions
# time = f.createDimension('time',1) 
time = f.dimensions['time']
lev = f.dimensions['lev']
lat = f.dimensions['lat']
lon = f.dimensions['lon']

#define variables

PASV = 'PASV1'

PASV2 = f.createVariable(PASV,np.float64,('time','lev','lat','lon'),zlib=True,fill_value=-1e+31)

PASV2.long_name = 'passive species'
PASV2.units = 'mol mol-1'
PASV2.averaging_method = 'instantaneous'

f.variables[PASV][:,:,:,:]   = 0.1000E-29
f.variables[PASV][0,43,23,36]   = 0.1000E-04
t = f.dimensions['time']
print(t)
#close ncfile
f.close()

