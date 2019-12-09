import numpy as np
#from scipy.io import netcdf
#import netCDF4
from netCDF4 import Dataset


filename = 'GEOSChem.Restart.20150101_0000z.nc4'
f = Dataset(filename,'r+',format='NETCDF4_CLASSIC')

lats = f.variables['lat']
lons = f.variables['lon']
levs = f.variables['lev']
#OH = f.variables['SpeciesRst_OH']

#print(OH)
print('***')
print('*** ATTENTION: choose which lat and lon? ***')
print('***')
print(lats[35:45])
print(lons[1:38])
print(levs[38])

# quit()

# create dimensions
# time = f.createDimension('time',1) 
time = f.dimensions['time']
lev = f.dimensions['lev']
lat = f.dimensions['lat']
lon = f.dimensions['lon']

#define variables

PASV = 'SpeciesRst_PASV1'

PASV2 = f.createVariable(PASV,np.float32,('time','lev','lat','lon'),zlib=True,fill_value=-1e+31)

PASV2.long_name = 'passive species'
PASV2.units = 'mol mol-1'
PASV2.averaging_method = 'instantaneous'

f.variables[PASV][:,:,:,:]   = 0.000E-29
# f.variables[PASV][0,43,23,36]   = 0.1000E-04
# f.variables[PASV][0,38,15:30,7]   = 0.1000E-04
# f.variables[PASV][0,38,35:45,7]   = 0.1000E-04
# ( time, lev, lat, lon )

#close ncfile
f.close()
