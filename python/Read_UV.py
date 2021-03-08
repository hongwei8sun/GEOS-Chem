import numpy as np
#from scipy.io import netcdf
#import netCDF4
from netCDF4 import Dataset


filename1 = '/n/home12/hongwei/hongwei/data/goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T3NVASM.5.12.4/2015/01/MERRA2_400.tavg3_3d_asm_Nv.20150101.nc4'
filename2 = '/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/data/ExtData/GEOS_4x5/GEOS_FP/2015/01/GEOSFP.20150101.A3dyn.4x5.nc'

f1 = Dataset(filename1,'r',format='NETCDF4_CLASSIC')
f2 = Dataset(filename2,'r',format='NETCDF4_CLASSIC')

# lats = f.variables['lat']
# lons = f.variables['lon']
# levs = f.variables['lev']

uwnd1 = f1.variables['U']
uwnd2 = f2.variables['U']

#print(OH)
print(uwnd1[2,70,1,:])
print(uwnd2[2,70,1,:])
#  [ 8 <time> x 72 <lev> x 46 <lat> x 72 <lon> ]

# create dimensions
# time = f.createDimension('time',1) 

#time = f.dimensions['time']

#define variables


#PASV2 = f.createVariable(PASV,np.float32,('time','lev','lat','lon'),zlib=True,fill_value=-1e+31)

#PASV2.long_name = 'passive species'
#PASV2.units = 'mol mol-1'
#PASV2.averaging_method = 'instantaneous'

#f.variables[PASV][:,:,:,:]   = 0.1000E-29
# f.variables[PASV][0,43,23,36]   = 0.1000E-04
# f.variables[PASV][0,38,15:30,7]   = 0.1000E-04
#f.variables[PASV][0,38,35:45,7]   = 0.1000E-04
# ( time, lev, lat, lon )

#close ncfile
f.close()

