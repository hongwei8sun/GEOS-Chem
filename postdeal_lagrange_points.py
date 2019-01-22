from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math
#pandas

#------------------------------------------------------
# some input parameters:
#------------------------------------------------------

# input file:
# files location/directory:
FILEDIR = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_standard_tracer/'
# GEOS-Chem file:
geos_nc = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20160701_0000z.nc4','r',format='NETCDF4_CLASSIC')
# Lagrange file:
lagrange_txt=np.loadtxt(FILEDIR+'Lagrange_box_i_lon_lat_lev.txt')

# modify detailed setting in their own section

# for function find_lon(x, Xmin, Dx):
Xmin = -182.5
Dx = 5.0

# For function find_lat(y, Ymin, Dy):
Ymin = -92.0
Dy = 4.0

# For function find_lev(z, Pedge)
Press = np.loadtxt("Pedge_73levels.txt")
Pedge = Press[0:145:2]


#------------------------------------------------
# geos ------------------------------------------
#------------------------------------------------

pasv1 = geos_nc.variables['SpeciesConc_PASV1']
lon = geos_nc.variables['lon']
lat = geos_nc.variables['lat']
print(pasv1)

del geos_nc

news = pasv1[:,:,:,:]
news[:,:,:,:] = 0.0
print(news.shape)

del pasv1

#------------------------------------------------
# lagrange --------------------------------------
#------------------------------------------------

print(len(lagrange_txt)) # 1 hour
nbox = 80000
ntimes = math.floor(len(lagrange_txt)/nbox)
print(ntimes)

Ndt = 24  # 24 hour
Nt = math.floor(ntimes/Ndt)
print(Nt)
lagr = np.arange( Nt * nbox * 3 ).reshape(Nt, nbox, 3)

i = 0
ii = i*Ndt                             # plot in every 10 time steps

while i < Nt:
	print(i)
	lagr[i,:,0] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 1]  # there are 1000 not 1001 in a[i*1000:(i+1)*1000:1,1] 
	lagr[i,:,1] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 2]
	lagr[i,:,2] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 3]
	i=i+1
	ii = i*Ndt                             # plot in every 10 time steps
print(lagr.shape)
del lagrange_txt
#------------------------------------------------
# post deal to regrid data  ---------------------
#------------------------------------------------

#def find_lon( x, Xmin, Dx ):
#	i_x = int( (x - Xmin) / Dx )
#	return i_x;

#def find_lat( y, Ymin, Dy ):
#	i_y = int( (y - Ymin) / Dy )
#	return i_y;

def find_lev( z, Pedge ):
	delt_P = Pedge - z
	index_min = np.argmin(abs(delt_P))
	if delt_P[index_min]<=0:
		i_z = index_min
	if delt_P[index_min]>0:
		i_z = index_min - 1
	return i_z;
#------------------------------------------------

t = 0
n = 0
print(Nt)
while t < Nt:
	print(t)
	while n < nbox:
		i_lon = int( (lagr[t,n,0] - Xmin) / Dx )   # find_lon(lagr[t,n,0], Xmin, Dx)
		i_lat = int( (lagr[t,n,1] - Ymin) / Dy )   # find_lat(lagr[t,n,1], Ymin, Dy)
		i_lev = find_lev(lagr[t,n,2] , Pedge )
		
		news[t,i_lev,i_lon,i_lat] = news[t,i_lev,i_lon,i_lat] + 1.0/nbox
		n = n + 1
	t = t + 1

# write news variable into nc file ----------------
filename = 'GEOSChem.SpeciesConc_inst.20160701_0000z.nc4'
f = Dataset(filename,'r+',format='NETCDF4_CLASSIC')

# create dimensions
# time = f.createDimension('time',1) 
time = f.dimensions['time']
lev = f.dimensions['lev']
lat = f.dimensions['lat']
lon = f.dimensions['lon']

#define variables

LAGRG = 'LAGRG'

LAGRG2 = f.createVariable(LAGRG,np.float32,('time','lev','lat','lon'),zlib=True,fill_value=-1e+31)

LAGRG2.long_name = 'Lagrange parcels concentration (100% in total)'
LAGRG2.units = '100%'
LAGRG2.averaging_method = 'instantaneous'

f.variables[LAGRG][:,:,:,:]   = news[:,:,:,:]

#close ncfile
f.close()