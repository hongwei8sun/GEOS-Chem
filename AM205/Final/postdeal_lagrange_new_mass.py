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
nbox     = 6904224

nbox_day = 18864
N_boxes  = nbox_day

# for function find_lon(x, Xmin, Dx):
Xmin = -182.5
Dx   = 5.0

# For function find_lat(y, Ymin, Dy):
Ymin = -92.0
Dy   = 4.0


# input file:
# files location/directory:
FILEDIR = '/n/home12/hongwei/hongwei/merra2_4x5_standard_1year/'

#------------------------------------------------
# geos ------------------------------------------
#------------------------------------------------

# GEOS-Chem file:
geos_nc = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4','r',format='NETCDF4_CLASSIC')

pasv1 = geos_nc.variables['SpeciesConc_PASV1']
lon   = geos_nc.variables['lon']
lat   = geos_nc.variables['lat']
lev   = geos_nc.variables['lev']

Nx    = len(lon)
Ny    = len(lat)
Nz    = len(lev)
Nt    = len(geos_nc.variables['time'])
print('pasv1',pasv1)

del geos_nc


# For function find_lev(z, Pedge)
Press = np.loadtxt("P_73_72_levels.txt")
Pedge2 = Press[0:145:2]
Pedge = Pedge2[::-1]
print('Pedge.shape',Pedge.shape)

news = pasv1[:,:,:,:]
news[:,:,:,:] = 0.0
print('news.shape',news.shape)

Lagra_mass = pasv1[:,:,:,:]
Lagra_mass[:,:,:,:] = 0.0

Euler_mass = pasv1[:,:,:,:]
Euler_mass[:,:,:,:] = 0.0

# GC_AD = pasv1[0,:,:,:] # lev, lat, lon
# GC_AD = 0.0


#------------------------------------------------
# total air mass in each grid  ------------------
#------------------------------------------------

AD_file = open(FILEDIR+'State_Met_AD.txt','r')

GC_AD = np.zeros( (Nz, Ny, Nx), dtype=float )

print(Nx)
print(Ny)
print(Nz)
for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            line = AD_file.readline()
            GC_AD[iz,iy,ix] = float(line)

for it in range(Nt):
    Euler_mass[it,:,:,:] = pasv1[it,:,:,:] * 98.0/28.97 * GC_AD[:,:,:]

del pasv1
#------------------------------------------------
# lagrange --------------------------------------
#------------------------------------------------
# Lagrange file:
# lagrange_txt=np.loadtxt(FILEDIR+'Lagrange_1day_box_i_lon_lat_lev.txt')

input_file = open(FILEDIR+'Lagrange_1day_box_i_lon_lat_lev.txt','r')


lagr = np.zeros( (Nt,nbox,4), dtype=float )

for itime in range(Nt):
    print(itime)
    for ibox in range(nbox):
        lines = input_file.readline()
        data = lines.split() #split string into a list
        lagr[itime,ibox,0] = float(data[1])
        lagr[itime,ibox,1] = float(data[2])
        lagr[itime,ibox,2] = float(data[3])

print('len(lagr)',len(lagr)) # 1 hour


#ntimes = math.floor(len(lagrange_txt)/nbox)

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
	if delt_P[index_min]>=0:
		i_z = index_min
	if delt_P[index_min]<0:
		i_z = index_min - 1
	return i_z;
#------------------------------------------------

t = 0
print(Nt)
while t < (Nt-2):   # 30 number in total 
# minus 1 for python begin at 0, minus 1 for lagragne is 1 more time than GEOS 
	print('t',t)     # 0~29 -> 1~30
	n = 0
	tt = t + 1
	while n < N_boxes:
		i_lon = int( (lagr[tt,n,0] - Xmin) / Dx )   # find_lon(lagr[t,n,0], Xmin, Dx)
		if i_lon==72:
			i_lon = 0
		i_lat = int( (lagr[tt,n,1] - Ymin) / Dy )   # find_lat(lagr[t,n,1], Ymin, Dy)
		i_lev = find_lev(lagr[tt,n,2] , Pedge )
		news[t,i_lev,i_lat,i_lon]       = news[t,i_lev,i_lat,i_lon] + 110.0/GC_AD[i_lev,i_lat,i_lon]*(28.97/98.0)  	# unit mol/mol

		Lagra_mass[t,i_lev,i_lat,i_lon] = Lagra_mass[t,i_lev,i_lat,i_lon] + 110.0  	# unit [kg]
		
		
		n = n + 1
	print('n',n)
	t = t + 1
	N_boxes = N_boxes + nbox_day

# write news variable into nc file ----------------
filename = 'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4'
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

LAGRG2.long_name = 'Lagrange parcels concentration'
LAGRG2.units = 'mol/mol'
LAGRG2.averaging_method = 'instantaneous'

f.variables[LAGRG][:,:,:,:]   = news[:,:,:,:]

# for the mass variable:
LAGR  = 'Lagra_mass'
LAGR2 = f.createVariable(LAGR,np.float32,('time','lev','lat','lon'),zlib=True,fill_value=-1e+31)

LAGR2.long_name = 'Lagrange parcels concentration in mass'
LAGR2.units = 'kg'
LAGR2.averaging_method = 'instantaneous'

f.variables[LAGR][:,:,:,:]   = Lagra_mass[:,:,:,:]

#

Euler = 'Euler_mass'
Euler2 = f.createVariable(Euler,np.float32,('time','lev','lat','lon'),zlib=True,fill_value=-1e+31)

Euler2.long_name = 'Lagrange parcels concentration in mass'
Euler2.units = 'kg'
Euler2.averaging_method = 'instantaneous'

f.variables[Euler][:,:,:,:]   = Euler_mass[:,:,:,:]

#close ncfile
f.close()


