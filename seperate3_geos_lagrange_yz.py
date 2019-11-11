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
FILEDIR2 = '/n/home12/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES4x5_1yr/'
#NcFile = Dataset(FILEDIR2+'Lagrange_Geos_Concentration_12deg.nc','r',format='NETCDF4_CLASSIC')
NcFile   = Dataset(FILEDIR2+'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4','r',format='NETCDF4_CLASSIC')

lat             = NcFile.variables['lat']
lon             = NcFile.variables['lon']
geos            = NcFile.variables['SpeciesConc_PASV1']
lagrange        = NcFile.variables['LAGRG']

#------------------------------------------------
FILEDIR = '/n/home12/hongwei/hongwei/merra2_4x5_standard/'

#------------------------------------------------
# total air mass in each grid  ------------------
#------------------------------------------------

AD_file = open(FILEDIR+'State_Met_AD.txt','r')

GC_AD = geos[0,:,:,:]*0.0

Nx = len(geos[0,0,0,:])
Ny = len(geos[0,0,:,0])
Nz = len(geos[0,:,0,0])

for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            line = AD_file.readline()
            GC_AD[iz,iy,ix] = float(line)


AREA_file = open(FILEDIR+'State_Met_AREA_M2.txt','r')

GC_AREA = np.zeros( (Nz, Ny, Nx), dtype=float )

for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            line = AREA_file.readline()
            GC_AREA[iz,iy,ix] = float(line)*1.0e2	# m -> cm

# grid height
H = np.loadtxt("Height_73_72_levels_km.txt")

H2 = H[0:145:2]
Dh2 = np.zeros(len(H2)-1)

for i in range(len(H2)-1):
    Dh2[i] = H2[i]-H2[i+1]

Dh = Dh2[::-1]*1000*100		# km -> cm

Height = H[1:145:2]
z = Height[::-1]

# grid volume
GC_V = GC_AREA[:,:,:] * 0.0

for i in range(Nx):
    for j in range(Ny):
        GC_V[:,j,i] = GC_AREA[:,j,i]*Dh[:]	# [cm3]


###
geos2           = geos[:,:,:,:]*0.0
for i in range(len(geos[:,0,0,0])):
    geos2[i,:,:,:] = geos[i,:,:,:]*GC_AD[:,:,:]/GC_V[:,:,:]	# [mol/cm3]

geos2_Xmean = np.mean(geos2[:,:,:,:],axis=3)


lagrange2           = lagrange[:,:,:,:]*0.0
for i in range(len(lagrange[:,0,0,0])):
    lagrange2[i,:,:,:] = lagrange[i,:,:,:]*GC_AD[:,:,:]/GC_V[:,:,:]	# [mol/cm3]

lagrange2_Xmean = np.mean(lagrange2[:,:,:,:],axis=3)

print(z[37:40])
print(lagrange2_Xmean[359,30:50,:])

# plot  -----------------------------------------
y = lat
z = z

Y, Z = np.meshgrid(y, z)

i = 200

plt.figure(figsize=(16,8))

# lagrange
plt.subplot(121)
levels = np.linspace(0.1e-15,3.1e-15,11)
plt.contourf(Y, Z, lagrange2_Xmean[i,:,:], levels, cmap='Reds')
plt.colorbar();


plt.xlabel('Lat [deg]')
plt.ylabel('Height [km]')
plt.title('Lagrange: Tracer concentration in Day='+str(i+1)+' [mol/cm3]', fontsize=10)

# geos
plt.subplot(122)
plt.contourf(Y, Z, geos2_Xmean[i,:,:], levels, cmap='Reds')
plt.colorbar();


plt.xlabel('Lat [deg]')
plt.ylabel('Height [km]')
plt.title('Euler: Tracer concentration in Day='+str(i+1)+' [mol/cm3]', fontsize=10)

plt.savefig(str(i+1)+'_yz.png')
plt.clf()
plt.cla()

