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
NA	= 6.022e23

FILEDIR1	= '/n/home12/hongwei/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES4x5_1yr/'

NcFile1	= Dataset(FILEDIR1+'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4','r',format='NETCDF4_CLASSIC')

lat1	= NcFile1.variables['lat']

# read concentration
geos1		= NcFile1.variables['SpeciesConc_PASV1']
lagrange1	= NcFile1.variables['LAGRG']

# read total air mass
AD_file1	= open(FILEDIR1+'State_Met_AD.txt','r')
GC1_AD	= geos1[0,:,:,:]*0.0

Nx	= len(geos1[0,0,0,:])
Ny	= len(geos1[0,0,:,0])
Nz	= len(geos1[0,:,0,0])

for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            line	= AD_file1.readline()
            GC1_AD[iz,iy,ix]	= float(line)

# read grid box area
AREA_file1	= open(FILEDIR1+'State_Met_AREA_M2.txt','r')
GC1_AREA	= np.zeros( (Nz, Ny, Nx), dtype=float )

for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            line	= AREA_file1.readline()
            GC1_AREA[iz,iy,ix]	= float(line)

GC1_AREA_Zmean	= np.mean(GC1_AREA[:,:,:], axis=0)

# change mol/mol to molec/cm2
geos11	= geos1[:,:,:,:]*0.0

for i in range(len(geos1[:,0,0,0])):
    geos11[i,:,:,:]	= geos1[i,:,:,:]*(GC1_AD[:,:,:]*1000.0/28.97)*NA         # [molec]

geos1_Zsum	= np.sum(geos11[:,:,:,:], axis=1)

for i in range(len(geos1[:,0,0,0])):
    geos1_Zsum[i,:,:]	= geos1_Zsum[i,:,:]/GC1_AREA[0,:,:]/1e4            # [molec/cm2]


#change mol/mol to molec/cm2
lagrange11	= lagrange1[:,:,:,:]*0.0

for i in range(len(lagrange1[:,0,0,0])):
    lagrange11[i,:,:,:]	= lagrange1[i,:,:,:]*(GC1_AD[:,:,:]*1000.0/28.97)*NA

lagrange1_Zsum	= np.sum(lagrange11[:,:,:,:], axis=1)

for i in range(len(lagrange1[:,0,0,0])):
    lagrange1_Zsum[i,:,:]	= lagrange1_Zsum[i,:,:]/GC1_AREA[0,:,:]/1e4


Nt1	= len(geos1_Zsum[:,0,0])-5
Nlat1	= len(geos1_Zsum[0,:,0])
Nlon1	= len(geos1_Zsum[0,0,:])

geos1_Zsum_Nt1	= geos1_Zsum[Nt1,:,:]
geos1_Zsum_1D	= geos1_Zsum_Nt1.reshape(Nlat1*Nlon1)


lagrange1_Zsum_Nt1	= lagrange1_Zsum[Nt1,:,:]
lagrange1_Zsum_1D	= lagrange1_Zsum_Nt1.reshape(Nlat1*Nlon1)

GC1_AREA_Zmean_1D	= GC1_AREA_Zmean.reshape(Nlat1*Nlon1)

geos1_1D     = np.vstack((geos1_Zsum_1D,GC1_AREA_Zmean_1D))
lagrange1_1D = np.vstack((lagrange1_Zsum_1D,GC1_AREA_Zmean_1D))


# Sort 2D numpy array by 1st row
geos1_2		= geos1_1D[ :, geos1_1D[0].argsort()]
lagrange1_2	= lagrange1_1D[ :, lagrange1_1D[0].argsort()]

# Import data
geos1_x1	= geos1_2[0,:]

Num_bin = 2000

hist,bins = np.histogram(geos1_x1,bins = Num_bin)
geos1_concent = bins
geos1_area = np.zeros(len(hist))

geos1_area[0] = sum( geos1_2[1,0:hist[0]-1] )
for i in range(len(hist)):
    geos1_area[i] = sum( geos1_2[ 1 , 0 : sum(hist[0:i+1]) ] )


lagrange1_x1	= lagrange1_2[0,:]

hist,bins	= np.histogram(lagrange1_x1,bins = Num_bin)
lagrange1_concent	= bins
lagrange1_area	= np.zeros(len(hist))

lagrange1_area[0]	= sum( lagrange1_2[1,0:hist[0]-1] )
for i in range(len(hist)):
    lagrange1_area[i]	= sum( lagrange1_2[ 1 , 0 : sum(hist[0:i+1]) ] )


# ====================================================================

geos1_Zsum       = np.sum(geos1[:,:,:,:], axis=1)

geos1_mass       = NcFile1.variables['Euler_mass']
geos1_mass_Zsum  = np.sum(geos1_mass[:,:,:,:], axis=1)

lagrange1_Zsum       = np.sum(lagrange1[:,:,:,:], axis=1)

lagrange1_mass       = NcFile1.variables['Lagra_mass']
lagrange1_mass_Zsum  = np.sum(lagrange1_mass[:,:,:,:], axis=1)


Nt1   = len(geos1_Zsum[:,0,0])-5
Nlat1 = len(geos1_Zsum[0,:,0])
Nlon1 = len(geos1_Zsum[0,0,:])

geos1_Zsum_Nt1          = geos1_Zsum[Nt1,:,:]
geos1_Zsum_1D          = geos1_Zsum_Nt1.reshape(Nlat1*Nlon1)

geos1_mass_Zsum_Nt1     = geos1_mass_Zsum[Nt1,:,:]
geos1_mass_Zsum_1D     = geos1_mass_Zsum_Nt1.reshape(Nlat1*Nlon1)

lagrange1_Zsum_Nt1      = lagrange1_Zsum[Nt1,:,:]
lagrange1_Zsum_1D      = lagrange1_Zsum_Nt1.reshape(Nlat1*Nlon1)

lagrange1_mass_Zsum_Nt1 = lagrange1_mass_Zsum[Nt1,:,:]
lagrange1_mass_Zsum_1D = lagrange1_mass_Zsum_Nt1.reshape(Nlat1*Nlon1)


geos1_1D     = np.vstack((geos1_Zsum_1D,geos1_mass_Zsum_1D))
lagrange1_1D = np.vstack((lagrange1_Zsum_1D,lagrange1_mass_Zsum_1D))


# Sort 2D numpy array by 1st row
geos1_2    = geos1_1D[ :, geos1_1D[0].argsort()]
lagrange1_2 = lagrange1_1D[ :, lagrange1_1D[0].argsort()]

# Import data
geos1_x1     = geos1_2[0,:]

hist,bins = np.histogram(geos1_x1,bins = Num_bin)
geos1_concent = bins
geos1_mass = np.zeros(len(hist))

geos1_mass[0] = sum( geos1_2[1,0:hist[0]-1] )
for i in range(len(hist)):
    geos1_mass[i] = sum( geos1_2[ 1 , 0 : sum(hist[0:i+1]) ] )


lagrange1_x1      = lagrange1_2[0,:]

hist,bins        = np.histogram(lagrange1_x1,bins = Num_bin)
lagrange1_concent = bins
lagrange1_mass    = np.zeros(len(hist))

lagrange1_mass[0] = sum( lagrange1_2[1,0:hist[0]-1] )
for i in range(len(hist)):
    lagrange1_mass[i] = sum( lagrange1_2[ 1 , 0 : sum(hist[0:i+1]) ] )


# =======================================================================
# Plot
plt.figure(figsize=(10,7), dpi= 80)

lagrange1_area = np.append(0,lagrange1_area)
lagrange1_mass = np.append(0,lagrange1_mass)
 
geos1_area = np.append(0,geos1_area)
geos1_mass = np.append(0,geos1_mass)

lagrange1_area = lagrange1_area/lagrange1_area[len(lagrange1_area)-1]
lagrange1_mass = lagrange1_mass/lagrange1_mass[len(lagrange1_mass)-1]
geos1_area = geos1_area/geos1_area[len(geos1_area)-1]
geos1_mass = geos1_mass/geos1_mass[len(geos1_mass)-1]

diag = np.linspace(0,1,100)

Al = np.trapz(diag, diag)
La = np.trapz(lagrange1_mass, lagrange1_area)
Eu = np.trapz(geos1_mass, geos1_area)

Gini_La = (Al-La)/Al
Gini_Eu = (Al-Eu)/Al

print(Al)
print('Gini_La_Eu')
print(Gini_La)
print(Gini_Eu)

# plot = = = =
plt.plot(lagrange1_area[0:len(lagrange1_area)],lagrange1_mass[0:len(lagrange1_mass)],'b',label='Lagrange 4*5')
plt.plot(geos1_area[0:len(geos1_area)],geos1_mass[0:len(geos1_mass)],'r',label='Euler 4*5')

plt.plot(diag,diag,'k.')


#plt.legend(loc='upper left')
plt.title(' Cumulative Fraction from Lower to Higher Concentration')
plt.xlabel(' Cumulate Area Fraction ')
plt.ylabel(' Cumulate Mass Fraction ')

plt.savefig('Cumulate_Area_Mass_Frequency_Distribution.png')
plt.clf()
plt.cla()
