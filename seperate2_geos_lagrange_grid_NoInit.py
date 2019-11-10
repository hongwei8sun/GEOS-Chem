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
            GC_AREA[iz,iy,ix] = float(line)

###
geos2           = geos[:,:,:,:]*0.0
for i in range(len(geos[:,0,0,0])):
    geos2[i,:,:,:] = geos[i,:,:,:]*GC_AD[:,:,:]

geos_Zmean       = np.mean(geos2[:,:,:,:], axis=1)
for i in range(len(geos[:,0,0,0])):
    geos_Zmean[i,:,:] = geos_Zmean[i,:,:]/GC_AREA[0,:,:]/1e4		# [mol/cm2]

geos_Zmean_Xmean = np.mean(geos_Zmean[:,:,:], axis=2)

lagrange2           = lagrange[:,:,:,:]*0.0
for i in range(len(lagrange[:,0,0,0])):
    lagrange2[i,:,:,:] = lagrange[i,:,:,:]*GC_AD[:,:,:]

lagrange_Zmean       = np.mean(lagrange2[:,:,:,:], axis=1)
for i in range(len(lagrange[:,0,0,0])):
    lagrange_Zmean[i,:,:]   = lagrange_Zmean[i,:,:]/GC_AREA[0,:,:]/1e4*28.97/98.0


lagrange_Zmean_Xmean = np.mean(lagrange_Zmean[:,:,:], axis=2)

# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(geos_Zmean[:,0,0])
print(Nt)

i=359
while i<Nt:
	plt.figure(figsize=(14,8))
	
	plt.subplot(2, 2, 3)
	#ax = fig.add_axes([0.1,0.1,0.4,0.4])
	
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
	m.drawcoastlines()
	m.drawparallels(np.arange(-90.,91.,30.))
	m.drawmeridians(np.arange(-180.,181.,60.))
	m.drawmapboundary(fill_color='white')
	
	parallels = np.arange(-90.,90,30.)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
	meridians = np.arange(-180.,180.,60.)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
	#m.drawcountries()
	
	# for geos ==================================================================
	ny = geos_Zmean.shape[1]; nx = geos_Zmean.shape[2]
	lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
	x, y = m(lons, lats) # compute map proj coordinates.
	print(i)
	
	bounds = np.linspace(0.0E-14,30.0E-14,11)
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	cs = m.pcolormesh(x, y, geos_Zmean[i,:,:], norm=norm, cmap='Reds')
	cs.cmap.set_under('w')
	cs.set_clim(1.0E-20)

	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
	cbar.set_label('mol/cm2')
	
	plt.title('GEOS-Chem (Concentration)', fontsize=10)
	#plt.title('GEOS-Chem (Blue/Shaded) & Lagrange (Red/Scatter)')

	plt.suptitle('Day: '+str(i+1), fontsize=16)
	
	# for GOES distribution ===================================================
	plt.subplot(2,2,4)
	plt.plot(geos_Zmean_Xmean[i,:],lat)
#	plt.xscale('log')
	plt.xlim(5.0E-14,40.0E-14)
	plt.ylim(-90,90)
	plt.xlabel('Tracer Concentration (mol/cm2)')
	plt.ylabel('Latitude (deg)')

	# for Lagrange ============================================================
	plt.subplot(2, 2, 1)

	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
	m.drawcoastlines()
	m.drawparallels(np.arange(-90.,91.,30.))
	m.drawmeridians(np.arange(-180.,181.,60.))
	m.drawmapboundary(fill_color='white')
	
	parallels = np.arange(-90.,90,30.)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
	meridians = np.arange(-180.,180.,60.)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
	
	bounds = bounds
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	cs = m.pcolormesh(x, y, lagrange_Zmean[i,:,:], norm=norm, cmap='Reds')
	cs.cmap.set_under('w')
	cs.set_clim(1.0E-20)
	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
	cbar.set_label('mol/cm2')
	
	plt.title('Lagrange (air parcels percent)', fontsize=10)
		

	# for Lagrange Distribution: ==============================================
	plt.subplot(2,2,2)
	plt.plot(lagrange_Zmean_Xmean[i,:],lat)

#	plt.xscale('log')
	plt.xlim(5.0E-14,40.0E-14)
	plt.ylim(-90,90)
	plt.ylabel('Latitude (deg)')


	plt.savefig(str(i+1)+'_xy2.png')
	plt.clf()
	plt.cla()
	
	i = i + 1
