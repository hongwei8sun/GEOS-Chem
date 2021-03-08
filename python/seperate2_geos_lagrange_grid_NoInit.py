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
FILEDIR2 = '/n/home12/hongwei/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES4x5_1yr/'
#NcFile = Dataset(FILEDIR2+'Lagrange_Geos_Concentration_12deg.nc','r',format='NETCDF4_CLASSIC')
NcFile   = Dataset(FILEDIR2+'GEOSChem.SpeciesConc_inst.20150201_0000z.nc4','r',format='NETCDF4_CLASSIC')

NA = 6.022e+23

lat             = NcFile.variables['lat']
geos            = NcFile.variables['SpeciesConc_PASV1']
lagrange        = NcFile.variables['LAGRG']

#------------------------------------------------
FILEDIR = '/n/home12/hongwei/hongwei/merra2_4x5_standard_1year/'

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
    geos2[i,:,:,:] = geos[i,:,:,:]*(GC_AD[:,:,:]*1000.0/28.97)*NA	# [mol/mol] to [molec]

geos_Zsum       = np.sum(geos2[:,:,:,:], axis=1)
for i in range(len(geos[:,0,0,0])):
    geos_Zsum[i,:,:] = geos_Zsum[i,:,:]/GC_AREA[0,:,:]/1e4		# [molec/cm2]

geos_Zsum_Xmean = np.mean(geos_Zsum[:,:,:], axis=2)


lagrange2           = lagrange[:,:,:,:]*0.0
for i in range(len(lagrange[:,0,0,0])):
    lagrange2[i,:,:,:] = lagrange[i,:,:,:]*(GC_AD[:,:,:]*1000.0/28.97)*NA

lagrange_Zsum       = np.sum(lagrange2[:,:,:,:], axis=1)
for i in range(len(lagrange[:,0,0,0])):
    lagrange_Zsum[i,:,:]   = lagrange_Zsum[i,:,:]/GC_AREA[0,:,:]/1e4

lagrange_Zsum_Xmean = np.mean(lagrange_Zsum[:,:,:], axis=2)

# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(geos_Zsum[:,0,0])
print(Nt)

i=359
while i<Nt:
	total_geos   = sum(sum(geos_Zsum[i,:,:]))
	percent_geos  = geos_Zsum[i,:,:]/total_geos
	total_lagr   = sum(sum(lagrange_Zsum[i,:,:]))
	percent_lagr = lagrange_Zsum[i,:,:]/total_lagr
	print(total_geos-total_lagr)
	
	percent_geos_Xmean = geos_Zsum_Xmean[i,:]/sum(geos_Zsum_Xmean[i,:])
	percent_lagr_Xmean = lagrange_Zsum_Xmean[i,:]/sum(lagrange_Zsum_Xmean[i,:])
	
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
	ny = geos_Zsum.shape[1]; nx = geos_Zsum.shape[2]
	lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
	x, y = m(lons, lats) # compute map proj coordinates.
	print(i)
	
	bounds = np.linspace(0.0001E-3,1.5E-3,31)
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	cs = m.pcolormesh(x, y, percent_geos, norm=norm, cmap='bwr')
	cs.cmap.set_under('w')
	cs.set_clim(bounds[0])

	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
	cbar.set_label('Unit: (molec/cm2)/molec')
	
	plt.title('Eulerian', fontsize=10)

	plt.suptitle('Day: '+str(i+1), fontsize=16)
	
	
	# for GOES distribution ===================================================
	plt.subplot(2,2,4)
	plt.plot(percent_geos_Xmean[:],lat)

	X_max = 0.06
	plt.xlim(0,X_max)
	plt.ylim(-90,90)
	plt.xlabel('Tracer Concentration Percentage [(molec/cm2)/molec]')
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
	
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	cs = m.pcolormesh(x, y, percent_lagr, norm=norm, cmap='bwr')
	cs.cmap.set_under('w')
	cs.set_clim(bounds[0])
	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
#	cbar.set_label('(molec/cm2)/molec')
	
	plt.title('Lagrangian', fontsize=10)
		

	# for Lagrange Distribution: ==============================================
	plt.subplot(2,2,2)
	plt.plot(percent_lagr_Xmean,lat)

#	plt.xscale('log')
	plt.xlim(0,X_max)
	plt.ylim(-90,90)
	plt.ylabel('Latitude (deg)')


	plt.savefig(str(i+1)+'_xy2.png')
	plt.clf()
	plt.cla()
	
	i = i + 1
