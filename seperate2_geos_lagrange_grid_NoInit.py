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
FILEDIR = '/n/home12/hongwei/GC_Python/Postdeal_lagrange_points/LAGR_GOES4x5_1mon/'

#NcFile = Dataset(FILEDIR+'Lagrange_Geos_Concentration_12deg.nc','r',format='NETCDF4_CLASSIC')
NcFile   = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20150101_0000z.nc4','r',format='NETCDF4_CLASSIC')

lat             = NcFile.variables['lat']
print(lat)
geos            = NcFile.variables['SpeciesConc_PASV1']
geos_Zsum       = np.sum(geos[:,:,:,:], axis=1)
geos_Zsum_Xmean = np.sum(geos_Zsum[:,:,:], axis=2)

lagrange            = NcFile.variables['LAGRG']
lagrange_Zsum       = np.sum(lagrange[:,:,:,:], axis=1)
lagrange_Zsum_Xmean = np.sum(lagrange_Zsum[:,:,:], axis=2)

# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(geos_Zsum[:,0,0])
print(Nt)

i=0
while i<Nt:
	plt.figure(figsize=(10,8))
	
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
	# for GEOS **********
	print(i)
	#clevs = [0.1E-10,0.1E-06,0.5E-06,0.1E-05,0.5E-05,0.1E-04]
	#cmap=plt.cm.get_cmap('Blues', 7)
	
	# uneven bounds changes the colormapping:
#	bounds = np.array([1.0E-12,1.0E-11,1.0E-10,1.0E-09,1.0E-08,1.0E-07])
#	bounds = np.array([5.0E-11,10.0E-11,15.0E-11,20.0E-11,25.0E-11,30.0E-11])
	bounds = np.linspace(1.0E-12,1.0E-06,10)
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	#pcm = ax[1].pcolormesh(X, Y, Z, norm=norm, cmap='RdBu_r')
	#cs = m.contourf(x, y, pasv_Zsum[i,43,:,:], norm=norm, cmap='Blues')
	#pasv_ave = pasv1.sum(axis=1)
	cs = m.pcolormesh(x, y, geos_Zsum[i,:,:], norm=norm, cmap='Reds')
	cs.cmap.set_under('w')
	cs.set_clim(1.0E-20)
	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
	cbar.set_label('mol/mol')
	
	plt.title('GEOS-Chem (Concentration)', fontsize=10)
	#plt.title('GEOS-Chem (Blue/Shaded) & Lagrange (Red/Scatter)')

	plt.suptitle('Day: '+str(i+1), fontsize=16)
	
	# for GOES distribution ===================================================
	plt.subplot(2,2,4)
	plt.plot(geos_Zsum_Xmean[i,:],lat)
#	plt.xscale('log')


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
	
#	bounds = np.array([5.0E-06,1.0E-05,5.0E-05,1.0E-04,5.0E-04,1.0E-03,5.0E-03,1.0E-02,5.0E-02,1.0E-01,5.0E-01])
#	bounds = np.array([1.0E-12,5.0E-12,1.0E-11,5.0E-11,1.0E-10,5.0E-10,1.0E-09,5.0E-09,1.0E-08,5.0E-08,1.0E-07,5.0E-07,1.0E-06])
#	bounds = np.array([1.0E-12,1.0E-11,1.0E-10,1.0E-09,1.0E-08,1.0E-07])
	bounds = np.linspace(1.0E-12,1.0E-06,10)
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	cs = m.pcolormesh(x, y, lagrange_Zsum[i,:,:], norm=norm, cmap='Reds')
	cs.cmap.set_under('w')
	cs.set_clim(1.0E-20)
	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
	cbar.set_label('mol/mol')
	
	plt.title('Lagrange (air parcels percent)', fontsize=10)
		
	# for Lagrange Distribution: ==============================================
	plt.subplot(2,2,2)
	plt.plot(lagrange_Zsum_Xmean[i,:],lat)
#	plt.xscale('log')

	plt.savefig(str(i+1)+'_xy2.png')
	plt.clf()
	plt.cla()
	
	i = i + 1
