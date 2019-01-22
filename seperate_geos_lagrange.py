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
# geos ------------------------------------------
#------------------------------------------------
#FILEDIR = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_standard_tracer/'
FILEDIR = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_standard_12deg/'

geos_nc = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20160701_0000z.nc4','r',format='NETCDF4_CLASSIC')

pasv1 = geos_nc.variables['SpeciesConc_PASV1']
print(pasv1)

pasv_mean = np.sum(pasv1[:,:,:,:], axis=1)

del geos_nc
print(pasv_mean.shape)
del pasv1
#------------------------------------------------
# lagrange --------------------------------------
#------------------------------------------------

lagrange_txt=np.loadtxt(FILEDIR+'Lagrange_1hr_box_i_lon_lat_lev.txt')

print(len(lagrange_txt)) # 1 hour
nbox = 80000
ntimes = math.floor(len(lagrange_txt)/nbox)


x1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)
y1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)

i = 0
Ndt = 24  # 24 hour
ii = i*Ndt                             # plot in every 10 time steps

while ii < (ntimes-1):
	print(i)
	x1[i,:] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 1]  # there are 1000 not 1001 in a[i*1000:(i+1)*1000:1,1] 
	y1[i,:] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 2]
	i=i+1
	ii = i*Ndt                             # plot in every 10 time steps

del lagrange_txt
#------------------------------------------------
# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(pasv_mean[:,0,0])
print(Nt)
sValue = x1[0,:]*0.0+0.5

i=0
while i<Nt:
	plt.figure(figsize=(7,9))

	plt.subplot(2, 1, 2)
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
	
	# for geos =====================
	ny = pasv_mean.shape[1]; nx = pasv_mean.shape[2]
	lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
	x, y = m(lons, lats) # compute map proj coordinates.
	# for GEOS **********
	print(i)
	#clevs = [0.1E-10,0.1E-06,0.5E-06,0.1E-05,0.5E-05,0.1E-04]
	#cmap=plt.cm.get_cmap('Blues', 7)
	
	# uneven bounds changes the colormapping:
	bounds = np.array([1.0E-11,5.0E-11,1.0E-10,5.0E-10,1.0E-09,5.0E-09,1.0E-08,5.0E-08,1.0E-07,5.0E-07,1.0E-06])
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	#pcm = ax[1].pcolormesh(X, Y, Z, norm=norm, cmap='RdBu_r')
	#cs = m.contourf(x, y, pasv_mean[i,43,:,:], norm=norm, cmap='Blues')
	#pasv_ave = pasv1.sum(axis=1)
	cs = m.pcolormesh(x, y, pasv_mean[i,:,:], norm=norm, cmap='Blues')
	cs.cmap.set_under('w')
	cs.set_clim(1.0E-11)
	# add colorbar.
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	
	cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
	cbar.set_label('mol mol-1')
	
	plt.title('GEOS-Chem (Concentration)', fontsize=10)
	#plt.title('GEOS-Chem (Blue/Shaded) & Lagrange (Red/Scatter)')
	
	plt.suptitle('Day: '+str(i), fontsize=16)

	# for Lagrange ===================
	plt.subplot(2, 1, 1)
	
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
	m.drawcoastlines()
	m.drawparallels(np.arange(-90.,91.,30.))
	m.drawmeridians(np.arange(-180.,181.,60.))
	m.drawmapboundary(fill_color='white')
	
	parallels = np.arange(-90.,90,30.)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
	meridians = np.arange(-180.,180.,60.)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
	
	
	i=i+1   # because there is no origin data in GEOS file
	plt.scatter(x1[i,:],y1[i,:],s=sValue,c='r',marker='.',zorder=10)
	
	plt.title('Lagrange (air parcels)', fontsize=10)
	
	
	plt.savefig(str(i)+'_xy.png')
	plt.clf()
	plt.cla()
