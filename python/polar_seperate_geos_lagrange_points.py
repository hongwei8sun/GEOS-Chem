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
FILEDIR = '/n/home12/hongwei/GC_lagrange/rundirs/geosfp_4x5_gc_timing/'
nbox = 360*4

geos_nc = Dataset(FILEDIR+'GEOSChem.SpeciesConc_inst.20160701_0000z.nc4','r',format='NETCDF4_CLASSIC')

pasv1 = geos_nc.variables['SpeciesConc_PASV1']
print(pasv1)

pasv_mean = np.sum(pasv1[:,:,:,:], axis=1)    # sum for the whole vertical levels

del geos_nc
print(pasv_mean.shape)
del pasv1
#------------------------------------------------
# lagrange --------------------------------------
#------------------------------------------------

lagrange_txt=np.loadtxt(FILEDIR+'Lagrange_1day_box_i_lon_lat_lev.txt')

print(len(lagrange_txt)) # 
# nbox = 6000
ntimes = math.floor(len(lagrange_txt)/nbox)

x1 = np.zeros( (ntimes, nbox) )
y1 = np.zeros( (ntimes, nbox) )

#x1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)
#y1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)

i = 0
Ndt = 1  #
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
	meridians = np.arange(-180.,-180.,60.)
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
	
	plt.suptitle('Day: '+str(i+1), fontsize=16)

	# for Lagrange ===================
	plt.subplot(2, 1, 1)
	
	# setup north polar stereographic basemap.
	# The longitude lon_0 is at 6-o'clock, and the
	# latitude circle boundinglat is tangent to the edge
	# of the map at lon_0. Default value of lat_ts
	# (latitude of true scale) is pole.
	m = Basemap(projection='nplaea',boundinglat=70,lon_0=270,resolution='l')
	m.drawcoastlines()
	# m.fillcontinents(color='coral',lake_color='aqua')
	# draw parallels and meridians.
	m.drawparallels(np.arange(70.,90.,10.))
	m.drawmeridians(np.arange(-180.,181.,20.))
	# m.drawmapboundary(fill_color='aqua')

	print(m.xmax,m.ymax)

	i = i + 1
	# draw tissot's indicatrix to show distortion.
	ax = plt.gca()
	lon1, lat1 = m(x1[i,:],y1[i,:])

	sValue = x1[0,:]*0.0+0.1
	plt.scatter(lon1,lat1,sValue,c='r',marker='.',zorder=10)
	
	plt.title('Lagrange (air parcels)', fontsize=10)
	
	plt.savefig(str(i)+'_xy.png')
	plt.clf()
	plt.cla()


# initial status:

geos0_nc = Dataset(FILEDIR+'GEOSChem.Restart.20160701_0000z.nc4','r',format='NETCDF4_CLASSIC')

pasv1 = geos0_nc.variables['SpeciesRst_PASV1']
print(pasv1)

pasv0_mean = np.sum(pasv1[:,:,:,:], axis=1)    # sum for the whole vertical levels

del geos0_nc
print(pasv0_mean.shape)
del pasv1


plt.figure(figsize=(7,9))

plt.subplot(2, 1, 2)
#ax = fig.add_axes([0.1,0.1,0.4,0.4])

m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
m.drawcoastlines()
m.drawparallels(np.arange(-2.,91.,4.))
m.drawmeridians(np.arange(-150.,181.,10.))
m.drawmapboundary(fill_color='white')

parallels = np.arange(-90.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
meridians = np.arange(-180.,-180.,60.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
#m.drawcountries()

# for geos =====================
ny = pasv0_mean.shape[1]; nx = pasv0_mean.shape[2]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
# for GEOS **********
#clevs = [0.1E-10,0.1E-06,0.5E-06,0.1E-05,0.5E-05,0.1E-04]
#cmap=plt.cm.get_cmap('Blues', 7)

# uneven bounds changes the colormapping:
bounds = np.array([1.0E-11,5.0E-11,1.0E-10,5.0E-10,1.0E-09,5.0E-09,1.0E-08,5.0E-08,1.0E-07,5.0E-07,1.0E-06])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

cs = m.pcolormesh(x, y, pasv0_mean[0,:,:], norm=norm, cmap='Blues')
cs.cmap.set_under('w')
cs.set_clim(1.0E-11)

# add colorbar.
fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0, 0))

cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
cbar.set_label('mol mol-1')

plt.title('GEOS-Chem (Concentration)', fontsize=10)
#plt.title('GEOS-Chem (Blue/Shaded) & Lagrange (Red/Scatter)')

plt.suptitle('Day: 0', fontsize=16)


# for Lagrange ===================
plt.subplot(2, 1, 1)

# setup north polar stereographic basemap.
# The longitude lon_0 is at 6-o'clock, and the
# latitude circle boundinglat is tangent to the edge
# of the map at lon_0. Default value of lat_ts
# (latitude of true scale) is pole.
m = Basemap(projection='nplaea',boundinglat=70,lon_0=270,resolution='l')
m.drawcoastlines()
# m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,80.,10.))
m.drawmeridians(np.arange(-180.,181.,20.))
# m.drawmapboundary(fill_color='aqua')

print(m.xmax,m.ymax)

# draw tissot's indicatrix to show distortion.
ax = plt.gca()
lon1, lat1 = m(x1[0,:],y1[0,:])

sValue = x1[0,:]*0.0+0.8
plt.scatter(lon1,lat1,sValue,c='r',marker='.',zorder=10)

plt.title('Lagrange (air parcels)', fontsize=10)

plt.savefig('0_xy.png')
plt.clf()
plt.cla()

