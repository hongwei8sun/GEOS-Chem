from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
#pandas
#------------------------------------------------
# geos ------------------------------------------
#------------------------------------------------

geos_nc = Dataset('./xy_times_geos/GEOSChem.SpeciesConc_inst.20160701_0000z.nc4','r',format='NETCDF4_CLASSIC')

pasv1 = geos_nc.variables['SpeciesConc_PASV1']
time = geos_nc.variables['time']
print(pasv1)

#------------------------------------------------
# lagrange --------------------------------------
#------------------------------------------------

lagrange_txt=np.loadtxt('./xy_times_lagrange/Lagrange_0deg_box_i_lon_lat_lev.txt')

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

#------------------------------------------------
# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(pasv1[:,0,0,0])
print(Nt)
sValue = x1[0,:]*0.0+0.5

i=0
while i<Nt:
	
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_axes([0.1,0.1,0.8,0.8])

	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
	m.drawcoastlines()
	m.drawparallels(np.arange(-90.,91.,30.))
	m.drawmeridians(np.arange(-180.,181.,60.))
	m.drawmapboundary(fill_color='white')
	
	parallels = np.arange(-90.,90,30.)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
	meridians = np.arange(-180.,180.,60.)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
	#m.drawcountries()
	
	# for geos =====================
	ny = pasv1.shape[2]; nx = pasv1.shape[3]
	lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
	x, y = m(lons, lats) # compute map proj coordinates.
	# for GEOS **********
	print(i)
	#clevs = [0.1E-10,0.1E-06,0.5E-06,0.1E-05,0.5E-05,0.1E-04]
	#cmap=plt.cm.get_cmap('Blues', 7)
	cs = m.contourf(x,y,pasv1[i,43,:,:],norm=colors.LogNorm(vmin=pasv1[:,43,:,:].min(), vmax=pasv1[:,43,:,:].max()),cmap='Blues')
	print(cs.cmap)
	cs.cmap.set_under('w')
	cs.set_clim(0.1E-10)
	i=i+1
	# add colorbar.
	cbar = m.colorbar(cs,location='bottom',pad="5%")
	cbar.set_label('mol mol-1')
	
	# for Lagrange *****
	plt.scatter(x1[i,:],y1[i,:],s=sValue,c='r',marker='.',zorder=10)
	# add title
	plt.title('GEOS-Chem (Blue/Shaded) & Lagrange (Red/Scatter)')
	plt.savefig(str(i)+'_xy.png')
	plt.clf()
	plt.cla()
