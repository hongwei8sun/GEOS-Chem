from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


nc = Dataset('GEOSChem.SpeciesConc_inst.20160701_0000z.nc4','r',format='NETCDF4_CLASSIC')

pasv1 = nc.variables['SpeciesConc_PASV1']
time = nc.variables['time']
print(pasv1)

# create figure and axes instances
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

ny = pasv1.shape[2]; nx = pasv1.shape[3]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

ntimes = len(pasv1[:,0,0,0])
print(ntimes)
i=0
while i<ntimes:
	print(i)
	#clevs = [0.1E-09,0.1e-08,0.1E-06,0.1E-05,0.1E-04,0.1E-03]
	cs = m.contourf(x,y,pasv1[i,43,:,:])
	i=i+1
	# add colorbar.
	cbar = m.colorbar(cs,location='bottom',pad="5%")
	cbar.set_label('mol mol-1')
	# add title
	plt.title('GEOS-Chem tracer transportation')
	plt.savefig(str(i)+'_xy.png')
plt.show()


