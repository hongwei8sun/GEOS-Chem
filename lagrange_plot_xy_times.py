from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt  
import math


#---------------------------------------------------------------

FILENAME='Lagrange_0deg_box_i_lon_lat_lev'
#FILENAME='Lagrange_10deg_box_i_lon_lat_lev'
# FILENAME='Lagrange_location_15day_10deg'
a=np.loadtxt(FILENAME+'.txt')

print(len(a))
nbox = 80000
ntimes = math.floor(len(a)/nbox)

x1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)
y1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)


sValue = x1[0,:]*0.0+4.0
i = 0
Ndt = 50

while i < (ntimes-1):
	print(i)
	
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
	m.drawcoastlines()
	m.drawparallels(np.arange(-90.,91.,30.))
	m.drawmeridians(np.arange(-180.,181.,60.))
	m.drawmapboundary(fill_color='white')
	parallels = np.arange(-90.,90,30.)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
	meridians = np.arange(-180.,180.,60.)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
	
	ii = i*Ndt                             # plot in every 10 time steps
	x1[i,:] = a[ii*nbox : (ii+1)*nbox : 1, 1]  # there are 1000 not 1001 in a[i*1000:(i+1)*1000:1,1] 
	y1[i,:] = a[ii*nbox : (ii+1)*nbox : 1, 2]
	plt.scatter(x1[i,:],y1[i,:],s=sValue,c='r',marker='.',zorder=10)
	plt.title(str(ii))
	plt.savefig(str(i)+'_xy.png')
	plt.clf()
	plt.cla()
	i=i+1
