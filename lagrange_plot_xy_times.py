from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt  
import math



#print(x1)
#x=np.arange(20,350)

#--------------------------------------------------------------
# sValue = x1*0.0+7.0


# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180)


m.drawcoastlines()
#m.fillcontinents(color='white',lake_color='white')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,91.,30.))
m.drawmeridians(np.arange(-180.,181.,60.))
m.drawmapboundary(fill_color='white')

parallels = np.arange(-90.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-180.,180.,60.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

#m.scatter(x1,y1,s=sValue,c='r',marker='.')
#l1=plt.plot(x1,y1,'r--',label='type1')
#plt.plot(x1,y1,s=sValue,c='r',marker='.')

#---------------------------------------------------------------

FILENAME='Lagrange_box_i_lon_lat_lev'
#FILENAME='Lagrange_location_15day_0deg'
a=np.loadtxt(FILENAME+'.txt')

print(len(a))
ntimes = math.floor(len(a)/1000.0)
x1 = np.arange(ntimes*1000).reshape(ntimes, 1000)
y1 = np.arange(ntimes*1000).reshape(ntimes, 1000)


sValue = x1[0,:]*0.0+4.0
i = 0
Ndt = 50

while i < (ntimes-1):
	ii = i*Ndt                             # plot in every 10 time steps
	x1[i,:] = a[ii*1000:(ii+1)*1000:1,1]  # there are 1000 not 1001 in a[i*1000:(i+1)*1000:1,1] 
	y1[i,:] = a[ii*1000:(ii+1)*1000:1,2]
	plt.scatter(x1[i,:],y1[i,:],s=sValue,c='r',marker='.',zorder=10)
	plt.title('Time(s): '+ str(ii))
	plt.savefig(str(i)+'_xy.png')
	i=i+1
