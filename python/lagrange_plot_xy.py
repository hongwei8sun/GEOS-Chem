from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt  

FILENAME='Lagrange_location_15day_10deg'
#FILENAME='Lagrange_location_15day_0deg'
a=np.loadtxt(FILENAME+'.txt')

x1=a[1:2000:20,1]
y1=a[1:2000:20,2]

print(x1)
#x=np.arange(20,350)

sValue = x1*0.0+7.0


# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180)


m.scatter(x1,y1,s=sValue,c='r',marker='.')


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

m.scatter(x1,y1,s=sValue,c='r',marker='.')
#l1=plt.plot(x1,y1,'r--',label='type1')
#plt.plot(x1,y1,s=sValue,c='r',marker='.')

plt.title('Location of a box calculated by Lagrange Module')
#plt.xlabel('lon')
#plt.ylabel('lat')

#plt.xlim(left=-4,right=1)
#plt.ylim(top=2,bottom=-1)

#plt.legend()

plt.savefig(FILENAME+'_xy.png')
