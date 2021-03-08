from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt  

FILENAME='Lagrange_location_15day_10deg'
#FILENAME='Lagrange_location_15day_0deg'
a=np.loadtxt(FILENAME+'.txt')

z1 = a[0:2000:1,3]

t1 = np.arange(1,len(z1)+1,1)
print(z1)
print(t1)
#x=np.arange(20,350)

sValue = z1*0.0+4.0


# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# resolution = 'c' means use crude resolution coastlines.

plt.scatter(t1,z1,s=sValue,c='r',marker='.')

plt.title('Vertical location of a box calculated by Lagrange Module')
plt.xlabel('time')
plt.ylabel('pressure (hPa)')

#plt.xlim(left=-4,right=1)
#plt.ylim(top=2,bottom=-1)

#plt.legend()

plt.savefig(FILENAME+'_z.png')
