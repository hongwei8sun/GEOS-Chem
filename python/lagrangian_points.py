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
# lagrange --------------------------------------
#------------------------------------------------
FILEDIR2 = '/n/home12/hongwei/hongwei/merra2_2x25_standard_10months/'

for i in range(1,29,1):
	lagrange_txt = np.loadtxt(FILEDIR2+'Lagrange_xyz_2015-1-'+str(i)+'-23:50:0.txt')
	
	print(len(lagrange_txt)) # 
	nbox = 6904224	 
	x1 = lagrange_txt[:,2]
	y1 = lagrange_txt[:,3]
	
	del lagrange_txt
	
	#------------------------------------------------
	# plot  -----------------------------------------
	#------------------------------------------------
	plt.figure(figsize=(7,9))
	
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
	

	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
	m.drawcoastlines()
	m.drawparallels(np.arange(-90.,91.,30.))
	m.drawmeridians(np.arange(-180.,181.,60.))
	m.drawmapboundary(fill_color='white')
	
	parallels = np.arange(-90.,90,30.)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
	meridians = np.arange(-180.,180.,60.)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
	
	sValue = x1*0.0+0.5
	plt.scatter(x1,y1,s=sValue,c='r',marker='.',zorder=10)
	
	plt.title('Lagrange (air parcels)', fontsize=10)
	
	plt.savefig(str(i)+'_xy.png')
	plt.clf()
	plt.cla()


