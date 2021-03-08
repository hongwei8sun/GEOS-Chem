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
FILEDIR = '/n/home12/hongwei/Python_script/Postdeal_lagrange_points/'

# NcFile = Dataset(FILEDIR+'Lagrange_Geos_Concentration_0deg.nc','r',format='NETCDF4_CLASSIC')
NcFile = Dataset(FILEDIR+'Lagrange_Geos_Concentration_12deg.nc','r',format='NETCDF4_CLASSIC')

geos = NcFile.variables['SpeciesConc_PASV1']
geos_Ymean = np.sum(geos[:,:,:,:], axis=3)
print(geos)

lagrange = NcFile.variables['LAGRG']
lagrange_Ymean = np.sum(lagrange[:,:,:,:], axis=3)
print(lagrange)


lat = NcFile.variables['lat']
Pcenter = np.loadtxt("P_73_72_levels.txt")
Pz = Pcenter[1:145:2]
lev = Pz[::-1]
#Pz = np.linspace(1,72,72)
print(lat.shape)
print(geos_Ymean.shape)
print(lagrange_Ymean.shape)

#------------------------------------------------
# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(geos_Ymean[:,0,0])
print(Nt)

i=0
while i<Nt:
	print(i)
	fig, (ax1, ax2) = plt.subplots(nrows=2)

	# GEOS-Chem #
	ax2.set_yscale("log") 
	bounds = np.array([1.0E-11,5.0E-11,1.0E-10,5.0E-10,1.0E-09,5.0E-09,1.0E-08,5.0E-08,1.0E-07,5.0E-07,1.0E-06])
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	fig2=ax2.pcolormesh(lat[:], lev[:], geos_Ymean[i,:,:], norm=norm, cmap="Blues")
	ax2.invert_yaxis()
	
	fig2.cmap.set_under('w')
	fig2.set_clim(1.0E-11)
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	fig.colorbar(fig2, ax=ax2, orientation='vertical', format=fmt)
	fig2.set_label('mol mol-1')
	
	ax2.set_title('GEOS-Chem')
	
	# lagrange #
	ax1.set_yscale("log")
	bounds = np.array([1.0E-06,5.0E-06,1.0E-05,5.0E-05,1.0E-04,5.0E-04,1.0E-03,5.0E-03,1.0E-02,5.0E-02,1.0E-01])
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	fig1=ax1.pcolormesh(lat[:], lev[:], lagrange_Ymean[i,:,:], norm=norm, cmap="Reds")
	ax1.invert_yaxis()
	
	fig1.cmap.set_under('w')
	fig1.set_clim(1.0E-11)
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	fig.colorbar(fig1, ax=ax1, orientation='vertical', format=fmt)
	fig1.set_label('100%')
	
	ax1.set_title('Lagrange')
	
	
	plt.subplots_adjust(hspace=0.5)	
	
	plt.suptitle('Day: '+str(i), fontsize=16)
	plt.savefig(str(i)+'_xy2.png')
	plt.clf()
	plt.cla()
	
	i = i + 1
