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

pasv_Ymean = np.sum(pasv1[:,:,:,:], axis=2)

lon = geos_nc.variables['lon']
Pcenter = np.loadtxt("P_73_72_levels.txt")
Pz = Pcenter[1:145:2]
lev = Pz[::-1]
print(Pz[:])
print(Pz.shape)

del geos_nc
print(pasv_Ymean.shape)
del pasv1
#------------------------------------------------
# lagrange --------------------------------------
#------------------------------------------------

#lagrange_txt=np.loadtxt(FILEDIR+'Lagrange_box_i_lon_lat_lev.txt')
lagrange_txt=np.loadtxt(FILEDIR+'Lagrange_1hr_box_i_lon_lat_lev.txt')

print(len(lagrange_txt)) # 1 hour
nbox = 80000
ntimes = math.floor(len(lagrange_txt)/nbox)


x1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)
z1 = np.arange( ntimes * nbox ).reshape(ntimes, nbox)

i = 0
Ndt = 24  # 24 hour
ii = i*Ndt                             # plot in every 10 time steps

while ii < (ntimes-1):
	print(i)
	x1[i,:] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 1]  # there are 1000 not 1001 in a[i*1000:(i+1)*1000:1,1] 
	z1[i,:] = lagrange_txt[ii*nbox : (ii+1)*nbox : 1, 3]
	i=i+1
	ii = i*Ndt                             # plot in every 10 time steps

del lagrange_txt
#------------------------------------------------
# plot  -----------------------------------------
#------------------------------------------------

# time step for ploting is 24 hours (once every day)
Nt = len(pasv_Ymean[:,0,0])
print(Nt)
sValue = x1[0,:]*0.0+0.5

i=0
while i<Nt:
	print(i)

	plt.set_yscale("log")
	bounds = np.array([1.0E-11,5.0E-11,1.0E-10,5.0E-10,1.0E-09,5.0E-09,1.0E-08,5.0E-08,1.0E-07,5.0E-07,1.0E-06])
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
	fig1 = plt.pcolormesh(lon[:], lev[:], pasv_Ymean[i,:,:], norm=norm, cmap="Reds")
	plt.invert_yaxis()

	fig1.cmap.set_under('w')
	fig1.set_clim(1.0E-11)
	fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0))
	fig.colorbar(fig1, ax=plt, orientation='vertical', format=fmt)
	fig1.set_label('mol mol-1')

	
	# for Lagrange *****
	i=i+1   # because there is no origin data in GEOS file
	fig2 = plt.scatter(x1[i,:],z1[i,:],s=sValue,c='r',marker='.',zorder=10)
	# add title

	plt.set_title('GEOS-Chem (Blue/Shaded) & Lagrange (Red/Scatter)')

	plt.suptitle('Day: '+str(i), fontsize=16)

	plt.subplots_adjust(hspace=0.5)

	plt.savefig(str(i)+'_xy.png')
	plt.clf()
	plt.cla()
