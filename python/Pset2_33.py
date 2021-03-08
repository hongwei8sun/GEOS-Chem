#from mpl_toolkits.basemap import Basemap, cm
#from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math
from scipy.interpolate import interp1d
from scipy.fftpack import fft,ifft

#DILEDIR = '/Users/hongweisun/Desktop/EPS-237 Planetary Radiation and Climate/Problem Set/Set 2'
f_txt=np.loadtxt('REMSpressure.txt')

sol = f_txt[:,0]
ps  = f_txt[:,3]

# Linuear interpolation
f = interp1d(sol, ps)

# the new data's time intervel is every 
#sol_new = np.linspace(10, 1258, 12481)

sol_new = np.linspace(10, 1258, 124801)
ps_new  = f(sol_new)

xx = np.arange(len(sol_new))

# Fourier transform
ps_fft = np.fft.fft(ps_new)
f = abs(ps_fft)/len(ps_new)
print('ps_fft',f[1247:1249])
#print(sol)

plt.figure(figsize=(9,11))


# Scatter plot on top of lines
plt.subplot(311)
#plt.plot(sol2, ps2, 'C3', zorder=1, lw=3)
plt.plot(xx[1:], f[1:])
plt.yscale('log')
plt.title('Fourier transform')
#plt.tight_layout()
#plt.xlabel('Sol (a few martain sols)')
#plt.ylabel('PS1')

plt.subplot(312)
#plt.plot(sol2, ps2, 'C3', zorder=1, lw=3)
plt.plot(xx[1:2000], f[1:2000])
plt.scatter(xx[1248], f[1248],facecolor='red',zorder=2)
plt.title('Diurnal Frequencies (red point)')

plt.subplot(313)
#plt.plot(sol2, ps2, 'C3', zorder=1, lw=3)
plt.plot(xx[1:30], f[1:30])
plt.scatter(xx[4], f[4],facecolor='red',zorder=2)
plt.title('Annual Frequencies (red point)')

plt.savefig('set2_33.png')

