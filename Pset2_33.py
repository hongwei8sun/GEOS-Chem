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
sol_new = np.linspace(10, 1258, 12481)
#sol_new = np.linspace(sol[0], sol[len(sol)-1], num=len(sol), endpoint=True)
ps_new  = f(sol_new)
print(sol_new)

# Fourier transform
ps_fft = abs(np.fft.fft(ps_new))
xf = np.arange(len(ps_new))        # Frequency

#print(sol)

plt.figure(figsize=(9,5.5))

#plt.figure()
plt.subplot(211)
#plt.plot(sol1, ps1, 'C3', lw=3)
plt.scatter(sol_new, ps_new, s=3)
plt.title('Surface pressure vs. Time')
plt.xlabel('Sol (entire sampling period)')
plt.ylabel('PS1')

# Scatter plot on top of lines
plt.subplot(212)
#plt.plot(sol2, ps2, 'C3', zorder=1, lw=3)
plt.plot(xf[3:], ps_fft[3:])
#plt.tight_layout()
#plt.xlabel('Sol (a few martain sols)')
#plt.ylabel('PS1')

plt.savefig('set2_33.png')

