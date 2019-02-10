#from mpl_toolkits.basemap import Basemap, cm
#from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math

#DILEDIR = '/Users/hongweisun/Desktop/EPS-237 Planetary Radiation and Climate/Problem Set/Set 2'
f_txt=np.loadtxt('REMSpressure.txt')

sol = f_txt[:,0]
ps  = f_txt[:,3]

#print(sol)
plt.figure(figsize=(9,5.5))

#plt.figure()
plt.subplot(211)
#plt.plot(sol1, ps1, 'C3', lw=3)
plt.scatter(sol, ps, s=3)
plt.title('Surface pressure vs. Time')
plt.xlabel('Sol (entire sampling period)')
plt.ylabel('PS1')

# Scatter plot on top of lines
plt.subplot(212)
#plt.plot(sol2, ps2, 'C3', zorder=1, lw=3)
plt.scatter(sol[0:300], ps[0:300], s=3, zorder=2)
plt.tight_layout()
plt.xlabel('Sol (a few martain sols)')
plt.ylabel('PS1')

plt.savefig('set2_3.png')
#plt.axis([0, 6, 0, 20])
#plt.show()
