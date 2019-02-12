#from mpl_toolkits.basemap import Basemap, cm
#from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math

#DILEDIR = '/Users/hongweisun/Desktop/EPS-237 Planetary Radiation and Climate/Problem Set/Set 2'

P1   = np.arange(32.3,2500.0, 10)
T_P1 = 520.0 * (P1/2500.0)**0.286
P2   = np.arange(10.0,32.3, 1)
T_P2 = P2 * 0.0 + 150.0

Z1   = np.arange(0, 75, 1)
print(Z1)
T_Z1 = 520.0 - 5.0 * Z1
Z2   = np.arange(74, 84, 1)
T_Z2 = Z2 * 0.0 + 150.0

#print(sol)
plt.figure(figsize=(9,8))

#plt.figure()
plt.subplot(121)
plt.plot(T_P1, P1)
plt.plot(T_P2, P2)
plt.yscale('log')
#plt.title('Temperature vs. Pressure')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')

plt.gca().invert_yaxis()


# Scatter plot on top of lines
plt.subplot(122)
plt.plot(T_Z1, Z1)
plt.plot(T_Z2, Z2)
#plt.title('Temperature vs. Pressure')
#plt.tight_layout()
plt.xlabel('Temperature (K)')
plt.ylabel('Height (km)')

plt.savefig('set2_2.png')
#plt.axis([0, 6, 0, 20])
#plt.show()
