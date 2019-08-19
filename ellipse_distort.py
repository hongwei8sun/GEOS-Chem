from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math

#from basic_units import cm
#from matplotlib.patches import Rectangle, Ellipse
from matplotlib import patches

#
xcenter, ycenter = 0.0, 0.0
width, height = 4, 2
angle = -30

theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
x = 0.9 * width * np.cos(theta)
y = 0.9 * height * np.sin(theta)

rtheta = np.radians(angle)
R = np.array([
    [np.cos(rtheta), -np.sin(rtheta)],
    [np.sin(rtheta),  np.cos(rtheta)],
    ])


x, y = np.dot(R, np.array([x, y]))
x += xcenter
y += ycenter

#
fig = plt.figure(figsize=(30, 10))

ax = fig.add_subplot(111)

ax.fill(x, y, alpha=0.5, facecolor='white',
        edgecolor='white', linewidth=1, zorder=1)

e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=0, linewidth=2, fill=False, zorder=2)

ax.add_patch(e1)


e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=45, linewidth=2, fill=False, zorder=2)

ax.add_patch(e1)


e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=90, linewidth=2, fill=False, zorder=2)

ax.add_patch(e1)

plt.ylim(-5, 5)
plt.xlim(-15, 15)

fig.savefig('ellipse_compare')
