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
widths  = [200, 269.69, 413.43, 578.44, 750.83, 926.48]
heights = [200, 148.32,  96.75, 69.151, 53.275, 43.174]

widths2   = [400, 539.38, 826.87, 1159.9, 1501.7, 1853.0]
heights2  = [400, 296.64, 193.50,  138.3, 106.55, 86.348]

angles   = [0.0, 0.73535, 1.0659, 1.2178, 1.3012, 1.3532]

width  = 200
height = 200
angle  = 0.0

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
for i in range(6):
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111)

    ax.fill(x, y, alpha=0.5, facecolor='white',
            edgecolor='white', linewidth=1, zorder=1)

    e1 = patches.Ellipse((xcenter, ycenter), widths[i], heights[i],
                angle=90-angles[i]/3.14*180.0, linewidth=2, fill=False, zorder=2)
    
    e2 = patches.Ellipse((xcenter, ycenter), widths2[i], heights2[i],
                angle=90-angles[i]/3.14*180.0, linewidth=2, fill=False, zorder=2)
    
    print(widths[i])
    
    ax.add_patch(e1)
    ax.add_patch(e2)
    
    plt.ylim(-500, 500)
    plt.xlim(-1000, 1000)
    
    fig.savefig(str(i)+'_ellipse')
