from netCDF4 import Dataset
import numpy as np

# Read UV wind speed data and lat/lon
NcFile  = Dataset('UVnc.nc4','r',format='NETCDF4_CLASSIC')

lon     = NcFile.variable[lon]
lat     = NcFile.variable[lat]
time    = NcFile.variable[time]
Dt = 60		# Unit: s

U       = NcFile.variable[U]	# U[times,lat,lon]
V       = NcFile.variable[V]


# define the initial location of particles
box_lat = [15,16]		# 15N
box_lon = [120,121]		# 120E


# integrate ODE
for it in range(len(time)):
    for i in range(len(box_lat)):
        curr_lon = box_lon[i]
        curr_lat = box_lat[i]

        uwnd = find_wind(U[it,:,:],curr_lon,curr_lat,lon,lat)	# find_wind() is a function calculating the wind speed at a specific location
        vwnd = find_wind(V[it,:,:],curr_lon,curr_lat,lon,lat)
        
        Dx = uwnd * Dt
        Dy = vwnd * Dt
        
        D_lon = m2lon(box_lat[i], box_lon[i], Dx)	# m2lon() is a function changes distance from m to corresponding lon degree
        D_lat = m2lat(box_lat[i], box_lon[i], Dy)	# m2lat() is a function change distance from m to corresponding lat degree
        
        box_lon[i] = box_lon[i] + D_lon
        box_lat[i] = box_lat[i] + D_lat

# save data box_lon, box_lat into a file, so that we can use this data to plot figures.


def find_wind(wind,it,curr_lat,curr_lon,lon,lat)
    i_lon = find_i(curr_lon,lon)	# find the nearest index of lon that less than curr_lon
    i_lat = find_i(curr_lat,lat)	# find the nearest index of lat that less than curr_lat
    
    for i in range(2):
        ii = i + i_lon
        if ii=IIPAR+1: ii = 0	# once the lon is bigger than 360 degree, change lon to 0 degree
        for j in range(2):
            jj   = j + i_lat
            num  = ii*2 + jj
            v[num] = wind[jj,ii]
            w[num] = IDW(curr_lon,curr_lat,lon[ii],lat[jj])	# IDW() is a function calculating the Inverse Distance Weight between 2 points
    
    Interplt_wind = w[:]*v[:]
    return Interplt_wind


def IDW(lon1,lat1,lon2,lat2)
    
    return w


def m2lon(lat, lon, dx)

    return D_lon


def m2lat(lat, lon, dy)

    return D_lat

