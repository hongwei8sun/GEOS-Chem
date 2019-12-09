from netCDF4 import Dataset
import numpy as np
import math

Re = 6.3781e+6 	# radius of earth, unit: [m]

#=================================================
# Read UV wind speed data and lat/lon
NcFile  = Dataset('./data/UV_wind.nc','r')

lon     = NcFile.variables['lon']
lat     = NcFile.variables['lat']

IIPAR	= len(lon)-1
JJPAR	= len(lat)-1

time    = NcFile.variables['time']
Dt	= 600		# Unit: s

U       = NcFile.variables['uwnd']	# U[times,lat,lon]
V       = NcFile.variables['vwnd']

# define the initial location of particles
nn = 9
n_box_max	= nn*nn
box_lon	= np.zeros(n_box_max)
box_lat	= np.zeros(n_box_max)

for i in range(nn):
    for j in range(nn):
        box_lon[i*nn+j]	= 120.0 + 0.1*(i+1)
        box_lat[i*nn+j]	= 15.0 + 0.1*(j+1)


res_lat = np.zeros((len(time)+1,len(box_lat)))
res_lon = np.zeros((len(time)+1,len(box_lat)))

res_lat[0,:] = np.array(box_lat)
res_lon[0,:] = np.array(box_lon)


#=====================================================
# functions
def find_i(a,b):
    idx = np.argmin( abs(np.array(b)-a) )
    if b[idx]>a: idx = idx-1
    return idx


def Calc_dist(lon1,lat1,lon2,lat2):
    lon1 = math.radians(lon1)
    lat1 = math.radians(lat1)

    lon2 = math.radians(lon2)
    lat2 = math.radians(lat2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (math.sin(dlat/2))**2 + math.cos(lat1) * math.cos(lat2) * (math.sin(dlon/2))**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = Re * c
    return d

# functions
def find_wind(wind,curr_lon,curr_lat,lon,lat):
    i_lon = find_i(curr_lon,lon)        # find the nearest index of lon that less than curr_lon
    i_lat = find_i(curr_lat,lat)        # find the nearest index of lat that less than curr_lat
    uv	= []
    distance	= []
    w	= []
    for i in range(2):
        ii = i + int(i_lon)
        if ii==IIPAR+1: ii = 0	# once the lon is bigger than 360 degree, change lon to 0 degree
        for j in range(2):
            jj   = j + i_lat
            num  = ii*2 + jj
            uv.append(wind[jj,ii])
            distance.append( Calc_dist(curr_lon,curr_lat,lon[ii],lat[jj]) )	# calculate the distance between two points

    for i in range(len(distance)):
        w.append( (1/distance[i]) / ( sum(1/np.array(distance[:])) ) )       

    Interplt_wind = np.mean( np.array(w[:]) * np.array(uv[:]) )

    return Interplt_wind


def m2lon(lat, lon, dx):
    lat = math.radians(lat)
    lon = math.radians(lon)

    D_lon = dx/(Re*math.cos(lat))
    D_lon = math.degrees(D_lon)

    return D_lon


def m2lat(lat, lon, dy):
    lat = math.radians(lat)
    lon = math.radians(lon)

    D_lat = math.degrees(dy/Re)

    return D_lat


#======================================================
# integrate ODE
for it in range(len(time)):
    print(it)
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
        
        res_lon[it+1,i] = box_lon[i]
        res_lat[it+1,i] = box_lat[i]


# save data box_lon, box_lat into a file, so that we can use this data to plot figures.

lons	= np.reshape( res_lon , len(res_lon[:,0])*len(res_lon[0,:]) )
lats	= np.reshape( res_lat , len(res_lat[:,0])*len(res_lat[0,:]) )

res	= np.vstack((lons, lats)).T
np.savetxt("lat_lon.txt", res, fmt="%f")




