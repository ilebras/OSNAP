from aux_funcs import *

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

#make a map that depicts the relevant measurement sites
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.basemap import cm as bcm
import os, sys, datetime, string

# Based off of Trond Kristiansen's createMapsEtopo1.py
# -*- coding: utf-8 -*-

# Laplace filter a 2D field with mask

# ROMS type maske,
#        M = 1 pÃ¥ grid-celler som brukes
#        M = 0 pÃ¥ grid-celler som maskeres vekk
# Default er uten maske
import numpy as np

osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

def laplace_X(F,M):
    """1D Laplace Filter in X-direction (axis=1)"""

    jmax, imax = F.shape

    # Add strips of land
    F2 = np.zeros((jmax, imax+2), dtype=F.dtype)
    F2[:, 1:-1] = F
    M2 = np.zeros((jmax, imax+2), dtype=M.dtype)
    M2[:, 1:-1] = M

    MS = M2[:, 2:] + M2[:, :-2]
    FS = F2[:, 2:]*M2[:, 2:] + F2[:, :-2]*M2[:, :-2]

    return np.where(M > 0.5, (1-0.25*MS)*F + 0.25*FS, F)

def laplace_Y(F,M):
    """1D Laplace Filter in Y-direction (axis=1)"""

    jmax, imax = F.shape

    # Add strips of land
    F2 = np.zeros((jmax+2, imax), dtype=F.dtype)
    F2[1:-1, :] = F
    M2 = np.zeros((jmax+2, imax), dtype=M.dtype)
    M2[1:-1, :] = M

    MS = M2[2:, :] + M2[:-2, :]
    FS = F2[2:, :]*M2[2:, :] + F2[:-2, :]*M2[:-2, :]

    return np.where(M > 0.5, (1-0.25*MS)*F + 0.25*FS, F)


# The mask may cause laplace_X and laplace_Y to not commute
# Take average of both directions

def laplace_filter(F, M=None):
    if M == None:
        M = np.ones_like(F)
    return 0.5*(laplace_X(laplace_Y(F, M), M) +
                laplace_Y(laplace_X(F, M), M))


import matplotlib
from numpy import ma
import pylab as pl



"""A set of utility functions for matplotlib"""

# ------------------
# Plot land mask
# ------------------

def landmask(M, color='0.8'):

   # Make a constant colormap, default = grey
   constmap = pl.matplotlib.colors.ListedColormap([color])

   jmax, imax = M.shape
   # X and Y give the grid cell boundaries,
   # one more than number of grid cells + 1
   # half integers (grid cell centers are integers)
   X = -0.5 + pl.arange(imax+1)
   Y = -0.5 + pl.arange(jmax+1)

   # Draw the mask by pcolor
   M = ma.masked_where(M > 0, M)
   pl.pcolor(X, Y, M, shading='flat', cmap=constmap)

# -------------
# Colormap
# -------------

# Colormap, smlgn. med Rob Hetland

def LevelColormap(levels, cmap=None):
    """Make a colormap based on an increasing sequence of levels"""

    # Start with an existing colormap
    if cmap == None:
        cmap = pl.get_cmap()

    # Spread the colours maximally
    nlev = len(levels)
    S = pl.arange(nlev, dtype='float')/(nlev-1)
    A = cmap(S)

    # Normalize the levels to interval [0,1]
    levels = pl.array(levels, dtype='float')
    L = (levels-levels[0])/(levels[-1]-levels[0])

    # Make the colour dictionary
    R = [(L[i], A[i,0], A[i,0]) for i in xrange(nlev)]
    G = [(L[i], A[i,1], A[i,1]) for i in xrange(nlev)]
    B = [(L[i], A[i,2], A[i,2]) for i in xrange(nlev)]
    cdict = dict(red=tuple(R),green=tuple(G),blue=tuple(B))

    # Use
    return matplotlib.colors.LinearSegmentedColormap(
        '%s_levels' % cmap.name, cdict, 256)

predir='/home/isabela/DATA/bathymetry/etopo1/'

def findSubsetIndices(min_lat,max_lat,min_lon,max_lon,lats,lons):

    """Array to store the results returned from the function"""
    res=np.zeros((4),dtype=np.float64)
    minLon=min_lon; maxLon=max_lon

    distances1 = []; distances2 = []
    indices=[]; index=1

    for point in lats:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    distances1 = []; distances2 = []; index=1

    for point in lons:
        s1 = maxLon-point # (vector subtract)
        s2 = minLon-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
    minJ=indices[1][2]
    maxJ=indices[0][2]
    minI=indices[3][2]
    maxI=indices[2][2]

    res[0]=minI; res[1]=maxI; res[2]=minJ; res[3]=maxJ;
    return res







# def makeMap(zoomlev):
#     if zoomlev=='CF':
#         lat_start=58.5
#         lat_end  =61.5
#         lon_start=-45
#         lon_end  =-39
#     elif zoomlev=='fullsouth':
#         lat_start=58
#         lat_end  =62
#         lon_start=-50
#         lon_end  =-40
#     elif zoomlev=='GR':
#         lat_start=58
#         lat_end  =70
#         lon_start=-55
#         lon_end  =-30
#     else:
#         lat_start=53
#         lat_end  =66
#         lon_start=-60
#         lon_end  =-20
#
#     figure(figsize=(8,8))
#
#     """Get the etopo1 data"""
#     etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
#     etopo1 = Dataset(etopo1name,'r')
#
#     lons = etopo1.variables["x"][:]
#     lats = etopo1.variables["y"][:]
#
#     res = findSubsetIndices(lat_start-5,lat_end+5,lon_start-40,lon_end+10,lats,lons)
#
#     lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
#     bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
#     bathySmoothed = laplace_filter(bathy,M=None)
#
#     levels=array([-7000,-6000,-5000,-3000, -2000, -1500, -1000,-500, -250])
#
#     if lon_start< 0 and lon_end < 0:
#         lon_0= - (abs(lon_end)+abs(lon_start))/2.0
#     else:
#         lon_0=(abs(lon_end)+abs(lon_start))/2.0
#
#     # print 'Center longitude ',lon_0
#
#     map = Basemap(llcrnrlat=lat_start,urcrnrlat=lat_end,
#                 llcrnrlon=lon_start,urcrnrlon=lon_end,
#                 resolution=None,projection='stere',
#                 lat_0=(lat_start+lat_end)/2,lon_0=lon_0)
#
#     x, y = map(lon,lat)
#
#     CS2 = map.contourf(x,y,bathySmoothed,arange(0,2.5e3,100),
#                       cmap=cm.Oranges_r)
#
#     CS1 = map.contourf(x,y,bathySmoothed,arange(-4e3,0,100),
#                       cmap=cm.Blues_r)
#
#
#
#     if (zoomlev=='CF') | (zoomlev=='fullsouth'):
#         CS0 = map.contour(x,y,bathySmoothed,levels,
#                            colors='grey')
#         clabel(CS0,fmt='%1.0f')
#
#
#         CS0.axis='tight'
#
#         map.drawmeridians(range(lon_start,lon_end,2),linewidth=0.001,labels=[0,0,0,1])
#         map.drawparallels(arange(59,lat_end,1),linewidth=0.001,labels=[1,0,0,0])
#         continent = map.contour(x,y,bathySmoothed,[0],colors='grey')
#     else:
#         map.drawmeridians(range(lon_start,lon_end,5),linewidth=0.001,labels=[0,0,0,1])
#         map.drawparallels(arange(54,lat_end+1,2),linewidth=0.001,labels=[1,0,0,0])
#
#
#     if zoomlev=='GR':
#         CS0 = map.contour(x,y,bathySmoothed,[-1e3,-500,-200,-100],
#                            colors='black')
#         clabel(CS0,fmt='%1.0f')
#     return map
#
# map=makeMap('GR')
# savefig('../figures/Map/Greenlandbathy.pdf')

#
# map=makeMap('zoomout')
# savefig('../figures/Map/Map_zoomout.pdf',bbox_inches='tight')
# lat_start=58.5
# lat_end  =61.5
# lon_start=-45
# lon_end  =-39
# map.plot([lon_start,lon_start,lon_end,lon_end,lon_start],[lat_start,lat_end,lat_end,lat_start,lat_start],latlon=True,linewidth=2,color='k')
# map.plot(CFlon,CFlat,latlon=True,linewidth=2,color='k')
# savefig('../figures/Map/Map_zoomout_wbox.pdf',bbox_inches='tight')
