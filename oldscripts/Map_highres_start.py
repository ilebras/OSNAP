from aux_funcs import *

#make a map that depicts the relevant measurement sites
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.basemap import cm as bcm
import os, sys, datetime, string
import laplaceFilter
import mpl_util


# Based off of Trond Kristiansen's createMapsEtopo1.py

# OOI mooring location: 59.93°N, 39.47°W

predir='/home/isabela/Documents/bathymetry/etopo1/'

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


lat_start=58
lat_end  =64

lon_start=-46
lon_end  =-40

name='Zoom-out'

def makeMap():
    figure(figsize=(8,8))

    """Get the etopo1 data"""
    etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]

    res = findSubsetIndices(lat_start-5,lat_end+5,lon_start-40,lon_end+10,lats,lons)

    lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
    # print "Extracted data for area %s : (%s,%s) to (%s,%s)"%(name,lon.min(),lat.min(),lon.max(),lat.max())
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplaceFilter.laplace_filter(bathy,M=None)

    levels=[-6000,-5000,-3000, -2000, -1500, -1000,-500, -400, -300,
             -200,  -100,  -50, 0, 100,200,300,400,500, 1e3, 2e3]

    if lon_start< 0 and lon_end < 0:
        lon_0= - (abs(lon_end)+abs(lon_start))/2.0
    else:
        lon_0=(abs(lon_end)+abs(lon_start))/2.0

    # print 'Center longitude ',lon_0

    map = Basemap(llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution=None,projection='lcc',
                lat_1=lat_start,lon_0=lon_0)

    x, y = map(lon,lat)
    # map.drawcoastlines()
    # map.drawcountries()
    # #     map.fillcontinents(color='grey')
    # map.drawmeridians(np.arange(lons.min(),lons.max(),2),labels=[0,0,0,1])
    # map.drawparallels(np.arange(lats.min(),lats.max(),1),labels=[1,0,0,0])

    # map.bluemarble()
    # CS1 = map.contourf(x,y,bathySmoothed,levels,
    #                   cmap=cm.Blues_r,
    #                   extend='upper',
    #                   alpha=1.0,
    #                   origin='lower')
    CS0 = map.contour(x,y,bathySmoothed,levels,
                       colors='k')



    map.plot(CFlon[:-1],CFlat[:-1],'ro',markersize=8,latlon=True,label='Cape Farewell Array')
    map.plot(CFlon[-1],CFlat[-1],'mo',markersize=8,latlon=True,label='NOC M1')
    # map.plot(-39.47,59.93,'o',color='orange',markersize=8,latlon=True,label='OOI mooring')#ooi mooring


    CS0.axis='tight'

    # plt.legend(loc=(1.05,0.5),numpoints=1)
    # plt.title(name)
    # plotfile='../figures/overview/map_'+str(name)+'.pdf'
    # plt.savefig(plotfile,dpi=300,orientation='portrait',bbox_inches='tight')




makeMap()
