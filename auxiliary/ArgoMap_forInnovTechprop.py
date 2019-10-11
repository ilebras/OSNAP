from map_funcs import *

def mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed):

    lon_0= (lon_end + lon_start)/2.0
    lat_0= - (abs(lat_end)+abs(lat_start))/2.0


    map = Basemap(lat_0=lat_0,lon_0=lon_0,llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='h',projection='stere')

    x, y = map(lon,lat)
    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000],colors='grey')
    # map.drawcoastlines()
    map.fillcontinents('darkgrey')

    return map,x,y,CS0

import glob
datadir='/home/isabela/Documents/proposals/2019Oct_InternalInnovation_floats/Argo_data/DataSelection_20191002_184329_8701800/'
import xarray as xr

datlist=glob.glob(datadir+'*trajectory*')

def SimpleMap():

    lat_start=60
    lat_end  =75
    lon_start=-55
    lon_end  =-10

    """Get the etopo1 data"""
    etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]

    res = findSubsetIndices(lat_start-5,lat_end+5,lon_start-40,lon_end+10,lats,lons)

    lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)

    map,x,y,CS0=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)
    map.drawmeridians(range(lon_start+3,lon_end+10,10),labels=[0,0,0,1],linewidth=0.0001,fontsize=12)
    map.drawparallels(arange(lat_start,lat_end+2,5),labels=[1,0,0,0],linewidth=0.0001,fontsize=12)


    for dd in datlist:
        dat=xr.open_dataset(dd)
        map.plot(dat.LONGITUDE.values,dat.LATITUDE.values,'.-',latlon=True,alpha=0.4)

    savefig('/home/isabela/Documents/proposals/2019Oct_InternalInnovation_floats/Argomap.pdf')
    savefig('/home/isabela/Documents/proposals/2019Oct_InternalInnovation_floats/Argomap.png',dpi=300)
    return map


SimpleMap()

# Interesting floats
i1='6901911'
i2='6902728'

glob.glob(datadir+'*'+i1+'*')

proflist=glob.glob(datadir+'*profiles*')


dtest=xr.open_dataset(proflist[0])


for dd in datlist:
    dat=xr.open_dataset(dd)
    figure()
    plot(dat.LONGITUDE.values,dat.LATITUDE.values,'.-')
    title(dd[-30:])
    xlim(-55,-10)
    ylim(60,75)
