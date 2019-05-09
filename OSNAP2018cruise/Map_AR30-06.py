from AR30_funcs import *

datadir='/home/isabela/Documents/cruises/OSNAP2018_AR30-06/Chipods/AR30/Data/'

from map_funcs import *

ctd=pickle.load(open('OSNAP2018cruise/pickles/CTD_1mbin.pickle','rb'))

def mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed):

    if lon_start< 0 and lon_end < 0:
        lon_0= - (abs(lon_end)+abs(lon_start))/2.0
        lat_0= - (abs(lat_end)+abs(lat_start))/2.0
    else:
        lon_0=(abs(lon_end)+abs(lon_start))/2.0
        lat_0= - (abs(lat_end)+abs(lat_start))/2.0

    map = Basemap(lat_0=lat_0,lon_0=lon_0,llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='h',projection='stere')

    x, y = map(lon,lat)

    # CS2 = map.contourf(x,y,bathySmoothed,[landlev,5e3],colors='lightgrey')
    map.fillcontinents(color='lightgrey')
    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000,-500,-200],colors='lightgrey')
    # NA = map.contour(x,y,bathySmoothed,[0],colors='black')



    CS0.axis='tight'

    return map,x,y,CS0

from matplotlib import gridspec

landlev=0
def makeMap():

    lat_start=57
    lat_end  =68
    lon_start=-52
    lon_end  =-20
    # lon_start=-62
    # lon_end  =-23

    figure(figsize=(10,5))


    """Get the etopo1 data"""
    etopo1name='/home/isabela/Documents/projects/bathymetry/etopo1/ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]

    res = findSubsetIndices(lat_start-5,lat_end+5,lon_start-40,lon_end+10,lats,lons)

    lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)

    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    plt.subplot(gs[0])

    map,x,y,CS0=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)

    map.drawmeridians(range(-60,lon_end+10,10),linewidth=0.001,labels=[0,0,0,1])
    map.drawparallels(arange(56,70,2),linewidth=0.001,labels=[1,0,0,0])

    cmapo=cm.get_cmap('rainbow', len(seclab.keys()))

    for ii,ss in enumerate(seclab):
        if (('sec' in ss) | ('cal' in ss)) & ('b' not in ss) & ('moor' not in ss) & ('cal 2' not in ss):
            map.plot(ctd.lon[seclab[ss]].values,ctd.lat[seclab[ss]].values,'.',label=ss,latlon=True,color=cmapo(ii))



    ############# Zoom-in location

    lat_start=58
    lat_end  =61
    lon_start=-49.5
    lon_end  =-40


    map.plot([lon_start,lon_start,lon_end,lon_end,lon_start],[lat_start,lat_end,lat_end,lat_start,lat_start],latlon=True,linewidth=2,color='k')

    plt.subplot(gs[1])

    map2,x,y,CS0=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)

    for ii,ss in enumerate(seclab):
        if 'moorings' not in ss:
            map2.plot(ctd.lon[seclab[ss]].values,ctd.lat[seclab[ss]].values,'o',label=ss,latlon=True,color=cmapo(ii))

    clabels=clabel(CS0,fmt='%1.0f',colors='k')

    map2.drawmeridians(range(-49,-39,4),linewidth=0.001,labels=[0,0,0,1])
    map2.drawparallels(arange(lat_start,lat_end,1),linewidth=0.001,labels=[1,0,0,0])

    legend(loc=(1.05,-0.1))
    savefig(figdir+'AR30/Map_AR30-06_sections.pdf',bbox_inches='tight')
    savefig(figdir+'AR30/Map_AR30-06_sections.png',bbox_inches='tight')

makeMap()
