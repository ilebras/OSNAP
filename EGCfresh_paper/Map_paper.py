from map_funcs import *

dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_'+versname+'.pickle','rb'))

u=dat['across track velocity']*cos(theta)-dat['along track velocity']*sin(theta)

v=dat['across track velocity']*sin(theta)+dat['along track velocity']*cos(theta)

#plot mean over top 150m
umean=u[:,:76,:].mean(dim='date').mean(dim='depth')
vmean=v[:,:76,:].mean(dim='date').mean(dim='depth')

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
    map.contourf(x,y,bathySmoothed,linspace(-4500,0,50),cmap=cm.Blues_r)
    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000,-200],colors='grey')
    # CS2 = map.contour(x,y,bathySmoothed,[landlev],colors='k')
    # CS2 = map.contourf(x,y,bathySmoothed,[landlev,5e3],colors='darkgrey')
    CS2=0
    map.drawcoastlines()
    # map.drawcoastlines()

    CS0.axis='tight'

    return map,x,y,CS0,CS2

landlev=0
def SimplePaperMap():

    lat_start=55
    lat_end  =67
    lon_start=-52
    lon_end  =-18
    # lon_start=-62
    # lon_end  =-23

    figure(figsize=(10,5))


    """Get the etopo1 data"""
    etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]

    res = findSubsetIndices(lat_start-5,lat_end+5,lon_start-40,lon_end+10,lats,lons)

    lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)

    # print 'Center longitude ',lon_0
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.75, 1])
    plt.subplot(gs[0])

    map,x,y,CS0,CS2=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)

    # subplot(121)

    map.drawmeridians(range(-60,lon_end+10,10),linewidth=0.001,labels=[0,0,0,1])
    map.drawparallels(arange(56,67,2),linewidth=0.001,labels=[1,0,0,0])

    ############# Zoom-in location

    lat_start=59
    lat_end  =61
    lon_start=-45
    lon_end  =-40.5


    map.plot([lon_start,lon_start,lon_end,lon_end,lon_start],[lat_start,lat_end,lat_end,lat_start,lat_start],latlon=True,linewidth=2,color='k')
    # map.plot(CFlon,CFlat,latlon=True,linewidth=2,color='k')
    plt.subplot(gs[1])
    # subplot(122)

    map2,x,y,CS0,CS2=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)

    map2.plot(CFlon,CFlat,'ko', markersize=7,zorder=98,latlon=True)
    qq=map2.quiver(CFlon,CFlat,umean,vmean,latlon=True,linewidth=0.5,scale=2,zorder=100)
    # colorbar(CS0)
    clabels=clabel(CS0,fmt='%1.0f',colors='k')
    [txt.set_backgroundcolor('white') for txt in clabels]
    [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.1,alpha=1)) for txt in clabels]

    map2.drawmeridians(range(lon_start,-39,2),linewidth=0.001,labels=[0,0,0,1])
    map2.drawparallels(arange(59.5,lat_end,0.5),linewidth=0.001,labels=[1,0,0,0])

    savefig('../figures/paperfigs/Map_rev.pdf',bbox_inches='tight')


SimplePaperMap()
