from map_funcs import *


[FLA,FLB]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/FLMAB_xrays.pickle','rb'))

modlon=xr.open_dataset(datadir+'OSNAP2016recovery/ItalianModel/lon_ItalianModel.nc')
modlat=xr.open_dataset(datadir+'OSNAP2016recovery/ItalianModel/lat_ItalianModel.nc')
LOCOlon=-39.5
LOCOlat=59

ssh=xr.open_dataset(datadir+'OSNAP2016recovery/ItalianModel/SSH_2014_2016_ItalianModel.nc')

mlon=modlon.nav_lon.values
mlat=modlat.nav_lat.values
mssh=ssh.mean(dim='time_counter').zos.values
mssh[mssh>=0]=NaN

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
    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000],colors='grey')

    CS2=map.contourf(mlon,mlat,mssh,20,cmap=cm.PuOr,latlon='True')

    map.drawcoastlines()


    return map,x,y,CS0,CS2

def SimplePaperMap():

    lat_start=57.5
    lat_end = 66
    lon_start=-44.5
    lon_end  =-20

    fig,gs=subplots(1,2,figsize=(10,5),gridspec_kw = {'width_ratios':[1.25,1]})


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
    # gs = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1])
    plt.subplot(gs[0])

    map,x,y,CS0,CS2=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)

    map.drawmeridians(range(-60,lon_end+10,10),linewidth=0.001,labels=[0,0,0,1])
    map.drawparallels(arange(56,67,2),linewidth=0.001,labels=[1,0,0,0])

    ############# Zoom-in location

    lat_start=58.5
    lat_end  =61
    lon_start=-44
    lon_end  =-38


    map.plot([lon_start,lon_start,lon_end,lon_end,lon_start],[lat_start,lat_end,lat_end,lat_start,lat_start],latlon=True,linewidth=2,color='k')

    plt.subplot(gs[1])


    map2,x,y,CS0,CS2=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)


    map2.plot(CFlon[4:],CFlat[4:],'ko', markersize=7,zorder=98,latlon=True)

    map2.plot(FLA.lon.mean(),FLA.lat.mean(),'ko', markersize=7,zorder=98,latlon=True)
    map2.plot(FLB.lon.mean(),FLB.lat.mean(),'ko', markersize=7,zorder=98,latlon=True)

    map2.plot(LOCOlon,LOCOlat,'ko', markersize=7,zorder=98,latlon=True)
    # map2.plot(-39.5,59,'ko', markersize=7,zorder=98,latlon=True)
    x1,x2=map2(-42.2,60.25)
    gs[1].text(x1,x2,'OSNAP',fontsize=12,bbox=dict(facecolor='white'))
    x3,x4=map2(-39.8,60)
    gs[1].text(x3,x4,'OOI',fontsize=12,bbox=dict(facecolor='white'))
    x5,x6=map2(-39.8,59.2)
    gs[1].text(x5,x6,'LOCO',fontsize=12,bbox=dict(facecolor='white'))
    # colorbar(CS0)
    # clabels=clabel(CS0,fmt='%1.0f',colors='k')
    # [txt.set_backgroundcolor('white') for txt in clabels]
    # [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.1,alpha=1)) for txt in clabels]

    map2.drawmeridians(range(-44,-39,2),linewidth=0.001,labels=[0,0,0,1])
    map2.drawparallels(arange(59,lat_end,1),linewidth=0.001,labels=[1,0,0,0])

    cbaxes = fig.add_axes([0.95, 0.3, 0.02, 0.4])
    cb=colorbar(CS2,ticks=arange(-1.2,0,0.2),label='SSH [m]',cax=cbaxes)

    savefig(figdir+'MixedLayer/pres/Map_ssh.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/pres/Map_ssh.png',bbox_inches='tight',dpi=300)

SimplePaperMap()
