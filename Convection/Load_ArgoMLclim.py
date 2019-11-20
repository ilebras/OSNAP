from aux_funcs import *
from map_funcs import *

datadir

alldat=xr.open_dataset(glob.glob(datadir+'aux_data/ArgoMLclim/*.nc')[0])

alldat['date']=date_from_matlab(alldat.profiledate.values)

alldat=alldat.rename({'iNPROF':'date'})

dat=alldat.where(alldat.profilelat>50).where(alldat.profilelat<70).where(alldat.profilelon>-60).where(alldat.profilelon<-10)

dat=dat.sortby('date')

def mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed):

    lon_0= (lon_end + lon_start)/2.0
    lat_0= - (abs(lat_end)+abs(lat_start))/2.0


    map = Basemap(lat_0=lat_0,lon_0=lon_0,llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='l',projection='stere')

    x, y = map(lon,lat)
    # lon4,lat4=np.meshgrid(en4.lon.values,en4.lat.values)
    # x4,y4=map(lon4,lat4)
    # CS1 = map.contourf(x4,y4,en1416_janmar.values,51,cmap=cm.YlGn,vmin=27.3,vmax=27.8,extend='both')

    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000],colors='grey')

    # map.drawcoastlines()
    map.fillcontinents('darkgrey')

    # eastind=osnap.LONGITUDE>-43
    # map.plot(osnap.LONGITUDE.values[eastind],osnap.LATITUDE.values[eastind],color='k',latlon=True,linewidth=4)

    return map,x,y,CS0#,CS1


dmin=27
dmax=28
dencmap=cm.Dark2_r
bounds = array([dmin,d1,d2,d3,dmax])
norm = mpl.colors.BoundaryNorm(bounds, dencmap.N)

MLcmap=cm.rainbow
MLbounds = array([0,250,500,1000,1500,2500])
MLnorm = mpl.colors.BoundaryNorm(MLbounds, MLcmap.N)

salcmap=cm.Dark2
salbounds = array([34,34.8,34.93,34.98,35.1])
salnorm = mpl.colors.BoundaryNorm(salbounds, salcmap.N)

tmpcmap=cm.Set1
tmpbounds = array([-1,3,4.2,4.7,6,10])
tmpnorm = mpl.colors.BoundaryNorm(tmpbounds, tmpcmap.N)


month1='-2-'
month2='-4-'

def SimpleMap(which='None'):

    lat_start=50
    lat_end  =66
    lon_start=-63
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
    # [map.plot(ooi_lon[ll],ooi_lat[ll],'o',latlon=True) for ll in ooi_lat]
    map.plot(CFlon,CFlat,'k*',latlon=True,zorder=101,markersize=20)
    if which=='depth':
            map.scatter(dat.profilelon.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,dat.profilelat.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,
            c=dat.da_mld.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,latlon=True,norm=MLnorm,cmap=MLcmap,zorder=100)
            colorbar(label='ML depth [m]')
    elif which=='dens':
            map.scatter(dat.profilelon.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,dat.profilelat.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,
            c=dat.da_mlpd.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,latlon=True,norm=norm,cmap=dencmap,zorder=100)
            colorbar(label='ML pot. density')
    elif which=='sal':
            map.scatter(dat.profilelon.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,dat.profilelat.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,
            c=dat.da_mls.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,latlon=True,norm=salnorm,cmap=salcmap,zorder=100)
            colorbar(label='ML salinity')
    elif which=='tmp':
            map.scatter(dat.profilelon.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,dat.profilelat.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,
            c=dat.da_mlt.sel(date=slice(str(yy)+month1+'1',str(yy)+month2+'1')).values,latlon=True,norm=tmpnorm,cmap=tmpcmap,zorder=100)
            colorbar(label='ML pot. temperature')
    return map


for yy in range(2013,2019):
    figure(figsize=(13,5))
    subplot(121)
    map=SimpleMap(which='sal')
    title(yy)
    subplot(122)
    map=SimpleMap(which='tmp')
    title(yy)
    savefig(figdir+'MixedLayer/ArgoMLclim/ArgoML_saltmp_FebMar_'+str(yy)+'.png',bbox_inches='tight')



for yy in range(2013,2019):
    figure(figsize=(13,5))
    subplot(121)
    map=SimpleMap(which='depth')
    title(yy)
    subplot(122)
    map=SimpleMap(which='dens')
    title(yy)
    savefig(figdir+'MixedLayer/ArgoMLclim/ArgoML_dpthden_FebMar_'+str(yy)+'.png',bbox_inches='tight')
