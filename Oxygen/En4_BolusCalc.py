from firstfuncs_1618 import *
from map_funcs import *



def mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed):

    lon_0= (lon_end + lon_start)/2.0
    lat_0= - (abs(lat_end)+abs(lat_start))/2.0


    map = Basemap(lat_0=lat_0,lon_0=lon_0,llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='h',projection='stere')

    x, y = map(lon,lat)
    # lon4,lat4=np.meshgrid(en4.lon.values,en4.lat.values)
    # x4,y4=map(lon4,lat4)
    # CS1 = map.contourf(x4,y4,en1416_janmar.values,51,cmap=cm.YlGn,vmin=27.3,vmax=27.8,extend='both')

    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000],colors='grey')

    # map.drawcoastlines()
    map.fillcontinents('darkgrey')

    # eastind=osnap.LONGITUDE>-43
    # map.plot(osnap.LONGITUDE.values[eastind],osnap.LATITUDE.values[eastind],color='k',latlon=True,linewidth=4)

    return map,x,y,CS0

def SimpleMap():

    lat_start=55
    lat_end  =67
    lon_start=-53
    lon_end  =-18

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

    return map


##################################################################################
################## MAP OUT WHERE THE PROFILES ARE... #############################
##################################################################################
dat={}
for dd in range(1,10):
    alldat=xr.open_dataset('/home/isabela/DATA/EN4/profiles/EN.4.2.1.f.profiles.g10.20150'+str(dd)+'.nc')
    dat[dd]=alldat.where(alldat.LATITUDE>55).where(alldat.LONGITUDE>-50).where(alldat.LONGITUDE<-30).where(alldat.LATITUDE<68)

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/'


for dd in range(1,10):
    figure()
    map=SimpleMap()
    map.plot(dat[dd].LONGITUDE.values,dat[dd].LATITUDE.values,'o',latlon=True,color='C'+str(dd))
    title('2015-'+str(dd))
    savefig(figdir+'Oxygen/EN4/EN4_profile_locations_2015-'+str(dd)+'.png')


##################################################################################
############### MAP OUT HOW THICK LOW PPV DENSITY LAYER IS... ####################
##################################################################################


grd=xr.open_dataset('/home/isabela/DATA/EN4/analyses/EN.4.2.1.f.analysis.g10.201501.nc').sel(lon=slice(360-50,360-30)).sel(lat=slice(55,68))

for dd in range(2,10):
    alldat=xr.open_dataset('/home/isabela/DATA/EN4/analyses/EN.4.2.1.f.analysis.g10.20150'+str(dd)+'.nc')
    grd=xr.concat([grd,alldat.sel(lon=slice(360-50,360-30)).sel(lat=slice(55,68))],dim='time')
grd['temperature']=grd['temperature']-273

def add_PDEN(xray):
        pres=gsw.p_from_z(-xray.depth,60)
        presmat,latmat,lonmat=meshgrid(pres,xray.lat,xray.lon,indexing='ij')
        PD_out=NaN*xray['salinity']
        for ii,tt in enumerate(grd.time):
            SA_out=gsw.SA_from_SP(xray['salinity'].sel(time=tt),presmat,lonmat,latmat)
            CT_out=gsw.CT_from_pt(SA_out,xray['temperature'].sel(time=tt))
            PD_out[ii,:,:,:]=gsw.sigma0(SA_out,CT_out)
        xray['pden']=(('time', 'depth', 'lat', 'lon'),PD_out)
        return xray


grd=add_PDEN(grd)

#Interpolate onto fine depth grid to up resolution of thickness
fgrd=grd.interp(depth=range(0,3000,2))

# Find depth of density boundaries
dup=27.73
dlow=27.77

lbnd=fgrd.depth.where(fgrd.pden>dlow).min(dim='depth')
ubnd=fgrd.depth.where(fgrd.pden>dup).min(dim='depth')
ubnd.sel(lon=slice(320,325)).sel(lat=slice(60,61)).mean(dim='lon').mean(dim='lat').plot()

thick=lbnd-ubnd

ubnd.isel(time=2).plot()
lbnd.isel(time=5).plot()

thick.isel(time=5).plot(vmin=0,vmax=500)
