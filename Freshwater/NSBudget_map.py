from firstfuncs_1618 import *
from map_funcs import *
import cartopy
import cartopy.crs as ccrs

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Schematic/'

osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')

# # #### load osnap data and cut out the eastern portion
lonbnd=-44
dat=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.nc')
osnap_obs=dat.sel(LONGITUDE=slice(lonbnd,0))
osnap_obs['LATITUDE']=dat['LATITUDE'].values[dat.LONGITUDE>lonbnd]
fs_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_fs_xray_2001.nc')
bso_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_bso_xray_2001.nc')


central_lon, central_lat = -20, 70
extent = [-45, 5, 50, 86]

"""Get the etopo1 data"""
etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
etopo1 = Dataset(etopo1name,'r')

lons = etopo1.variables["x"][:]
lats = etopo1.variables["y"][:]
res = findSubsetIndices(extent[2]-20,extent[3]+20,extent[0]-80,extent[1]+80,lats,lons)
lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
bathySmoothed = laplace_filter(bathy,M=None)

import matplotlib.ticker as mticker

def map_setup():
    figure(figsize=(14,10))
    ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))
    ax.set_extent(extent)
    # g1=ax.gridlines()
    # g1.xlocator=mticker.FixedLocator([-120,-80,-50,-20,10,40,80])
    # g1.ylocator=mticker.FixedLocator([45,65,75,90])
    data_crs = ccrs.PlateCarree()
    ax.contourf(lon,lat,bathySmoothed,levels=arange(0,5000,1000),colors='lightgrey',transform=data_crs)
    bb=ax.contourf(lon,lat,bathySmoothed,levels=arange(-5000,1,250),cmap=cm.Greys_r,transform=data_crs)
    colorbar(bb,ticks=arange(-4000,1,1000))
    ax.contour(lon,lat,bathySmoothed,levels=0,colors='k',transform=data_crs)
    lwi=5
    c1='limegreen'
    ax.plot(fs_obs.LONGITUDE,fs_obs.LATITUDE,linewidth=lwi,color=c1,transform=data_crs)
    ax.plot(bso_obs.LONGITUDE,bso_obs.LATITUDE,linewidth=lwi,color=c1,transform=data_crs)
    ax.plot(osnap_obs.LONGITUDE,osnap_obs.LATITUDE,linewidth=lwi,color=c1,transform=data_crs)
    ax.plot(fs.LONGITUDE,fs.LATITUDE,linewidth=lwi,color='purple',transform=data_crs)
    ax.plot(ns.LONGITUDE,ns.LATITUDE,linewidth=lwi,color='purple',transform=data_crs)
    ax.plot(bso.LONGITUDE,bso.LATITUDE,linewidth=lwi,color='purple',transform=data_crs)
    ax.plot(osnap.LONGITUDE,osnap.LATITUDE,linewidth=lwi,color='purple',transform=data_crs)
    # ax.text(-25,85,'20$^{\circ}$W',transform=data_crs,fontsize=14,backgroundcolor='w')
    # ax.text(-51,60,'50$^{\circ}$W',transform=data_crs,fontsize=14,backgroundcolor='w')
    # ax.text(7,61,'10$^{\circ}$E',transform=data_crs,fontsize=14,backgroundcolor='w')
    # ax.text(-81,78,'80$^{\circ}$W',transform=data_crs,fontsize=14,backgroundcolor='w')
    # ax.text(35,78,'40$^{\circ}$E',transform=data_crs,fontsize=14,backgroundcolor='w')
    # ax.text(-57,64,'65$^{\circ}$N',transform=data_crs,fontsize=14,backgroundcolor='w')
    # ax.text(-73,74,'75$^{\circ}$N',transform=data_crs,fontsize=14,backgroundcolor='w')
    savefig(figdir+'NordicSea_Map_newcol.png',bbox_inches='tight')
    savefig(figdir+'NordicSea_Map_newcol.pdf',bbox_inches='tight')

map_setup()

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
