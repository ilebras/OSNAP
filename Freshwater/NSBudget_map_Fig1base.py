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
predir='/home/isabela/DATA/bathymetry/etopo1/'
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
    data_crs = ccrs.PlateCarree()
    ax.contourf(lon,lat,bathySmoothed,levels=arange(0,10000,5000),colors='k',transform=data_crs,zorder=2)
    bb=ax.contour(lon,lat,bathySmoothed,levels=[-400],colors='grey',transform=data_crs,zorder=2)

    lwi=5
    c1='limegreen'
    ax.plot(fs_obs.LONGITUDE,fs_obs.LATITUDE,linewidth=lwi,color=c1,transform=data_crs,zorder=1)
    ax.plot(bso_obs.LONGITUDE,bso_obs.LATITUDE,linewidth=lwi,color=c1,transform=data_crs,zorder=1)
    ax.plot(osnap_obs.LONGITUDE,osnap_obs.LATITUDE,linewidth=lwi,color=c1,transform=data_crs,zorder=1)
    nscol='darkgreen'
    ax.plot(fs.LONGITUDE,fs.LATITUDE,linewidth=lwi,color=nscol,transform=data_crs,zorder=1)
    ax.plot(ns.LONGITUDE,ns.LATITUDE,linewidth=lwi,color=nscol,transform=data_crs,zorder=1)
    ax.plot(bso.LONGITUDE,bso.LATITUDE,linewidth=lwi,color=nscol,transform=data_crs,zorder=1)
    ax.plot(osnap.LONGITUDE,osnap.LATITUDE,linewidth=lwi,color=nscol,transform=data_crs,zorder=1)
    ax.set_facecolor('#c6dbef')
    savefig(figdir+'NordicSea_Map_newcol.png',bbox_inches='tight')
    savefig(figdir+'NordicSea_Map_newcol.pdf',bbox_inches='tight')

map_setup()

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
