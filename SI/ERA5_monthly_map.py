from aux_funcs_2020 import *

dat=xr.open_dataset(glob.glob('/home/isabela/DATA/ERA5/*_monthly_windstress*')[0])

#interpolate dat to 1degree grid for visibility:
#before the land mask!
dat=dat.interp(latitude=arange(90,0,-0.5)).interp(longitude=arange(-100,30,0.75))

lmask=xr.open_dataset(glob.glob('/home/isabela/DATA/ERA5/*mask*')[0])
lind=where(lmask.longitude>180)[0][0]
lmask['longitude']=hstack((lmask['longitude'][:lind],lmask['longitude'][lind:]-360))
lmask=lmask.sortby('longitude')
lmask['lsm_better']=(('latitude','longitude'),lmask.lsm.values[0,:,:])
lmask=lmask.interp(latitude=arange(90,0,-0.5)).interp(longitude=arange(-100,30,0.75))
# lmask.lsm_better.plot()

dat=dat.where(lmask.lsm_better<0.4)

datlim=dat.sel(time=slice('2017-1-1','2017-1-31')).sel(longitude=slice(-100,20)).sel(latitude=slice(70,45))

# en4=xr.open_dataset('/home/isabela/DATA/EN4/analyses/EN.4.2.1.f.analysis.g10.201702.nc')

# ssh=xr.open_dataset('/home/isabela/Documents/projects/OSNAP/data/aux_data/SSH/CMEMS_monthly_subpolarSLA_2014-2018_dwnl200831.nc')

mdt=xr.open_dataset('/home/isabela/DATA/DTU_MDT/DTU10MDT_2min.nc')
lind=where(mdt.lon>180)[0][0]
mdt['lon_corr']=('lon',hstack((mdt['lon'][:lind],mdt['lon'][lind:]-360)))
mdt=mdt.swap_dims({'lon':'lon_corr'})
mdt_cut=mdt.mdt.sortby('lon_corr').sel(lat=slice(40,75)).sel(lon_corr=slice(-100,20))

CFlon=-42.5
CFlat=60
import cartopy.crs as ccrs
import cartopy

extent = [-60, -20, 50, 68]
central_lon, central_lat = (extent[0]+extent[1])/2, (extent[2]+extent[3])/2

def quiver_map():
    fig, ax = plt.subplots(1,1, figsize=(8,5),subplot_kw={'projection': ccrs.Orthographic(central_lon, central_lat)},constrained_layout=True)
    data_crs = ccrs.PlateCarree()
    ax.set_extent(extent)
    land_def=cartopy.feature.NaturalEarthFeature(category='physical',name='land',scale='10m',facecolor='grey')
    ax.add_feature(land_def,zorder=2)
    gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='grey',linestyle='--')
    lati=0.2
    ax.plot([CFlon-1,CFlon-1,CFlon+1,CFlon+1,CFlon-1],[CFlat-lati,CFlat+lati,CFlat+lati,CFlat-lati,CFlat-lati],'r-',transform=data_crs,zorder=3,linewidth=3)
    # col=ax.contourf(en4.lon.values,en4.lat.values,en4.salinity.isel(depth=0).isel(time=0).values,transform=data_crs,extend='both',levels=arange(33,35.3,0.05))
    col=ax.contourf(mdt_cut.lon_corr.values,mdt_cut.lat.values,mdt_cut.values,50,vmin=-0.7,vmax=0.7,transform=data_crs,cmap=cm.RdBu_r,extend='both',zorder=1)
    fax=fig.add_axes([1,0.25,0.02,0.5])
    colorbar(col,cax=fax,label='Mean dynamic topography [m]')
    ax.quiver(datlim.longitude.values,datlim.latitude.values,
              datlim.metss.mean(dim='time').values,datlim.mntss.mean(dim='time').values,
              transform=data_crs,zorder=4,headwidth=2,headlength=3,scale_units='inches',scale=1)
    ax.set_facecolor('k')
    savefig(figdir+'SI/paper_prep/Quivermap.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/Quivermap.pdf',bbox_inches='tight')

quiver_map()
