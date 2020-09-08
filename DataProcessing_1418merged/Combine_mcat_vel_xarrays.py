from firstfuncs_1618 import *

moornum=1
mcat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(moornum)+'_mcat_vertgrid_daily_2m.nc')
vel=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(moornum)+'_vel_vertgrid_daily_2m.nc')
dat=xr.merge([mcat,vel])
dat=dat.assign_coords(distance=(distvec[moornum-1]))

for moornum in range(2,9):
    mcat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(moornum)+'_mcat_vertgrid_daily_2m.nc')
    vel=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(moornum)+'_vel_vertgrid_daily_2m.nc')
    dat_tmp=xr.merge([mcat,vel])
    dat_tmp=dat_tmp.assign_coords(distance=(distvec[moornum-1]))
    dat=xr.concat([dat,dat_tmp],dim='distance')

dat.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_2m.nc','w',format='netCDF4')
