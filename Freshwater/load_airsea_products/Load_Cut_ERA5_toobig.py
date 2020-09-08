from firstfuncs_1618 import *

datadir

dat=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/adaptor.mars.internal-1585838857.1923296-2183-27-f2933893-82e9-4c4d-b51f-07fa5dbb6346.nc')
lind=where(dat.longitude>180)[0][0]
dat['longitude']=hstack((dat['longitude'][:lind],dat['longitude'][lind:]-360))

cut_dat=dat.sel(latitude=slice(80,55))
cut_dat=cut_dat.sortby('longitude').sel(longitude=slice(-45,25))
cut_dat=cut_dat.sel(expver=1)

cut_dat.to_netcdf(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/ERA5_2000-2020_regioncut_ilebras.nc','w')
