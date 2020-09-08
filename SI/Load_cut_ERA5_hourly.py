from aux_funcs_2020 import *

datlist=sort(glob.glob(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/hourly_downloaded_200618/adaptor*'))

datlist

dat=xr.open_dataset(datlist[0])
for dd in datlist[1:]:
    dattmp=xr.open_dataset(dd)
    dat=xr.concat([dat,dattmp],dim='time')


cut_dat=dat.sel(latitude=60).sel(longitude=slice(-43,-41)).sortby('time')

cut_dat['tau']=cut_dat.metss*cos(theta)+cut_dat.mntss*sin(theta)
cut_dat['hf']=cut_dat['msnlwrf']+cut_dat['msnswrf']+cut_dat['msshf']+cut_dat['mslhf']

def plotvar(var):
    cut_dat[var].mean(dim='longitude').resample(time='1D').mean(dim='time').plot()

plotvar('hf')
plotvar('tau')

cut_dat.to_netcdf(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/hourly_downloaded_200618/ERA5_hourly_SI_ilebras.nc','w',format='netCDF4')
