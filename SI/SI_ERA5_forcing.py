from aux_funcs_2020 import *

era5=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/hourly_downloaded_200618/ERA5_hourly_SI_ilebras.nc')

def plotvar(var):
    era5[var].mean(dim='longitude').plot()
    era5[var].mean(dim='longitude').resample(time='1D').mean(dim='time').plot()
    era5[var].mean(dim='longitude').resample(time='1W').mean(dim='time').plot()

plotvar('hf')
plotvar('tau')
