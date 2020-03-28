from firstfuncs_1618 import *

dat=xr.open_dataset(datadir+'aux_data/Air-Sea_Feili/ERA5/era5_surface_fluxes_monthly_averaged_for_NordicSeas.nc')

dat.e_0001.mean(dim='time').plot()

dat.e_0001.sum(dim='latitude').sum(dim='longitude').plot()

dat.tp_0001.sum(dim='latitude').sum(dim='longitude').plot()


dat.tp_0001.
