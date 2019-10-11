from firstfuncs_1618 import *

CF5=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF5_2018recovery_dailymerged.nc')
CF6=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF6_2018recovery_dailymerged.nc')
CF7=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF7_2018recovery_dailymerged.nc')
M1=xr.open_dataset(datadir+'OSNAP2018recovery/M1_netcdf/M1_mcat_0316_daily.nc')
CF5.DEPTH

def botplot(var,savename):
    CF6[var].sel(DEPTH=1800).plot(figsize=(14,4),label='CF6, 1800m')
    CF7[var].sel(DEPTH=1900).plot(label='CF7, 1900m')
    legend()
    title('')
    savefig(figdir+'merging_overview/'+savename+'.png',bbox_inches='tight')


botplot('PSAL','BotDriftComp_CF6CF7_PSAL')

botplot('PTMP','BotDriftComp_CF6CF7_PTMP')

botplot('PDEN','BotDriftComp_CF6CF7_PDEN')
