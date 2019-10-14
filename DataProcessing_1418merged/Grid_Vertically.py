from firstfuncs_1618 import *
###### Load microcat data:
moornum=4

def load_mcat(moornum):
    if moornum==1:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF'+str(moornum)+'_recon_2016recovery_dailymerged.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_recon_2018recovery_dailymerged.nc')
    elif moornum==2:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF'+str(moornum)+'_2016recovery_dailymerged.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_recon_2018recovery_dailymerged.nc')
    elif (moornum==4) | (moornum==7):
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF'+str(moornum)+'_corr_2016recovery_dailymerged.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_corr_2018recovery_dailymerged.nc')
    elif (moornum==5) | (moornum==7):
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF'+str(moornum)+'_2016recovery_dailymerged.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_corr_2018recovery_dailymerged.nc')
    elif moornum==8:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/M1_netcdf/M1_mcat_0114_0215_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/M1_netcdf/M1_mcat_0316_daily.nc')
    else:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF'+str(moornum)+'_2016recovery_dailymerged.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_2018recovery_dailymerged.nc')
    return dat16,dat18


d16,d18=load_mcat(1)
