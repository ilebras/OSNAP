from firstfuncs_1618 import *

lmask=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/ERA5/era5_landmask.nc')
lind=where(lmask.longitude>180)[0][0]
lmask['longitude']=hstack((lmask['longitude'][:lind],lmask['longitude'][lind:]-360))
lmask['lsm_better']=(('latitude','longitude'),lmask.lsm.values[0,:,:])
dat=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/ERA5/era5_surface_fluxes_monthly_averaged_for_NordicSeas.nc')

dat.msr_0001.mean(dim='time').plot()
dat=dat.where(lmask.lsm_better==0)

dat

#load NORESM boundaries and extract data within box
osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')

bnd_lon=hstack((fs.LONGITUDE.values,bso.LONGITUDE.values,ns.LONGITUDE.values[::-1],osnap.LONGITUDE.values[::-1],-45,-40,fs.LONGITUDE.values[0]))
bnd_lat=hstack((fs.LATITUDE.values,bso.LATITUDE.values,ns.LATITUDE.values[::-1],osnap.LATITUDE.values[::-1],60,75,fs.LATITUDE.values[0]))

dat.tp_0001.mean(dim='time').plot()
plot(bnd_lon,bnd_lat,linewidth=3,color='k')

import regionmask

bnd_tpl=array(tuple((zip(bnd_lon,bnd_lat))))


mybox = regionmask.Regions([bnd_tpl], name='mybox')
ax = mybox.plot()


#create a mask
mymask = mybox.mask(dat.longitude.values,dat.latitude.values)

mymask.where(mymask==0).plot()

mymask=mymask.rename({'lon':'longitude','lat':'latitude'})
dat.tp_0001.mean(dim='time').where(mymask==0).plot()
cut_dat=dat.where(mymask==0)

cut_dat.tp_0001.mean(dim='time').plot()



lonmat,latmat=meshgrid(dat.longitude,dat.latitude)
shape(lonmat-0.25/2)
shape(lonmat[:,-1]+0.25/2)
loni=vstack((lonmat.T-0.25/2,lonmat[:,-1].T+0.25/2)).T
lathold=vstack((latmat.T,latmat[:,-1].T)).T
lati=vstack((latmat+0.25/2,latmat[-1,:]+0.25/2))
lonhold=vstack((latmat,latmat[-1,:]))
plot(lonmat[:10,:10],latmat[:10,:10],'.');
plot(loni[:10,:10],lati[:10,:10],'o');


latdist=gsw.distance(lonhold,lati,axis=0)
londist=gsw.distance(loni,lathold)
londist
plot(latdist)
cut_dat['area']=(('latitude','longitude'),londist*latdist)
cut_dat['area'].plot()
for var in cut_dat:
    if '0005' in var:
        cut_dat=cut_dat.drop(var)

cut_dat

for var in cut_dat:
    if '0001' in var:
        if ('hf' in var) | ('wrf' in var):
            cut_dat[var[:-5]+'_int']=(cut_dat[var]*cut_dat['area']).sum(dim='longitude').sum(dim='latitude')/1e12
        elif ('tp_' in var) | ('e_' in var ):
            cut_dat[var[:-5]+'_int']=(cut_dat[var]*cut_dat['area']).sum(dim='longitude').sum(dim='latitude')/24/60/60/1e6


#From ERA5 documentation:
#The accumulations in monthly means of daily means (stream=moda/edmo) have been scaled to have units that include "per day",
#so for accumulations in these streams:
#The hydrological parameters are in units of "m of water per day"
#The energy (turbulent and radiative) and momentum fluxes should be divided by 86400 seconds (24 hours) to convert to the commonly used units of Wm-2 and Nm-2, respectively.
# seems not to apply to heat flux, and confused about radiative fluxes, but they are not what I want so ignoring for now...
(cut_dat.mslhf_int+cut_dat.msshf_int+cut_dat.msnlwrf_int+cut_dat.msnswrf_int).plot()


(cut_dat.tp_int).plot()
cut_dat.e_int.plot()
((cut_dat.tp_int-cut_dat.e_int)).plot()

cut_dat.to_netcdf(datadir+'aux_data/Air-Sea_FW/ERA5/ERA5_ephf_cutint_ilebras.nc','w')

# 'evaporation', 'mean_convective_precipitation_rate', 'mean_evaporation_rate',
# 'mean_large_scale_precipitation_rate', 'mean_potential_evaporation_rate', 'mean_runoff_rate',
# 'mean_snow_evaporation_rate', 'mean_snowfall_rate', 'mean_surface_latent_heat_flux',
# 'mean_surface_net_long_wave_radiation_flux', 'mean_surface_net_short_wave_radiation_flux', 'mean_surface_runoff_rate',
# 'mean_surface_sensible_heat_flux', 'mean_total_precipitation_rate', 'runoff',
# 'total_precipitation',
