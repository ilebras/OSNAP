from firstfuncs_1618 import *
ooi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/OOI_HYPM_xray_fromFemke.nc')

ooi
d1='2014-10-1'
d2='2014-10-15'
plot(ooi.salinity.sel(date=slice(d1,d2)),ooi.temperature.sel(date=slice(d1,d2)),'k.');
d1='2015-3-1'
d2='2015-3-10'
plot(ooi.salinity.sel(date=slice(d1,d2)),ooi.temperature.sel(date=slice(d1,d2)),'y.');
d1='2015-4-1'
d2='2015-4-10'
plot(ooi.salinity.sel(date=slice(d1,d2)),ooi.temperature.sel(date=slice(d1,d2)),'r.');
d1='2015-6-1'
d2='2015-6-10'
plot(ooi.salinity.sel(date=slice(d1,d2)),ooi.temperature.sel(date=slice(d1,d2)),'b.');
d1='2015-9-1'
d2='2015-9-10'
plot(ooi.salinity.sel(date=slice(d1,d2)),ooi.temperature.sel(date=slice(d1,d2)),'g.');
