from aux_funcs import *

dat1=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/deployment0001_GI03FLMA-RIM01-02-ADCPSL003-recovered_inst-adcp_velocity_earth_20140912T210000-20150818T100000.nc')

dat1.bin_depths[90:,0].plot()

plot(dat1.time[90:],dat1.northward_seawater_velocity[90:,0])

dat2=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/deployment0002_GI03FLMA-RIM01-02-ADCPSL003-recovered_inst-adcp_velocity_earth_20150819T003714.970000-20160716T233714.970000.nc')



dat2.bin_depths[:,0].plot()

plot(dat2.time,dat2.northward_seawater_velocity[:,0])


dat2.time[:10]

dat_500=xr.Dataset({'u':(['date'],hstack((dat1.eastward_seawater_velocity[90:,0].values,dat2.eastward_seawater_velocity[:,2].values))),
                    'v':(['date'],hstack((dat1.northward_seawater_velocity[90:,0].values,dat2.northward_seawater_velocity[:,2].values))),
                    'dpths':(['date'],hstack((dat1.bin_depths[90:,0].values,dat2.bin_depths[:,2].values))),},
                    coords={'date':hstack((dat1.time[90:].values,dat2.time.values))},)

dat_500.dpths.plot()
dat_500.u.plot()
dat_500.v.plot()
dat2.time[:10]

dat_500.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/ADCP_OOI_FLA_500m.nc','w',format='netCDF4')
