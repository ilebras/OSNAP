from firstfuncs_1618 import *

###### Direct to the velocity data
def load_vel(moornum):
    if moornum==8:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/M1_netcdf/M1_vel_2016recovery_hourly.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/M1_netcdf/M1_vel_2018recovery_hourly.nc')
    else:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Hourly_netcdf/CF'+str(moornum)+'_vel_2016recovery_hourly_dropadcpbin.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Hourly_netcdf/CF'+str(moornum)+'_vel_2018recovery_hourly_dropadcpbin.nc')

    return dat16,dat18

# Load all and combine into a single xarray
moornum=5
dat16,dat18=load_vel(moornum)
dat=xr.concat([dat16,dat18],dim='TIME').drop('LATITUDE').drop('LONGITUDE').resample(TIME='1H').mean(dim='TIME')
dat=dat.assign_coords(distance=(distvec[moornum-1]))
for moornum in range(6,9):
    dat16,dat18=load_vel(moornum)
    dat_tmp=xr.concat([dat16,dat18],dim='TIME').drop('LATITUDE').drop('LONGITUDE').resample(TIME='1H').mean(dim='TIME')
    dat_tmp=dat_tmp.assign_coords(distance=(distvec[moornum-1]))
    if moornum==8:
        dat_tmp=dat_tmp.drop('WCUR')
    dat=xr.concat([dat,dat_tmp],dim='distance')

dat=dat.transpose('DEPTH','distance','TIME')

dat.UCUR[:,1,:].plot()


dat['across track velocity']=dat.UCUR*cos(theta)+dat.VCUR*sin(theta)
dat['along track velocity']=-dat.UCUR*sin(theta)+dat.VCUR*cos(theta)
DVmat=-gsw.z_from_p(dat['PRES'],60)

datadic={}
datadic['U']=dat['along track velocity'].values
datadic['V']=dat['across track velocity'].values
datadic['DV']=DVmat

for kk in datadic:
    print(kk)
    print(shape(datadic[kk]))

io.savemat(open(datadir+'OSNAP_CFgridded_2014-2018/hourly/CF_M_hourlyvel_2014-2018.mat','wb'),datadic)
