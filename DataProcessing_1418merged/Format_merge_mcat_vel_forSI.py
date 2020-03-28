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

def load_mcat(moornum):
    if moornum==8:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/M1_netcdf/M1_mcat_2016recovery_hourly.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/M1_netcdf/M1_mcat_2018recovery_hourly.nc')
    else:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Hourly_netcdf/CF'+str(moornum)+'_mcat_2016recovery_hourly.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Hourly_netcdf/CF'+str(moornum)+'_mcat_2018recovery_hourly.nc')

    return dat16,dat18

##########################################################################################################
################################   Start with vel ########################################################
#########################################################################################################
moornum=5
dat16,dat18=load_vel(moornum)

dat16.DEPTH

dat16.DEPTH=[int(dd) for]


dat=xr.concat([dat16,dat18],dim='TIME').drop('LATITUDE').drop('LONGITUDE')

dat=dat.assign_coords(distance=(distvec[moornum-1]))
for moornum in range(6,9):
    dat16,dat18=load_vel(moornum)
    dat_tmp=xr.concat([dat16,dat18],dim='TIME').drop('LATITUDE').drop('LONGITUDE')
    dat_tmp=dat_tmp.assign_coords(distance=(distvec[moornum-1]))
    if moornum==8:
        dat_tmp=dat_tmp.drop('WCUR')
    dat=xr.concat([dat,dat_tmp],dim='distance')

dat=dat.transpose('DEPTH','distance','TIME')

dat['across track velocity']=dat.UCUR*cos(theta)+dat.VCUR*sin(theta)
dat['along track velocity']=-dat.UCUR*sin(theta)+dat.VCUR*cos(theta)


##########################################################################################################
################################   Add mcat       ########################################################
#########################################################################################################
moornum=5
d16,d18=load_mcat(moornum)

d=xr.concat([d16.drop('LATITUDE').drop('LONGITUDE'),d18.drop('LATITUDE').drop('LONGITUDE')],dim='TIME')#
d=d.assign_coords(distance=(distvec[moornum-1]))
for moornum in range(6,9):
    d16,d18=load_mcat(moornum)
    d_tmp=xr.concat([d16.drop('LATITUDE').drop('LONGITUDE'),d18.drop('LATITUDE').drop('LONGITUDE')],dim='TIME')
    d_tmp=d_tmp.assign_coords(distance=(distvec[moornum-1]))
    d=xr.concat([d,d_tmp],dim='distance')
d=d.transpose('DEPTH','distance','TIME')

dat

d=d.drop('CNDC').drop('TEMP')
d['depth_var']=(('DEPTH','distance','TIME'),-gsw.z_from_p(d['PRES'],60))
dat['depth_var']=(('DEPTH','distance','TIME'),-gsw.z_from_p(dat['PRES'],60))

dat


datadir
d.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/hourly/CF_M_hourly_2014-2018_MCATonly.nc','w',format='netCDF4')
dat.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/hourly/CF_M_hourly_2014-2018_VELonly.nc','w',format='netCDF4')
