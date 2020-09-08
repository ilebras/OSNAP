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
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Hourly_netcdf/CF'+str(moornum)+'_mcat_2018recovery_hourlymerged.nc')

    return dat16,dat18

##########################################################################################################
################################   Start with vel ########################################################
#########################################################################################################
moornum=4
dat16,dat18=load_vel(moornum)

dat=xr.concat([dat16,dat18],dim='TIME').drop('LATITUDE').drop('LONGITUDE')

# dat=dat.assign_coords(distance=(distvec[moornum-1]))
# for moornum in range(6,9):
#     dat16,dat18=load_vel(moornum)
#     dat_tmp=xr.concat([dat16,dat18],dim='TIME').drop('LATITUDE').drop('LONGITUDE')
#     dat_tmp=dat_tmp.assign_coords(distance=(distvec[moornum-1]))
#     if moornum==8:
#         dat_tmp=dat_tmp.drop('WCUR')
#     dat=xr.concat([dat,dat_tmp],dim='distance')

dat

# dat=dat.transpose('DEPTH','distance','TIME')

# dat['across track velocity']=dat.UCUR*cos(theta)+dat.VCUR*sin(theta)
# dat['along track velocity']=-dat.UCUR*sin(theta)+dat.VCUR*cos(theta)


##########################################################################################################
################################   Add mcat       ########################################################
#########################################################################################################
moornum=4
d16,d18=load_mcat(moornum)

help(resample)
d16_h=d16.resample(TIME='1H').nearest()
d16_h=d16_h.rename({'PRES':'PRES_CTD'})
dat16_h=dat16.resample(TIME='1H').nearest()
dat16_h=dat16_h.rename({'PRES':'PRES_VEL'})
bothdat16=xr.merge([d16_h,dat16_h])
bothdat16

bothdat16.PRES_CTD.plot()
datadir
bothdat16.to_netcdf(datadir+'OSNAP2016recovery/Hourly_netcdf/CF4_mcat_vel_2016recovery_forMonica_200515.nc','w')

d=xr.concat([d16.drop('LATITUDE').drop('LONGITUDE'),d18.drop('LATITUDE').drop('LONGITUDE')],dim='TIME')#
d=d.assign_coords(distance=(distvec[moornum-1]))
for moornum in range(6,9):
    d16,d18=load_mcat(moornum)
    d_tmp=xr.concat([d16.drop('LATITUDE').drop('LONGITUDE'),d18.drop('LATITUDE').drop('LONGITUDE')],dim='TIME')
    d_tmp=d_tmp.assign_coords(distance=(distvec[moornum-1]))
    d=xr.concat([d,d_tmp],dim='distance')
d=d.transpose('DEPTH','distance','TIME')

d['PSAL'][:,2,15:70].plot()
d

d.TIME[15:70]

Dmat=-gsw.z_from_p(d['PRES'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values,60)
DVmat=-gsw.z_from_p(dat['PRES'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values,60)

datadic={}
datadic['U']=dat['along track velocity'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values
datadic['V']=dat['across track velocity'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values
datadic['DV']=DVmat
datadic['D']=Dmat
datadic['T']=d['PTMP'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values
datadic['S']=d['PSAL'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values
datadic['R']=d['PDEN'].sel(TIME=slice(d.TIME[0],dat.TIME[-1])).values

for kk in datadic:
    print(kk)
    print(shape(datadic[kk]))

io.savemat(open(datadir+'OSNAP_CFgridded_2014-2018/hourly/CF_M_hourly_2014-2018_200130.mat','wb'),datadic)
