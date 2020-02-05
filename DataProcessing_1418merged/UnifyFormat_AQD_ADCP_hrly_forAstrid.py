from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'

def quickpresplot(xray):
    figure(figsize=(14,6))
    plot(xray.TIME,xray.PRES,'.');

# Not going to filter these, so that Astrid can apply same exact filtering!

################## CF moorings
for yr in ['2016','2018']:
    for moornum in range(7,8):
        aqdir=glob.glob(datadir+'OSNAP'+yr+'recovery/AQD_Data_CF/*CF'+str(moornum)+'*')
        for ii,qq in enumerate(aqdir):
            if ii==0:
                aqd=xr.open_dataset(aqdir[0]).drop('MAGVAR')
                aqd=aqd.resample(TIME='1H').mean()
            else:
                aqd_tmp=xr.open_dataset(qq).drop('MAGVAR')
                aqd_tmp=aqd_tmp.resample(TIME='1H').mean()
                aqd=xr.concat([aqd,aqd_tmp],dim='DEPTH')
        # aqd_sm=lowpassfilt(aqd,0.5)
        if yr=='2016':
            adcpdir=glob.glob(datadir+'OSNAP'+yr+'recovery/ADCP_Data_CF/*cf'+str(moornum)+'*.nc')
            adcp=xr.open_dataset(adcpdir[0]).rename({'BINDEPTH':'DEPTH'}).drop('PITCH').drop('ROLL').drop('HEADING').drop('MAGVAR')
            adcp=adcp.rename({'PRES':'INSTPRES'})
            adcp['PRES']=adcp.INSTPRES.values*adcp.VCUR/adcp.VCUR-hstack((0,abs(diff(adcp.DEPTH)))).cumsum()
            adcp=adcp.drop('INSTPRES').drop('INSTRDEPTH')
        elif yr=='2018':
            adcpdir=glob.glob(datadir+'OSNAP'+yr+'recovery/ADCP_Data_CF/*CF'+str(moornum)+'*.nc')
            adcp=xr.open_dataset(adcpdir[0]).drop('DEPTH').rename({'BINDEPTH':'DEPTH'}).drop('INSTRDEPTH').drop('PITCH').drop('ROLL').drop('HEADING').drop('MAGVAR').drop('PREST')
        # adcp=adcp.where(adcp.DEPTH<adcp.DEPTH.max())
        adcp=adcp.resample(TIME='1H').nearest()
        # adcp_sm=lowpassfilt(adcp,1)
        vel=xr.concat([aqd,adcp],dim='DEPTH').sortby('DEPTH')
        if moornum==7:
            vel.UCUR[:,vel.DEPTH==500]=NaN
        quickpresplot(vel)
        title('CF'+str(moornum)+', '+yr)
        vel.to_netcdf(datadir+'OSNAP'+yr+'recovery/Hourly_netcdf/CF'+str(moornum)+'_vel_'+yr+'recovery_hourly_dropadcpbin.nc','w',format='netCDF4')

############################# M1
# What does M1 2016-2018 data look like? Might be what I'm trying to match.
m1dir18=glob.glob(datadir+'OSNAP2018recovery/M1_netcdf/*Nortek*')
m1_18_aqd=xr.open_dataset(m1dir18[0]).resample(TIME='1H').mean()

m1dir16=glob.glob(datadir+'OSNAP2016recovery/M1_netcdf/*Nortek*')
m1_14=xr.open_dataset(m1dir16[0]).resample(TIME='1H').mean()
m1_15=xr.open_dataset(m1dir16[1]).resample(TIME='1H').mean()

m1_16_aqd=xr.concat([m1_14,m1_15],dim='DEPTH').sortby('DEPTH').resample(TIME='1H').mean()

m1dir18_adcp=glob.glob(datadir+'OSNAP2018recovery/M1_netcdf/*ADCP*')
m1_18_adcp=xr.open_dataset(m1dir18_adcp[0])
m1_18_adcp=m1_18_adcp.resample(TIME='1H').nearest()
m1_18_adcp=m1_18_adcp.rename({'BINDEPTH':'DEPTH'})
m1_18_adcp=m1_18_adcp.drop('INSTRDEPTH')

m1_18=xr.concat([m1_18_aqd,m1_18_adcp],dim='DEPTH')


m1dir16_adcp=glob.glob(datadir+'OSNAP2016recovery/M1_netcdf/*ADCP*')

m1_14_adcp=xr.open_dataset(m1dir16_adcp[0]).resample(TIME='1H').nearest()
m1_14_adcp=m1_14_adcp.rename({'BINDEPTH':'DEPTH'})

m1_15_adcp=xr.open_dataset(m1dir16_adcp[1]).resample(TIME='1H').nearest()
m1_15_adcp=m1_15_adcp.rename({'BINDEPTH':'DEPTH'})

m1_16_adcp=xr.concat([m1_14_adcp,m1_15_adcp],dim='DEPTH').sortby('DEPTH')
m1_16_adcp=m1_16_adcp.drop('INSTRDEPTH')

m1_16=xr.concat([m1_16_aqd,m1_16_adcp],dim='DEPTH')

m1_16.to_netcdf(datadir+'OSNAP2016recovery/M1_netcdf/M1_vel_2016recovery_hourly.nc','w',format='netCDF4')

m1_18.to_netcdf(datadir+'OSNAP2018recovery/M1_netcdf/M1_vel_2018recovery_hourly.nc','w',format='netCDF4')
