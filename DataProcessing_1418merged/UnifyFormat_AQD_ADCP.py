from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'


def quickpresplot(xray):
    figure(figsize=(14,6))
    plot(xray.TIME,xray.PRES,'.');


#here is a (2nd order butterworth) 40hour lowpass filter:
#enter the data frequency in hours
def lowpassfilt(xray,tstep):
    xray_sm=NaN*xray
    Z,X = sig.butter(2,tstep/40, output='ba')
    for var in ['PRES','UCUR','VCUR']:
        for ii,dd in enumerate(xray.DEPTH):
            nanind=~isnan(xray[var][:,ii])
            if sum(nanind)>10:
                xray_sm[var][nanind,ii]=sig.filtfilt(Z,X,xray[var][nanind,ii])
            else:
                xray_sm[var][:,ii]=NaN
    return xray_sm

################## CF moorings
for yr in ['2016','2018']:
    for moornum in range(1,7):
        aqdir=glob.glob(datadir+'OSNAP'+yr+'recovery/AQD_Data_CF/*CF'+str(moornum)+'*')
        aqdir
        for ii,qq in enumerate(aqdir):
            if ii==0:
                aqd=xr.open_dataset(aqdir[0]).drop('MAGVAR')
            else:
                aqd_tmp=xr.open_dataset(qq).drop('MAGVAR')
                aqd=xr.concat([aqd,aqd_tmp],dim='DEPTH')
        aqd_sm=lowpassfilt(aqd,0.5)
        if yr=='2016':
            adcpdir=glob.glob(datadir+'OSNAP'+yr+'recovery/ADCP_Data_CF/*cf'+str(moornum)+'*.nc')
            adcp=xr.open_dataset(adcpdir[0]).rename({'BINDEPTH':'DEPTH'}).drop('PITCH').drop('ROLL').drop('HEADING').drop('MAGVAR')
            adcp=adcp.rename({'PRES':'INSTPRES'})
            adcp['PRES']=adcp.INSTPRES.values*adcp.VCUR/adcp.VCUR-hstack((0,abs(diff(adcp.DEPTH)))).cumsum()
            adcp=adcp.drop('INSTPRES').drop('INSTRDEPTH')
        elif yr=='2018':
            adcpdir=glob.glob(datadir+'OSNAP'+yr+'recovery/ADCP_Data_CF/*CF'+str(moornum)+'*.nc')
            adcp=xr.open_dataset(adcpdir[0]).drop('DEPTH').rename({'BINDEPTH':'DEPTH'}).drop('INSTRDEPTH').drop('PITCH').drop('ROLL').drop('HEADING').drop('MAGVAR').drop('PREST')
        adcp=adcp.where(adcp.DEPTH<adcp.DEPTH.max())
        adcp_sm=lowpassfilt(adcp,1)
        vel=xr.concat([aqd_sm,adcp_sm],dim='DEPTH').sortby('DEPTH')
        vel_daily=vel.resample(TIME='1D').first()
        if moornum==7:
            vel_daily.UCUR[:,vel_daily.DEPTH==500]=NaN
        quickpresplot(vel_daily)
        title('CF'+str(moornum)+', '+yr)

        vel_daily.to_netcdf(datadir+'OSNAP'+yr+'recovery/Daily_netcdf/CF'+str(moornum)+'_vel_'+yr+'recovery_daily_dropadcpbin_drop','w',format='netCDF4')


############################# M1
# What does M1 2016-2018 data look like? Might be what I'm trying to match.
m1dir18=glob.glob(datadir+'OSNAP2018recovery/M1_netcdf/*Nortek*')
m1_18=xr.open_dataset(m1dir18[0])

m1dir16=glob.glob(datadir+'OSNAP2016recovery/M1_netcdf/*Nortek*')
m1_14=xr.open_dataset(m1dir16[0])
m1_15=xr.open_dataset(m1dir16[1])

m1_16=xr.concat([m1_14,m1_15],dim='DEPTH').sortby('DEPTH')

m1_16.PRES.plot()

m1_18_sm=lowpassfilt(m1_18,0.5)
m1_16_sm=lowpassfilt(m1_16,1)

m1_18.VCUR[:,-1].plot(figsize=(14,3))
m1_18_sm.VCUR[:,-1].plot()

m1_18_daily_aqd=m1_18_sm.resample(TIME='1D').first()

m1_16_daily_aqd=m1_16_sm.resample(TIME='1D').first()

m1dir18_adcp=glob.glob(datadir+'OSNAP2018recovery/M1_netcdf/*ADCP*')
m1_18_adcp=xr.open_dataset(m1dir18_adcp[0])
m1_18_adcp=m1_18_adcp.rename({'BINDEPTH':'DEPTH'})

m1_18_adcp_sm=lowpassfilt(m1_18_adcp,1)
m1_18_daily_adcp=m1_18_adcp_sm.resample(TIME='1D').first()
m1_18_daily_adcp=m1_18_daily_adcp.drop('INSTRDEPTH')
m1_18_daily=xr.concat([m1_18_daily_aqd,m1_18_daily_adcp],dim='DEPTH')

plot(m1_18_daily.TIME,m1_18_daily.VCUR);
m1dir16_adcp=glob.glob(datadir+'OSNAP2016recovery/M1_netcdf/*ADCP*')

m1dir16_adcp
m1_14_adcp=xr.open_dataset(m1dir16_adcp[0])
m1_14_adcp=m1_14_adcp.rename({'BINDEPTH':'DEPTH'})

m1_15_adcp=xr.open_dataset(m1dir16_adcp[1])
m1_15_adcp=m1_15_adcp.rename({'BINDEPTH':'DEPTH'})

m1_16_adcp=xr.concat([m1_14_adcp,m1_15_adcp],dim='DEPTH').sortby('DEPTH')
m1_16_adcp_sm=lowpassfilt(m1_16_adcp,1)

m1_16_adcp_sm.VCUR.plot()

m1_16_daily_adcp=m1_16_adcp_sm.resample(TIME='1D').first()
m1_16_daily_adcp=m1_16_daily_adcp.drop('INSTRDEPTH')
m1_16_daily=xr.concat([m1_16_daily_aqd,m1_16_daily_adcp],dim='DEPTH')

plot(m1_16_daily.TIME,m1_16_daily.PRES)


m1_16_daily.to_netcdf(datadir+'OSNAP2016recovery/M1_netcdf/M1_vel_2016recovery_daily.nc','w',format='netCDF4')

m1_18_daily.to_netcdf(datadir+'OSNAP2018recovery/M1_netcdf/M1_vel_2018recovery_daily.nc','w',format='netCDF4')
