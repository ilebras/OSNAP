from firstfuncs_1618 import *

dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF4_2016recovery_dailymerged.nc')


dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF4_2018recovery_dailymerged.nc')


def plot_overview(dat1,dat2,savename):
    f,[ax1,ax2,ax3]=subplots(3,1,figsize=(20,15),sharex=True)
    ax1.plot(dat1.TIME,dat1.PRES)
    ax1.plot(dat2.TIME,dat2.PRES)
    ax1.invert_yaxis()
    ax1.set_ylabel('pressure [db]')
    ax2.plot(dat1.TIME,dat1.PTMP)
    ax2.plot(dat2.TIME,dat2.PTMP)
    ax2.set_ylabel('pot. temperature [$^\circ$C]')
    ax3.plot(dat1.TIME,dat1.PSAL)
    ax3.plot(dat2.TIME,dat2.PSAL)
    ax3.set_ylabel('practical salinity []')
    savefig(figdir+'merging_overview/'+savename+'.png')

plot_overview(dat16,dat18,'CF4_overview_1618')

def TSplot(dat1,dat2,savename):
    figure(figsize=(9,8))
    plot(dat1.PSAL,dat1.PTMP,'o',alpha=0.3)
    plot(dat2.PSAL,dat2.PTMP,'o',alpha=0.3)
    xlabel('practical salinity []')
    ylabel('pot. temperature [$^\circ$C]')
    savefig(figdir+'merging_overview/'+savename+'.png')


TSplot(dat16,dat18,'CF4_TS_1618')

#############################################################################
###### Make a few comparisons in light of what was done for first deployment
#### Compare salinity and density of instruments directly above one another in both deployments
#############################################################################

def comp_densal(d1,d2,tit,savename):
    dat1=dat16
    dat2=dat18
    if d2==400:
        d2_18=380
    else:
        d2_18=d2
    f,[ax1,ax2]=subplots(2,1,figsize=(10,5),sharex=True)
    ax1.plot(dat1.TIME,dat1.PSAL.sel(DEPTH=d1)-dat1.PSAL.sel(DEPTH=d2))
    ax1.plot(dat2.TIME,dat2.PSAL.sel(DEPTH=d1)-dat2.PSAL.sel(DEPTH=d2_18))
    ax1.axhline(0,color='k')
    ax1.set_title(tit)
    ax1.set_ylabel('salinity difference')
    ax2.plot(dat1.TIME,dat1.PDEN.sel(DEPTH=d1)-dat1.PDEN.sel(DEPTH=d2))
    ax2.plot(dat2.TIME,dat2.PDEN.sel(DEPTH=d1)-dat2.PDEN.sel(DEPTH=d2_18))
    ax2.axhline(0,color='k')
    ax2.set_ylabel('pot. density difference')
    savefig(figdir+'merging_overview/CF4_'+savename+'.png')

comp_densal(50,100,'50m-100m instrument','saldendiff_50-100m')

comp_densal(350,400,'350m-400m instrument','saldendiff_350-400m')

comp_densal(200,400,'200m-400m instrument','saldendiff_200-400m')

def plotonevar(seldat,selvar,dpth1,dpth2):
    figure(figsize=(15,4))
    plot(seldat.TIME,seldat[selvar].sel(DEPTH=dpth1));
    plot(seldat.TIME,seldat[selvar].sel(DEPTH=dpth2));


plotonevar(dat16,'PSAL',200,350)

plotonevar(dat16,'PSAL',350,400)

plotonevar(dat16,'PSAL',200,400)
plotonevar(dat16,'PRES',350,400)

plotonevar(dat18,'PDEN',200,350)


#############################################################################
#############################################################################
###### Conclusion after poking around: 350m instrument which is on the ADCP frame is the issue: remove in both years
###### Still some inversions between 200 and 400, but unclear what its coming from
###### 350 instrument has strange patterns in it
###### Note that all instruments on this mooring were hairy on 2018 recovery
#############################################################################
#############################################################################

dat16_corr=dat16.where(dat16.DEPTH!=350).dropna(dim='DEPTH')
dat18_corr=dat18.where(dat18.DEPTH!=350).dropna(dim='DEPTH')

dat16_corr.to_netcdf(datadir+'OSNAP2016recovery/mcat_nc/CF4_mcat_corr_2016recovery_daily.nc','w',format='netCDF4')
dat18_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF4_mcat_corr_2018recovery_daily.nc','w',format='netCDF4')
