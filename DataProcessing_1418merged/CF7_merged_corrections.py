from firstfuncs_1618 import *

dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF7_mcat_2016recovery_daily.nc')

dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF7_mcat_2018recovery_daily.nc')


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
    savefig(figdir+'merging_overview_mcats/'+savename+'.png')

plot_overview(dat16,dat18,'CF7_overview_1618')

def TSplot(dat1,dat2,savename):
    figure(figsize=(8,6))
    plot(dat1.PSAL,dat1.PTMP,'o',alpha=0.3)
    plot(dat2.PSAL,dat2.PTMP,'o',alpha=0.3)
    xlabel('practical salinity []')
    ylabel('pot. temperature [$^\circ$C]')
    savefig(figdir+'merging_overview_mcats/'+savename+'.png')


TSplot(dat16,dat18,'CF7_TS_1618')

#############################################################################
###### Make a few comparisons in light of what was done for first deployment
#### Compare salinity and density of instruments directly above one another in both deployments
#############################################################################

dat16.DEPTH
dat18.DEPTH
def zoominbot(d1,d2):
    f,[ax1,ax2,ax3,ax4]=subplots(4,1,figsize=(20,12),sharex=True)
    for dd in [1000,1500,1900]:
        ax1.plot(d1.TIME,d1.PRES.sel(DEPTH=dd));
        ax1.plot(d2.TIME,d2.PRES.sel(DEPTH=dd));
        ax1.set_ylabel('pressure [db]')
        ax2.plot(d1.TIME,d1.PSAL.sel(DEPTH=dd));
        ax2.plot(d2.TIME,d2.PSAL.sel(DEPTH=dd));
        ax2.set_ylabel('salinity []')
        ax3.plot(d1.TIME,d1.PTMP.sel(DEPTH=dd));
        ax3.plot(d2.TIME,d2.PTMP.sel(DEPTH=dd));
        ax3.set_ylabel('pot. temperature [$^\circ$C]')
        ax4.plot(d1.TIME,d1.PDEN.sel(DEPTH=dd));
        ax4.plot(d2.TIME,d2.PDEN.sel(DEPTH=dd));
        ax4.set_ylabel('pot. density anomaly')

zoominbot(dat16,dat18)
savefig(figdir+'merging_overview_mcats/CF7_precorr_overview.png')

#############################################
### A closer look at mid-depth salinity time series
#############################################

dat16.DEPTH

def sal_checkout():
    figure(figsize=(16,12))
    subplot(311)
    dat16.sel(DEPTH=750).PSAL.plot()
    dat18.sel(DEPTH=750).PSAL.plot()
    gca().set_xticklabels('')
    xlabel('')
    ylim(34.88,34.98)
    grid('on')
    subplot(312)
    dat16.sel(DEPTH=1000).PSAL.plot()
    dat18.sel(DEPTH=1000).PSAL.plot()
    gca().set_xticklabels('')
    xlabel('')
    ylim(34.88,34.98)
    grid('on')
    subplot(313)
    dat16.sel(DEPTH=1500).PSAL.plot()
    dat18.sel(DEPTH=1500).PSAL.plot()
    ylim(34.88,34.98)
    grid('on')
    savefig(figdir+'merging_overview_mcats/CF7_middepth_sals_sep.png',bbox_inches='tight')

    figure(figsize=(16,4))
    dat16.sel(DEPTH=250).PSAL.plot()
    dat18.sel(DEPTH=250).PSAL.plot()
    dat18.sel(DEPTH=500).PSAL.plot()

    dat16.sel(DEPTH=750).PSAL.plot()
    dat18.sel(DEPTH=750).PSAL.plot()

    dat16.sel(DEPTH=1000).PSAL.plot()
    dat18.sel(DEPTH=1000).PSAL.plot()

    dat16.sel(DEPTH=1500).PSAL.plot()
    dat18.sel(DEPTH=1500).PSAL.plot()
    ylim(34.88,34.98)
    title('')
    grid('on')
    savefig(figdir+'merging_overview_mcats/CF7_middepth_sals_tog.png',bbox_inches='tight')

    figure(figsize=(16,4))

    dat16.sel(DEPTH=250).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    dat18.sel(DEPTH=250).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    # dat16.sel(DEPTH=500).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    dat18.sel(DEPTH=500).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    dat16.sel(DEPTH=750).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    dat18.sel(DEPTH=750).PSAL.resample(TIME='M').mean(dim='TIME').plot()

    dat16.sel(DEPTH=1000).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    dat18.sel(DEPTH=1000).PSAL.resample(TIME='M').mean(dim='TIME').plot()

    dat16.sel(DEPTH=1500).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    dat18.sel(DEPTH=1500).PSAL.resample(TIME='M').mean(dim='TIME').plot()
    ylim(34.88,34.98)
    title('')
    grid('on')
    savefig(figdir+'merging_overview_mcats/CF7_middepth_sals_smo.png',bbox_inches='tight')

    tstd='2W'
    figure(figsize=(16,4))


    dat16.sel(DEPTH=750).PSAL.resample(TIME=tstd).std(dim='TIME').plot()
    dat18.sel(DEPTH=750).PSAL.resample(TIME=tstd).std(dim='TIME').plot()

    dat16.sel(DEPTH=1000).PSAL.resample(TIME=tstd).std(dim='TIME').plot()
    dat18.sel(DEPTH=1000).PSAL.resample(TIME=tstd).std(dim='TIME').plot()

    dat16.sel(DEPTH=1500).PSAL.resample(TIME=tstd).std(dim='TIME').plot()
    dat18.sel(DEPTH=1500).PSAL.resample(TIME=tstd).std(dim='TIME').plot()
    # ylim(34.88,34.98)
    title('')
    grid('on')
    savefig(figdir+'merging_overview_mcats/CF7_middepth_sals_std.png',bbox_inches='tight')



sal_checkout()



dat18_corr=dat18.copy()
dat18_corr['PSAL'][:,-1]=dat18['PSAL'].sel(DEPTH=1500)

def add_SA_CT_PT(xray):
    SA_out=gsw.SA_from_SP(xray['PSAL'],xray['PRES'],xray.LONGITUDE,xray.LATITUDE)
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],xray['PRES'])
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['ASAL']=(('TIME','DEPTH'),SA_out)
    xray['PTMP']=(('TIME','DEPTH'),PT_out)
    xray['CTMP']=(('TIME','DEPTH'),CT_out)
    xray['PDEN']=(('TIME','DEPTH'),PD_out)

add_SA_CT_PT(dat18_corr)


zoominbot(dat16,dat18_corr)
savefig(figdir+'merging_overview/CF7_deepcorr_overview.png')

# This will do, Pig

plot(dat18_corr.PSAL,dat18_corr.PTMP,'.')
savefig(figdir+'merging_overview/CF7_deepcorr_TS.png')

dat18_corr.DEPTH
import calendar
def toTime(d):
  return [calendar.timegm(dd.timetuple()) for dd in d]

def np64ToDatetime(DA):
  return [datetime.datetime.utcfromtimestamp((dd-np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for dd in DA]



tvec18=toTime(np64ToDatetime(dat18.TIME.values))

def removeline(tvec,svec):
    linefunc=poly1d(polyfit(tvec,svec, 1))
    liney=linefunc(tvec)
    corrected=svec-liney+liney[0]
    return liney,corrected

na,psalcorr_100=removeline(tvec18,dat18.PSAL.sel(DEPTH=100))

dat16.DEPTH
dat18.DEPTH

def comp_salcorr100():
    figure(figsize=(20,4))
    for dd in [50,100,250]:
        dat16.PSAL.sel(DEPTH=dd).plot(label=dd)
        dat18.PSAL.sel(DEPTH=dd).plot(label=dd)
    # psalcorr_100.plot(label='100 corrected')
    ylim(34.6,35.1)
    title('')
    ylabel('salinity')
    legend()


comp_salcorr100()
savefig(figdir+'merging_overview/CF7_100m_zoomin.png')

#############################################################################
#############################################################################
###### Since we have the 50m instrument the whole time, just going to nan the 100m instrument data out
###### Note that the temperature is good, so might want to refer back to that.
###### Also going to take out the 50m instrument on first deployment - might be causing issues.
#############################################################################
#############################################################################
dat16_corr=dat16.where(dat16.DEPTH!=50).dropna(dim='DEPTH')
dat18_corr=dat18_corr.where(dat18_corr.DEPTH!=100).dropna(dim='DEPTH')

plot(dat18_corr.PSAL,dat18_corr.PTMP,'.')
savefig(figdir+'merging_overview/CF7_corrfinal_TS.png')

dat16_corr.to_netcdf(datadir+'OSNAP2016recovery/mcat_nc/CF7_mcat_corr_2016recovery_daily.nc','w',format='netCDF4')
dat18_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF7_mcat_corr_2018recovery_daily.nc','w',format='netCDF4')
