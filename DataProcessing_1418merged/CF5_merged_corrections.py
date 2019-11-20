from firstfuncs_1618 import *

dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF5_2016recovery_dailymerged.nc')

dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF5_2018recovery_dailymerged.nc')


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

plot_overview(dat16,dat18,'CF5_overview_1618')

def TSplot(dat1,dat2,savename):
    figure(figsize=(8,6))
    plot(dat1.PSAL,dat1.PTMP,'o',alpha=0.3)
    plot(dat2.PSAL,dat2.PTMP,'o',alpha=0.3)
    xlabel('practical salinity []')
    ylabel('pot. temperature [$^\circ$C]')
    savefig(figdir+'merging_overview/'+savename+'.png')


TSplot(dat16,dat18,'CF5_TS_1618')

#############################################################################
###### Make a few comparisons in light of what was done for first deployment
#### Compare salinity and density of instruments directly above one another in both deployments
#############################################################################

dat16.DEPTH

dat18.DEPTH
def zoominbot(d1,d2,dlist):
    f,[ax1,ax2,ax3,ax4]=subplots(4,1,figsize=(20,12),sharex=True)
    for dd in dlist:
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

zoominbot(dat16,dat18,[750,1000,1300])
savefig(figdir+'merging_overview/CF5_precorr_bottom_overview.png')

#Troubling drift in the 1000m instrument -- take our that data, will interpolate over it.
dat18_corr=dat18.copy()
dat18_corr=dat18_corr.where(dat18_corr.DEPTH!=1000).dropna(dim='DEPTH')


zoominbot(dat16,dat18,[250,500,750])
savefig(figdir+'merging_overview/CF5_precorr_mid_overview.png')
# 500m instrument has a strange salinity anomaly that I will remove.
dat18_corr.PSAL.sel(DEPTH=500).plot(figsize=(14,3))
dat18_corr.PSAL.sel(DEPTH=500).where(dat18_corr.sel(DEPTH=500).PSAL>=34.805).plot()
dat18_corr.PSAL[dat18_corr.PSAL[:,-2]<=34.805,-2]=NaN

dat18_corr.PSAL.sel(DEPTH=500).plot(figsize=(14,3))

zoominbot(dat16,dat18,[100,250])
savefig(figdir+'merging_overview/CF5_precorr_top_overview.png')


dat18_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF5_mcat_corr_2018recovery_daily.nc','w',format='netCDF4')
