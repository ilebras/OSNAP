from firstfuncs_1618 import *

dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF3_2016recovery_dailymerged.nc')

dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF3_2018recovery_dailymerged.nc')


def plot_overview(dat1,dat2,savename):
    f,[ax1,ax2,ax3,ax4]=subplots(4,1,figsize=(20,20),sharex=True)
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
    ax4.plot(dat1.TIME,dat1.PDEN)
    ax4.plot(dat2.TIME,dat2.PDEN)
    ax4.set_ylabel('pot. density anomaly')
    savefig(figdir+'merging_overview/'+savename+'.png')

plot_overview(dat16,dat18,'CF3_overview_1618')

def TSplot(dat1,dat2,savename):
    figure(figsize=(8,6))
    plot(dat1.PSAL,dat1.PTMP,'o',alpha=0.3)
    plot(dat2.PSAL,dat2.PTMP,'o',alpha=0.3)
    xlabel('practical salinity []')
    ylabel('pot. temperature [$^\circ$C]')
    savefig(figdir+'merging_overview/'+savename+'.png')


TSplot(dat16,dat18,'CF3_TS_1618')
