from aux_funcs import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'

dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF2_2016recovery_dailymerged.nc')

dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF2_2018recovery_dailymerged.nc')


def plot_overview(dat1,dat2):
    f,[ax1,ax2,ax3]=subplots(3,1,figsize=(12,15),sharex=True)
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

plot_overview(dat16,dat18)
savefig(figdir+'merging_overview/CF2_overview_1618.png')


def TSplot(dat1,dat2):
    figure(figsize=(9,8))
    plot(dat1.PSAL,dat1.PTMP,'o',alpha=0.3)
    plot(dat2.PSAL,dat2.PTMP,'o',alpha=0.3)
    xlabel('practical salinity []')
    ylabel('pot. temperature [$^\circ$C]')


TSplot(dat16,dat18)
savefig(figdir+'merging_overview/CF2_TS_1618.png')

#############################################################################
###### For this case, I'm going to use the first deployment as a baseline
### See if it makes sense to use the same technique as CF1, but informed by first deployment somehow...
## reminder: constant t offset, monthly s offset
#############################################################################

def plot_diff_manyways(dat,axx,thediff,colch,labch):
    axx.plot(dat.TIME,thediff,label='',alpha=0.5,color=colch)
    axx.plot(dat.resample(TIME='1M').mean(dim='TIME').TIME,thediff.resample(TIME='1M').mean(dim='TIME'),'o-',color=colch,linewidth=3,label=labch)
    axx.axhline(mean(thediff),color=colch)

# Plot the difference between 50db instrument temp and sal with the two other instruments.
def plot_saltmp_diff():
    f,[ax1,ax2]=subplots(2,1,figsize=(12,10),sharex=True)
    plot_diff_manyways(dat16,ax1,dat16.PSAL.sel(DEPTH=50)-dat16.PSAL.sel(DEPTH=100),'C0','50m-100m')
    plot_diff_manyways(dat16,ax1,dat16.PSAL.sel(DEPTH=50)-dat16.PSAL.sel(DEPTH=200),'C1','50m-200m')
    plot_diff_manyways(dat16,ax1,dat16.PSAL.sel(DEPTH=100)-dat16.PSAL.sel(DEPTH=200),'C2','100m-200m')
    plot_diff_manyways(dat18,ax1,dat18.PSAL.sel(DEPTH=100)-dat18.PSAL.sel(DEPTH=200),'C3','100m-200m')
    ax1.axhline(0,color='k')
    ax1.legend()
    ax1.set_ylabel('salinity difference')

    plot_diff_manyways(dat16,ax2,dat16.PTMP.sel(DEPTH=50)-dat16.PTMP.sel(DEPTH=100),'C0','50m-100m')
    plot_diff_manyways(dat16,ax2,dat16.PTMP.sel(DEPTH=50)-dat16.PTMP.sel(DEPTH=200),'C1','50m-200m')
    plot_diff_manyways(dat16,ax2,dat16.PTMP.sel(DEPTH=100)-dat16.PTMP.sel(DEPTH=200),'C2','100m-200m')
    plot_diff_manyways(dat18,ax2,dat18.PTMP.sel(DEPTH=100)-dat18.PTMP.sel(DEPTH=200),'C3','100m-200m')
    ax2.axhline(0,color='k')
    ax2.set_ylabel('temperature difference')

plot_saltmp_diff()

#############################################################################
###### No apparent seasonality, going to go with constant for both I think.
##### Kind of interesting in its own right...
#############################################################################
#### Refer to CF1 script for guidance on next steps if necessary
