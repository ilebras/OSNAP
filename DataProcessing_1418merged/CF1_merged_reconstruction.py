from firstfuncs_1618 import *

dat16=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF1_2016recovery_dailymerged.nc')

dat16_recon=xr.open_dataset(datadir+'OSNAP2016recovery/mcat_nc/CF1_recon_2016recovery_dailymerged.nc')

dat18=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF1_2018recovery_dailymerged.nc')


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

plot_overview(dat16_recon,dat18)
savefig(figdir+'merging_overview/CF1_overview_16recon_w18.png')


def TSplot(dat1,dat2):
    figure(figsize=(9,8))
    plot(dat1.PSAL,dat1.PTMP,'o',alpha=0.3)
    plot(dat2.PSAL,dat2.PTMP,'o',alpha=0.3)
    xlabel('practical salinity []')
    ylabel('pot. temperature [$^\circ$C]')


TSplot(dat16_recon,dat18)
savefig(figdir+'merging_overview/CF1_TS_16recon_w18.png')

#############################################################################
############ The first deployment reconstruction does not look bad,
####### Going to leave as is and try something similar for the second deployment
#############################################################################

# I actually only have to reconstruct the 50dbar data.

def plot_diff_manyways(axx,thediff,colch,labch):
    axx.plot(dat18.TIME,thediff,label='',alpha=0.5,color=colch)
    axx.plot(dat18.resample(TIME='1M').mean(dim='TIME').TIME,thediff.resample(TIME='1M').mean(dim='TIME'),'o-',color=colch,linewidth=3,label=labch)
    axx.axhline(mean(thediff),color=colch)



#############################################################################
############ Reconstruct using both instruments, see how well they agree
##### Use same method as first deployment: constant t offset, monthly s offset
#############################################################################

mtime=(dat18).resample(TIME='1M').mean(dim='TIME').TIME

sal_mdiff={}
sal_mdiff['from100']=(dat18.PSAL.sel(DEPTH=50)-dat18.PSAL.sel(DEPTH=100)).resample(TIME='1M').mean(dim='TIME')
sal_mdiff['from200']=(dat18.PSAL.sel(DEPTH=50)-dat18.PSAL.sel(DEPTH=200)).resample(TIME='1M').mean(dim='TIME')

sal_mdiff_fill_100=linspace(sal_mdiff['from100'][8],sal_mdiff['from100'][0],5)
sal_mdiff_fill_200=linspace(sal_mdiff['from200'][8],sal_mdiff['from200'][0],5)

sal_mdiff_int={}
sal_mdiff_int['from100']=hstack((sal_mdiff['from100'][:9],sal_mdiff_fill_100[1:-1],sal_mdiff['from100'][:9],sal_mdiff_fill_100[1:-1],sal_mdiff['from100'][:2]))
sal_mdiff_int['from200']=hstack((sal_mdiff['from200'][:9],sal_mdiff_fill_200[1:-1],sal_mdiff['from200'][:9],sal_mdiff_fill_200[1:-1],sal_mdiff['from200'][:2]))

# Plot the difference between 50db instrument temp and sal with the two other instruments.
def plot_saltmp_diff():
    f,[ax1,ax2]=subplots(2,1,figsize=(12,10),sharex=True)
    plot_diff_manyways(ax1,dat18.PSAL.sel(DEPTH=50)-dat18.PSAL.sel(DEPTH=100),'C0','50m-100m')
    ax1.plot(mtime,sal_mdiff_int['from100'],'o-',color='C0')
    plot_diff_manyways(ax1,dat18.PSAL.sel(DEPTH=50)-dat18.PSAL.sel(DEPTH=200),'C1','50m-200m')
    ax1.plot(mtime,sal_mdiff_int['from200'],'o-',color='C1')
    plot_diff_manyways(ax1,dat18.PSAL.sel(DEPTH=100)-dat18.PSAL.sel(DEPTH=200),'C2','100m-200m')
    ax1.axhline(0,color='k')
    ax1.legend()
    ax1.set_ylabel('salinity difference')

    plot_diff_manyways(ax2,dat18.PTMP.sel(DEPTH=50)-dat18.PTMP.sel(DEPTH=100),'C0','50m-100m')
    plot_diff_manyways(ax2,dat18.PTMP.sel(DEPTH=50)-dat18.PTMP.sel(DEPTH=200),'C1','50m-200m')
    plot_diff_manyways(ax2,dat18.PTMP.sel(DEPTH=100)-dat18.PTMP.sel(DEPTH=200),'C2','100m-200m')
    ax2.axhline(0,color='k')
    ax2.set_ylabel('temperature difference')

plot_saltmp_diff()
savefig(figdir+'merging_overview/CF1_saltmpdiff.png')

ptmp_r50={}
ptmp_r50['from100']=dat18.PTMP.sel(DEPTH=100)+mean(dat18.PTMP.sel(DEPTH=50)-dat18.PTMP.sel(DEPTH=100))
ptmp_r50['from200']=dat18.PTMP.sel(DEPTH=200)+mean(dat18.PTMP.sel(DEPTH=50)-dat18.PTMP.sel(DEPTH=200))


def toTimestamp(d):
  return calendar.timegm(d.timetuple())

def to_datetime_andordinal(date):
    """
    Converts a numpy datetime64 object to a python datetime object
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    dtimeobj=datetime.datetime.utcfromtimestamp(timestamp)
    return datetime.datetime.toordinal(dtimeobj)

time_month=[to_datetime_andordinal(ddd) for ddd in mtime.values]
time_all=[to_datetime_andordinal(ddd) for ddd in dat18.TIME.values]

f100=interp1d(time_month,sal_mdiff_int['from100'],bounds_error=False)
f200=interp1d(time_month,sal_mdiff_int['from200'],bounds_error=False)

sal_mdiff_fulltime={}
sal_mdiff_fulltime['from100']=f100(time_all)
sal_mdiff_fulltime['from200']=f200(time_all)

psal_r50={}
psal_r50['from100']=dat18.PSAL.sel(DEPTH=100)+sal_mdiff_fulltime['from100']
psal_r50['from200']=dat18.PSAL.sel(DEPTH=200)+sal_mdiff_fulltime['from200']



def comp_from100200():
    f,[ax1,ax2]=subplots(2,1,figsize=(12,10))
    ptmp_r50['from100'].plot(ax=ax1,label='from 100m')
    ptmp_r50['from200'].plot(ax=ax1,label='from 200m')
    dat18.PTMP.sel(DEPTH=50).plot(ax=ax1,label='directly measured')
    ax1.legend()
    ax1.set_title('')
    psal_r50['from100'].plot(ax=ax2,label='from 100m')
    psal_r50['from200'].plot(ax=ax2,label='from 200m')
    dat18.PSAL.sel(DEPTH=50).plot(ax=ax2,label='directly measured')
    ax2.set_title('')


comp_from100200()


###############################################################################
##### Save a reconstructed product which keeps recorder 50m data
#### And adds the 100m reconstruction beyond that
#### This is simply because the 100m is present,closer and noisier
#### It wouldn't make sense for 50m instrument to have less variability
################################################################################

dat18_recon=dat18.copy()
dat18_recon['PSAL'].sel(DEPTH=50)[isnan(dat18['PSAL'].sel(DEPTH=50))]=psal_r50['from100'][isnan(dat18['PSAL'].sel(DEPTH=50))].values
dat18_recon['PTMP'].sel(DEPTH=50)[isnan(dat18['PTMP'].sel(DEPTH=50))]=ptmp_r50['from100'][isnan(dat18['PTMP'].sel(DEPTH=50))].values
dat18_recon['PRES'].sel(DEPTH=50)[isnan(dat18['PRES'].sel(DEPTH=50))]=mean(dat18['PRES'].sel(DEPTH=50))

plot(dat18_recon['PRES']);

plot_overview(dat16_recon,dat18_recon)
savefig(figdir+'merging_overview/CF1_overview_1618recon.png')

dat18_recon.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF1_mcat_recon_2018recovery_daily.nc','w',format='netCDF4')
