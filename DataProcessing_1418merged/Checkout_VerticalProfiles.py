from firstfuncs_1618 import *


dat={}
for ii in range(1,9):
    dat['cf'+str(ii)]=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(ii)+'_mcat_vertgrid_daily.nc')


def plotprofs(dic,moorname,date1,date2):
    f,axx=subplots(1,3,figsize=(12,4),sharey=True)
    d1=0
    d2=2000
    axx[0].invert_yaxis()
    axx[0].set_ylabel('depth [m]')
    for ii,prop in enumerate(['PSAL','PTMP','PDEN']):
        axx[ii].plot(dic[prop].sel(depth=slice(d1,d2)).sel(date=slice(date1,date2)).T,dic.depth.sel(depth=slice(d1,d2)))
        axx[ii].plot(dic[prop].sel(depth=slice(d1,d2)).sel(date=slice(date1,date2)).mean(dim='date'),dic.depth.sel(depth=slice(d1,d2)),color='k',linewidth=3)
    axx[0].set_xlabel('salinity')
    axx[0].set_xlim(34.5,35.1)
    axx[1].set_xlabel('pot. temp [$^{\circ}$C]')
    axx[2].set_xlabel('pot. density [kg m$^{-3}$]')
    axx[1].set_title(date1+' to '+date2)
    savefig(figdir+'profiles_sal_investigation/'+moorname+'_Profiles_'+date1+'to'+date2+'.png',bbox_inches='tight')

# Loop through and get seasonal:
mdates=[str(ii)[:10] for ii in dat['cf1'].resample(date='3M').mean(dim='date').date.values]

for ii in range(len(mdates)-1):
    for jj in range(5,9):
        plotprofs(dat['cf'+str(jj)],'CF'+str(jj),mdates[ii],mdates[ii+1])

len(mdates)


colvec=4*['C0','C1','C2','C3']
lsty=array([4*['-'],4*['--'],4*['-.'],4*[':']]).flatten()
lsty

# For each mooring, plot seasonal mean on same figure:
def plotmeanonly(axx,dic,date1,date2,jj):
    d1=0
    d2=2000
    axx[0].set_ylabel('depth [m]')
    for ii,prop in enumerate(['PSAL','PTMP','PDEN']):
        axx[ii].plot(dic[prop].sel(depth=slice(d1,d2)).sel(date=slice(date1,date2)).mean(dim='date'),dic.depth.sel(depth=slice(d1,d2)),color=colvec[jj],linewidth=2,label=date1,linestyle=lsty[jj])
    axx[0].set_xlabel('salinity')
    axx[0].set_xlim(34.8,35.05)
    axx[1].set_xlabel('pot. temp [$^{\circ}$C]')
    axx[2].set_xlabel('pot. density [kg m$^{-3}$]')


for kk in range(5,9):
    f,axx=subplots(1,3,figsize=(12,4),sharey=True)
    for jj in range(16):
        plotmeanonly(axx,dat['cf'+str(kk)],mdates[jj],mdates[jj+1],jj)
    legend(loc=(1.05,-0.05))
    axx[0].invert_yaxis()
    axx[1].set_title('CF'+str(kk))
    savefig(figdir+'profiles_sal_investigation/MeanSeasProfiles_CF'+str(kk)+'.png',bbox_inches='tight')


ydates=['2014-8-31','2015-8-31','2016-8-31','2017-8-31','2018-8-31']

for kk in range(5,9):
    f,axx=subplots(1,3,figsize=(12,4),sharey=True)
    for jj in range(4):
        plotmeanonly(axx,dat['cf'+str(kk)],ydates[jj],ydates[jj+1],jj)
    legend(loc=(1.05,0.2))
    axx[0].invert_yaxis()
    axx[1].set_title('CF'+str(kk))
    savefig(figdir+'profiles_sal_investigation/MeanYrlyProfiles_CF'+str(kk)+'.png',bbox_inches='tight')



############### CF67_deepsal:

cf6=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF6_mcat_vertgrid_daily.nc')
cf7=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CF7_mcat_vertgrid_daily.nc')


def plotprofs(dic,moorname,date1,date2):
    f,axx=subplots(1,3,figsize=(12,4),sharey=True)
    d1=300
    d2=2000
    axx[0].invert_yaxis()
    axx[0].set_ylabel('depth [m]')
    for ii,prop in enumerate(['PSAL','PTMP','PDEN']):
        axx[ii].plot(dic[prop].sel(depth=slice(d1,d2)).sel(date=slice(date1,date2)).T,dic.depth.sel(depth=slice(d1,d2)))
        axx[ii].plot(dic[prop].sel(depth=slice(d1,d2)).sel(date=slice(date1,date2)).mean(dim='date'),dic.depth.sel(depth=slice(d1,d2)),color='k',linewidth=3)
    axx[0].set_xlabel('salinity')
    axx[0].set_xlim(34.88,34.98)
    axx[1].set_xlabel('pot. temp [$^{\circ}$C]')
    axx[2].set_xlabel('pot. density [kg m$^{-3}$]')
    axx[1].set_title(date1+' to '+date2)
    savefig(figdir+'merging_overview_mcats/CF67_deepsal/'+moorname+'_Profiles_'+date1+'to'+date2+'.png',bbox_inches='tight')


plotprofs(cf6,'CF6','2014-09-01','2018-08-01')
plotprofs(cf6,'CF6','2016-09-01','2018-08-01')
plotprofs(cf7,'CF7','2014-09-01','2018-08-01')
plotprofs(cf7,'CF7','2014-09-01','2016-08-01')
plotprofs(cf7,'CF7','2016-09-01','2018-08-01')

# Loop through and get monthly plots:
mdates=[str(ii)[:10] for ii in cf7.resample(date='M').mean(dim='date').date.values]

for ii in range(len(mdates)-1):
    plotprofs(cf6,'CF6',mdates[ii],mdates[ii+1])
    # plotprofs(cf7,'CF7',mdates[ii],mdates[ii+1])
