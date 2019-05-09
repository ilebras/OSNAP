#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1804shelf.pickle','rb'))

daily=newgrid.copy()

# mask the fields based on bathymetry file (this version does not have extra fields on either side, thats just for plotting)
#
# bathf=interpolate.interp1d(bathdist,bathbath)
# bathonmygrid=bathf(newgrid.distance)
#
# daily=newgrid.copy()
# for vv in newgrid:
#     if vv[0]!='d':
#         print(vv)
#         for dd,adist in enumerate(daily.distance):
#             daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=bathonmygrid[dd])

# daily.distance
# daily.depth
# dist1=45
# dist2=95
# depth1=150
# depth2=350
#
# daily['temperature'].sel(distance=slice(dist1,dist2)).isel(depth=slice(depth1,depth2)).mean(dim='distance').mean(dim='depth').plot()
# daily['temperature'].sel(distance=slice(dist1,dist2)).isel(depth=slice(depth1,depth2)).mean(dim='distance').mean(dim='depth').resample('M',dim='date').plot(linewidth=3)


#################################################################################
###################### (Freshwater) Transport #################################
#################################################################################

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))

mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))



# srefa=34.8
srefa=34.9189 # for lazy osnap recalc
srefb=34.9
# srefc=35
srefc=34.9682 # for lazy osnap recalc

##adding 10km width to CF1

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

onesxr=daily.salinity/daily.salinity

distvec
daily.distance[:9]

sep=9

cc={}
cc['trans']=daily.xport[:sep,:,:].sum('depth').sum('distance')
cc['trans plus']=daily['xport plus'][:sep,:,:].sum('depth').sum('distance')


cc['fresha']=(daily['xport'][:sep,:-1,:]*1e3*(daily.salinity[:sep,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
cc['freshb']=(daily['xport'][:sep,:-1,:]*1e3*(daily.salinity[:sep,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
cc['fresha plus']=(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
cc['freshb plus']=(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
cc['freshc plus']=(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-srefc)/srefc).sum('depth').sum('distance')

cc['freshc']=(daily['xport'][:sep,:-1,:]*1e3*(daily.salinity[:sep,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
cc['area']=(onesxr[:sep,:-1,:]*depthdiffmat[:sep,:,:]*middistmat[:sep,:,:]/1e3).sum('depth').sum('distance')
cc['sal']=(daily['xport'][:sep,:-1,:]*daily['salinity'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['tmp']=(daily['xport'][:sep,:-1,:]*daily['temperature'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['den']=(daily['xport'][:sep,:-1,:]*daily['potential density'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['meanvel']=cc['trans']/cc['area']

egic={}
egic['trans']=daily['xport over 27.8'][sep:,:,:].sum('depth').sum('distance')
egic['fresha']=(daily['xport over 27.8'][sep:,:-1,:]*1e3*(daily.salinity[sep:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
egic['freshb']=(daily['xport over 27.8'][sep:,:-1,:]*1e3*(daily.salinity[sep:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
egic['freshc']=(daily['xport over 27.8'][sep:,:-1,:]*1e3*(daily.salinity[sep:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
egic['area']=(onesxr.where(daily['potential density']<27.8)[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport over 27.8'][sep:,:-1,:]*daily['salinity'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport over 27.8'][sep:,:-1,:]*daily['temperature'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport over 27.8'][sep:,:-1,:]*daily['potential density'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']

ic={}
ic['trans']=daily['xport over 27.8'][sep:,:,:].where(daily.salinity>srefb).sum('depth').sum('distance')
ic['fresha']=(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
ic['freshb']=(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
ic['freshc']=(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
# ic['area']=(onesxr[6:,:-1,:].where(daily.salinity>srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
# ic['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['den']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['meanvel']=ic['trans']/ic['area']
47+19

eg={}
eg['trans']=daily['xport over 27.8'][sep:,:,:].where(daily.salinity<=srefb).sum('depth').sum('distance')
eg['fresha']=(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
eg['freshb']=(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
eg['freshc']=(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
# eg['area']=(onesxr[6:,:-1,:].where(daily.salinity<=srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
# eg['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['den']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['meanvel']=eg['trans']/eg['area']

egic['freshc'].sel(date=slice('2014-08-01','2016-04-30'))
(cc['freshc']+egic['freshc']).sel(date=slice('2014-08-01','2016-04-30')).mean()
(cc['freshc plus']+egic['freshc']).sel(date=slice('2014-08-01','2016-04-30')).mean()
80/141
107/141

pickle.dump([cc,egic,eg,ic],open('../pickles/transdic.pickle','wb'))

def decompose(cur,srefchoose):
    cur['abar']=nanmean(cur['area'])
    cur['vbar']=nanmean(cur['meanvel'])
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']
    cur['sbar']=nanmean((cur['sal']-srefchoose)/srefchoose)
    cur['sprime']=((cur['sal']-srefchoose)/srefchoose-cur['sbar'])
    cur['vpAS']=(cur['sbar']*cur['abar']*cur['vprime']+cur['sbar']*cur['abar']*cur['vbar'])*1e3
    cur['VAsp']=(cur['sprime']*cur['abar']*cur['vbar']+cur['sbar']*cur['abar']*cur['vbar'])*1e3
    return cur

cc=decompose(cc,srefb)
egic=decompose(egic,srefb)


def colorstripes():
    axvspan(datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,9,1),datetime.datetime(2015,12,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2015,12,1),datetime.datetime(2016,3,1),color=egcol,alpha=0.4)


#################################################################################
######## Make 3 panel (freshwater) transport decomp for each current
####################################################################################

def eachpanel(field,colo,axit,labit='',xr=daily,pnofilt=1,letlab=''):
    if pnofilt==1:
        field.plot(alpha=0.5,color=colo,label='',ax=axit)
    axit.plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit)
    axit.set_xlabel('')
    axit.set_ylabel('')
    axit.text(0.01, 0.85,letlab,transform = axit.transAxes,fontsize=15)

def plot3pan_coastal():
    f, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True, figsize=(8,7))
    ##2
    eachpanel(cc['trans'],ccol,ax1)
    eachpanel(cc['trans plus'],'grey',ax1,pnofilt=0,letlab='a)')
    ax1.set_ylim([-1.5,0.1])

    ##2
    eachpanel(cc['freshb'],ccol,ax2,labit='Between moorings')
    eachpanel(cc['freshb plus'],'grey',ax2,pnofilt=0,labit='Including shelf',letlab='b)')
    ax2.legend(loc=(0.27,0.95),framealpha=1)
    ax2.set_ylim([-5,80])
    ##3
    eachpanel(cc['freshb'],ccol,ax3,labit='vAS',pnofilt=0,letlab='c)')
    eachpanel(cc['vpAS'],'k',ax3,pnofilt=0,labit='$v\' \overline{AS}$')
    eachpanel(cc['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{vA}S\'$')
    ax3.legend()
    ax3.set_ylim(12,32)
    fs=14
    xlab=-0.01
    f.text(xlab, 0.85, 'Transport [Sv]', va='center', rotation='vertical',fontsize=fs)
    f.text(xlab, 0.4, 'Freshwater transport [mSv]', va='center', rotation='vertical',fontsize=fs)
    plt.tight_layout()
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
    ax1.set_title('Coastal current',fontsize=fs+4)
    savefig('../figures/paperfigs/coastal_tseries.pdf',bbox_inches='tight')

plot3pan_coastal()

def plot3pan_slope():
    f, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True, figsize=(8,7))
    ##2
    eachpanel(egic['trans'],egicol,ax1)
    eachpanel(eg['trans'],'k',ax1,pnofilt=0)
    eachpanel(ic['trans'],'grey',ax1,pnofilt=0,letlab='a)')
    ax1.set_ylim([-30,0])
    ##2
    eachpanel(egic['freshb'],egicol,ax2)
    eachpanel(eg['freshb'],'k',ax2,pnofilt=0,labit='East Greenland Current (S < 34.9)')
    eachpanel(ic['freshb'],'grey',ax2,pnofilt=0,labit='Irminger Current (S > 34.9)',letlab='b)')
    ax2.axhline(0,color='k')
    ax2.legend(loc=(0.3,0.95),framealpha=1)
    ax2.set_ylim([-35,120])
    ##3
    eachpanel(egic['freshb'],egicol,ax3,labit='vAS',pnofilt=0)
    eachpanel(egic['vpAS'],'k',ax3,pnofilt=0,labit='$v\' \overline{AS}$')
    eachpanel(egic['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{vA}S\'$',letlab='c)')
    ax3.legend()
    ax3.set_ylim(-15,70)
    fs=14
    xlab=-0.01
    f.text(xlab, 0.85, 'Transport [Sv]', va='center', rotation='vertical',fontsize=fs)
    f.text(xlab, 0.4, 'Freshwater transport [mSv]', va='center', rotation='vertical',fontsize=fs)
    plt.tight_layout()
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
    ax1.set_title('Slope current',fontsize=fs+4)
    savefig('../figures/paperfigs/slope_tseries.pdf',bbox_inches='tight')

plot3pan_slope()

def compfresh():
    fig=figure(figsize=(10,3.5))
    ax=subplot(111)
    freshtype='freshb'
    fill_between(daily.date.values, sig.filtfilt(B,A, cc[freshtype]), sig.filtfilt(B,A,cc[freshtype+' plus']),color=ccol,label='',alpha=0.6)
    eachpanel(cc[freshtype],ccol,ax,pnofilt=0,labit='Coastal current')
    eachpanel(eg[freshtype],egcol,ax,pnofilt=0,labit='Slope current < 34.9')
    eachpanel(ic[freshtype],icol,ax,pnofilt=0,labit='Slope current > 34.9')
    axhline(0,color='k')
    # axvline(datetime.datetime(2015,7,1,0),color='grey',linewidth=0.8)
    # pwf(cc[freshtype]+egic[freshtype],'grey',1)
    ylim([-30,95])
    ylabel('Freshwater transport [mSv]')
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
    colorstripes()
    fig.autofmt_xdate()
    dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]
    gca().set_xticks(dticks)
    # legend(loc=(1.05,0.5),framealpha=1)
    legend()
    labyy=80
    text(datetime.datetime(2014,10,15),labyy,'FALL',color=ccol,fontsize=13.5)
    text(datetime.datetime(2015,1,10),labyy,'WINTER',color=egcol,fontsize=13.5)
    text(datetime.datetime(2015,4,3),labyy,'SPRING/SUMMER',color='grey',fontsize=13.5)
    # text(datetime.datetime(2015,7,5),75,'SUMMER',color='grey',fontsize=15)
    savefig('../figures/paperfigs/freshcomp_symm.pdf',bbox_inches='tight')

compfresh()

def refcomp():
    f, ((ax1), (ax2), (ax3), (ax4)) = plt.subplots(4, 1, sharex=True, figsize=(8,9))
    eachpanel(cc['fresha'],'b',ax1,pnofilt=0,labit='34.8')
    eachpanel(cc['freshb'],'r',ax1,pnofilt=0,labit='34.9')
    eachpanel(cc['freshc'],'g',ax1,pnofilt=0,labit='35')
    ax1.set_title('Coastal Current')

    eachpanel(eg['fresha'],'b',ax2,pnofilt=0,labit='34.8')
    eachpanel(eg['freshb'],'r',ax2,pnofilt=0,labit='34.9')
    eachpanel(eg['freshc'],'g',ax2,pnofilt=0,labit='35')
    ax2.legend()
    ax2.set_title('Slope Current < 34.9')

    eachpanel(ic['fresha'],'b',ax3,pnofilt=0,labit='34.8')
    eachpanel(ic['freshb'],'r',ax3,pnofilt=0,labit='34.9')
    eachpanel(ic['freshc'],'g',ax3,pnofilt=0,labit='35')
    ax3.set_title('Slope Current > 34.9')
    ax3.set_ylabel('Freshwater transport [mSv]')

    eachpanel(egic['fresha']+cc['fresha'],'b',ax4,pnofilt=0,labit='34.8')
    eachpanel(egic['freshb']+cc['freshb'],'r',ax4,pnofilt=0,labit='34.9')
    eachpanel(egic['freshc']+cc['freshc'],'g',ax4,pnofilt=0,labit='35')
    ax4.set_title('Coastal + Slope')


refcomp()
savefig('../notes/section_details/transfigs/refsalcomp.pdf',bbox_inches='tight')



#################################################################################
#### TS diagrams
#################################################################################
def TSplot():
        hexbin(daily['salinity'].where(daily.distance>0).values.flatten(),daily['temperature'].where(daily.distance>0).values.flatten(),cmap=cm.hot_r,bins='log',gridsize=200)
        colorbar(label='[log # measurements]')
        xlim([32,35.2])
        ylim([-2,10])
        xlabel('Salinity')
        plot([32.5,32.5,34.25,34.25,32.5],[-1.8,0,0,-1.8,-1.8],color='k')
        den=contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(24,28,0.5))
        manual_locations=[(32,2),(33,2),(33,2),(34,2)]
        clabel(den,fmt='%1.1f')#,manual=manual_locations)
        text(34.95,7,'AW',color='k',fontsize=20)
        text(34.0,0.11,'PW',color='k',fontsize=20)
        ylabel('Potential temperature [$^\circ$C]')
        axvline(34.9,color='k')
        savefig('../figures/paperfigs/TSall.pdf',bbox_inches='tight')

TSplot()


def partTS(date1,date2,textit,ax2,velvers=0):
    # if velvers==1:
    #     ax2.hexbin(daily['salinity'].values.flatten(),daily['temperature'].values.flatten(),cmap=cm.Greys,bins='log',gridsize=200)
    #     d2=500
    #     axpart=ax2.hexbin(daily['salinity'].sel(date=slice(date1,date2),depth=slice(0,d2)).values.flatten(),
    #     daily['temperature'].sel(date=slice(date1,date2),depth=slice(0,d2)).values.flatten(),
    #     C=daily['across track velocity'].sel(date=slice(date1,date2),depth=slice(0,d2)).values.flatten(),
    #     vmin=-0.6,vmax=0,cmap=univec['uacross'][2],gridsize=200)
    # else:
    ax2.hexbin(daily['salinity'].where(daily.distance>0).values.flatten(),daily['temperature'].where(daily.distance>0).values.flatten(),cmap=cm.Greys,bins='log',gridsize=200,rasterized='True')
    axpart=ax2.hexbin(daily['salinity'].sel(date=slice(date1,date2)).values.flatten(),daily['temperature'].sel(date=slice(date1,date2)).values.flatten(),
    cmap=cm.hot_r,bins='log',mincnt=1,gridsize=200,rasterized='True')
    ax2.set_xlim([32,35.2])
    ax2.set_ylim([-2,8])
    ax2.plot([32.5,32.5,34.25,34.25,32.5],[-1.8,0,0,-1.8,-1.8],color='k')
    # ax2.contour(salvec,tmpvec,pdenmat,colors='grey',zorder=50,levels=arange(24,28,0.5))
    # ax2.text(34.9495,6,'AW',color='k',fontsize=16,zorder=60)
    # ax2.text(34.3,-1,'PW',color='k',fontsize=16,zorder=60)
    ax2.axvline(34.9,color='k')
    ax2.text(32.25,6,textit,fontsize=18)

    return axpart

def plot4TS_seas(vv=0):
    f, ((ax11, ax22), (ax33, ax44)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(9.5,7))
    partTS('2014-10-01','2015-1-1','FALL',ax11,vv)
    partTS('2015-01-01','2015-04-1','WINTER',ax22,vv)
    partTS('2015-04-01','2015-7-1','SPRING',ax33,vv)
    axts=partTS('2015-7-01','2015-10-1','SUMMER',ax44,vv)
    f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=16)
    plt.tight_layout()
    cbaxes = f.add_axes([1.01, 0.15, 0.025, 0.7])
    # if vv==1:
    #     cbar=colorbar(axts,label='Velocity in top 500m [m/s]',cax=cbaxes)
    #     savefig('../figures/paperfigs/TS_seasonal_vel.pdf',bbox_inches='tight')
    if vv==0:
        cbar=colorbar(axts,label='[log # measurements]',cax=cbaxes)
        savefig('../figures/paperfigs/TS_seasonal.pdf',bbox_inches='tight')

# plot4TS_seas(1)

plot4TS_seas(0)

def plot4TS_seas_symm(vv=0):
    f, (ax11, ax22, ax33) = plt.subplots(3, 1, sharex=True, sharey=True,figsize=(5,10))
    partTS('2014-10-01','2015-1-1','FALL',ax11,vv)
    partTS('2015-01-01','2015-04-1','WINTER',ax22,vv)
    ax22.text(34.3,-1,'PW',color='k',fontsize=16)
    axts=partTS('2015-04-01','2015-9-1','SPRING/SUMMER',ax33,vv)
    f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=16)
    plt.tight_layout()
    cbaxes = f.add_axes([1.01, 0.3, 0.04, 0.4])
    cbar=colorbar(axts,label='[log # measurements]',cax=cbaxes)
    savefig('../figures/paperfigs/TS_seasonal_symm.pdf',bbox_inches='tight',)

plot4TS_seas_symm()

#############################################################
############## Current vel profiles  ######################
#############################################################

def pltprf(dd,date1,date2,labit,field):
    if dd==0:
        plot(daily[field].sel(date=slice(date1,date2)).sel(distance=0).mean(dim='date'),daily.depth,label=labit,linewidth=3)
    else:
        skipna=10
        plot(daily[field][:,::skipna,:].sel(date=slice(date1,date2)).sel(distance=45).mean(dim='date'),daily.depth[::skipna],label=labit,linewidth=3)


def pltseas(dd,field):
    pltprf(dd,'2014-10-01','2015-1-1','Fall',field)
    pltprf(dd,'2015-1-01','2015-4-1','Winter',field)
    pltprf(dd,'2015-4-01','2015-9-1','Spring/Summer',field)

def velprfcmp():
    fs=14
    f=figure(figsize=(6,7))
    subplot(221)
    pltseas(0,'across track velocity')
    ylim([180,0])
    f.text(-0.03, 0.82, 'Coastal current', rotation='vertical',fontsize=fs+1)
    text(-0.88, 15, 'a)',fontsize=fs+1)
    xlim(-0.9,0)
    subplot(223)
    pltseas(1,'across track velocity')
    ylim([1100,0])
    f.text(0, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fs)
    text(-0.88, 100, 'c)',fontsize=fs+1)
    xlim(-0.9,0)
    xlabel('velocity [m/s]',fontsize=fs)
    f.text(-0.03, 0.35, 'Slope current', rotation='vertical',fontsize=fs+1)
    subplot(222)
    pltseas(0,'salinity')
    text(32.05, 15, 'b)',fontsize=fs+1)
    xlim(32,34.2)
    ylim([180,0])
    gca().set_yticklabels('')
    subplot(224)
    pltseas(1,'salinity')
    text(34.37, 100, 'd)',fontsize=fs+1)
    xlim(34.35,35)
    gca().set_yticklabels('')
    legend(loc=(1.05,0.9),fontsize=fs-1)
    ylim([1100,0])
    xlabel('salinity',fontsize=fs)
    savefig('../figures/paperfigs/Seas_velandsalprofs.pdf',bbox_inches='tight')


velprfcmp()

#### Thinking about error and reading off seasonal max/min below

def autocorr(x,t=1):
    y = x-mean(x)
    norm=sum(y**2).values
    print(norm)
    result=correlate(y,y,mode='full')/norm
    figure(figsize=(10,3))
    ylabel('Autocorrelation')
    plot(x.date,result[len(x)-1:])
    xlim([datetime.datetime(2014,8,17),datetime.datetime(2015,4,17)])
    axhline(0,color='k')
    figure(figsize=(10,3))
    ylabel('Integral time scale [days]')
    plot(x.date,cumsum(result[len(x)-1:]))
    # xlim([datetime.datetime(2014,8,17),datetime.datetime(2015,4,17)])
    # ylim([-10,10])
    axhline(0,color='k')

#Can use above function to look at all the
#I'm going to choose a 5day integral time scale to be semi-conservative.
#Also note that seasonality dominates some of these time series, but I'm going to say that's too long.
autocorr(egic['trans'])

len(daily.date)
#effective degrees of freedom
Neff=70
def geterr(field):
    serr=sqrt(sum([(ss-mean(field))**2 for ss in field])/Neff)
    return serr

geterr(egic['trans'])


XXXXXXXXX
# Report max and min vals:

## For EGC portion transport

egic['trans'].mean()
egic['trans'].std()
eg['trans'].mean()
eg['trans'].std()
ic['trans'].mean()
ic['trans'].std()
print('first max')
## First max
mdate1=datetime.datetime(2014,8,1)
mdate2=datetime.datetime(2015,11,1)
print(max(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)
print('second max')
## Second max
mdate1=datetime.datetime(2015,6,1)
mdate2=datetime.datetime(2016,8,1)
print(max(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)
print('min')
## Min
mdate1=datetime.datetime(2015,2,1)
mdate2=datetime.datetime(2016,1,1)
print(min(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)


srefb

### For all portions freshwater trans

for field in [egic,eg]:
    print('mean and std')
    print(field['freshb'].mean().values)
    print(field['freshb'].std().values)
    print('first max')
    ## First max
    mdate1=datetime.datetime(2014,8,1)
    mdate2=datetime.datetime(2015,11,1)
    print(max(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2))))))
    print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
    print('second max')
    ## Second max
    mdate1=datetime.datetime(2015,6,1)
    mdate2=datetime.datetime(2016,8,1)
    print(max(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2))))))
    print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
    print('min')
    ## Min
    mdate1=datetime.datetime(2014,12,1)
    mdate2=datetime.datetime(2016,1,1)
    print(min(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2)))))
    print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2))))].values)




print('mean and std')
print(ic['freshb'].mean().values)
print(ic['freshb'].std().values)
print('first max')
## First max
mdate1=datetime.datetime(2014,8,1)
mdate2=datetime.datetime(2015,11,1)
print(min(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2))))))
print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
print('second max')
## Second max
mdate1=datetime.datetime(2015,6,1)
mdate2=datetime.datetime(2016,8,1)
print(min(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2))))))
print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
print('min')
## Min
mdate1=datetime.datetime(2014,12,1)
mdate2=datetime.datetime(2016,1,1)
print(max(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2))))))
print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
################################################
#### Coastal section
################################################
for key in cc:
    if ('trans' in key) | ('fresh' in key):
            print(key)
            print('mean=',cc[key].mean().values)
            print('std=',cc[key].std().values)

date1=datetime.datetime(2014,8,15)
## Print max and min (freshwater) transports
middate=datetime.datetime(2015,5,1,0)
for key in cc:
    if ('trans' in key) | ('fresh' in key):
            print(key)
            print('mean=',cc[key].mean().values)
            print('std=',cc[key].std().values)
            print(max(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(date1,middate))))))
            print(daily.date[argmax(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(date1,middate)))))].values)


lastdate=datetime.datetime(2016,7,30)
for key in cc:
    if ('trans' in key) | ('fresh' in key):
        print(key)
        print(max(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(middate,lastdate))))))
        print(daily.date.sel(date=slice(middate,lastdate))[argmax(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(middate,lastdate)))))].values)

mdate1=datetime.datetime(2014,10,1)
mdate2=datetime.datetime(2015,9,1)
for key in cc:
    if ('trans' in key) | ('fresh' in key):
        print(key)
        print(min(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(mdate1,mdate2))))))
        print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(mdate1,mdate2)))))].values)



## Print (freshwater) transport for first two weeks of August

daily.date

date2=datetime.datetime(2014,9,1)
for key in cc:
    if ('trans' in key) | ('fresh' in key):
        print(key)
        print(cc[key].sel(date=slice(date1,date2)).mean().values)
