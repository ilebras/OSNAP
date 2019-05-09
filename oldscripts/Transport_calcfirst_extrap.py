#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1803extrap.pickle','rb'))

# mask the fields based on bathymetry file (this version does not have extra fields on either side, thats just for plotting)

bathf=interpolate.interp1d(bathdist,bathbath)
bathonmygrid=bathf(newgrid.distance)

daily=newgrid.copy()
for vv in newgrid:
    if vv[0]!='d':
        print(vv)
        for dd,adist in enumerate(daily.distance):
            daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=bathonmygrid[dd])




#################################################################################
###################### (Freshwater) Transport #################################
#################################################################################

mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
mid_dist[1]=1.5 ## being strict -- only on one side of mooring
mid_dist[-1]=2.5 ## being strict -- only on one side of mooring
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))


srefa=34.8
srefb=34.9
srefc=35

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3

##adding 10km width to CF1
mid_dist_plus=mid_dist.copy()
mid_dist_plus[0]=10
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

onesxr=daily.salinity/daily.salinity

cc={}
cc['trans']=daily.xport[:6,:,:].sum('depth').sum('distance')
cc['trans plus']=daily['xport plus'][:6,:,:].sum('depth').sum('distance')

cc['fresha']=(daily['xport'][:6,:-1,:]*1e3*(daily.salinity[:6,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
cc['freshb']=(daily['xport'][:6,:-1,:]*1e3*(daily.salinity[:6,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
cc['freshb plus']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat_plus[:6,:]*(daily.salinity[:6,:-1,:]-srefb)/srefb).sum('depth').sum('distance')

cc['freshc']=(daily['xport'][:6,:-1,:]*1e3*(daily.salinity[:6,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
cc['area']=(onesxr[:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('depth').sum('distance')
cc['sal']=(daily['xport'][:6,:-1,:]*daily['salinity'][:6,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['tmp']=(daily['xport'][:6,:-1,:]*daily['temperature'][:6,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['den']=(daily['xport'][:6,:-1,:]*daily['potential density'][:6,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['meanvel']=cc['trans']/cc['area']

egic={}
egic['trans']=daily['xport over 27.8'][6:,:,:].sum('depth').sum('distance')
egic['fresha']=(daily['xport over 27.8'][6:,:-1,:]*1e3*(daily.salinity[6:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
egic['freshb']=(daily['xport over 27.8'][6:,:-1,:]*1e3*(daily.salinity[6:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
egic['freshc']=(daily['xport over 27.8'][6:,:-1,:]*1e3*(daily.salinity[6:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
egic['area']=(onesxr.where(daily['potential density']<27.8)[6:,:-1,:]*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport over 27.8'][6:,:-1,:]*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport over 27.8'][6:,:-1,:]*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport over 27.8'][6:,:-1,:]*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']

ic={}
ic['trans']=daily['xport over 27.8'][6:,:,:].where(daily.salinity>srefb).sum('depth').sum('distance')
ic['fresha']=(daily['xport over 27.8'][6:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[6:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
ic['freshb']=(daily['xport over 27.8'][6:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[6:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
ic['freshc']=(daily['xport over 27.8'][6:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[6:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
# ic['area']=(onesxr[6:,:-1,:].where(daily.salinity>srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
# ic['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['den']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['meanvel']=ic['trans']/ic['area']


eg={}
eg['trans']=daily['xport over 27.8'][6:,:,:].where(daily.salinity<=srefb).sum('depth').sum('distance')
eg['fresha']=(daily['xport over 27.8'][6:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[6:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
eg['freshb']=(daily['xport over 27.8'][6:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[6:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
eg['freshc']=(daily['xport over 27.8'][6:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[6:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
# eg['area']=(onesxr[6:,:-1,:].where(daily.salinity<=srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
# eg['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['den']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['meanvel']=eg['trans']/eg['area']

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
    ax2.set_ylim([-10,70])
    ##3
    eachpanel(cc['freshb'],ccol,ax3,labit='vAS',pnofilt=0,letlab='c)')
    eachpanel(cc['vpAS'],'k',ax3,pnofilt=0,labit='$v\' \overline{AS}$')
    eachpanel(cc['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{vA}S\'$')
    ax3.legend()
    ax3.set_ylim(13,28)
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
    ax2.set_ylim([-35,100])
    ##3
    eachpanel(egic['freshb'],egicol,ax3,labit='vAS',pnofilt=0)
    eachpanel(egic['vpAS'],'k',ax3,pnofilt=0,labit='$v\' \overline{AS}$')
    eachpanel(egic['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{vA}S\'$',letlab='c)')
    ax3.legend()
    ax3.set_ylim(-20,60)
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
    ylim([-35,95])
    ylabel('Freshwater transport [mSv]')
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
    colorstripes()
    fig.autofmt_xdate()
    dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]
    gca().set_xticks(dticks)
    # legend(loc=(1.05,0.5),framealpha=1)
    legend()
    text(datetime.datetime(2014,10,1),75,'FALL',color=ccol,fontsize=13.5)
    text(datetime.datetime(2015,1,2),75,'WINTER',color=egcol,fontsize=13.5)
    text(datetime.datetime(2015,4,1),75,'SPRING/SUMMER',color='grey',fontsize=13.5)
    # text(datetime.datetime(2015,7,5),75,'SUMMER',color='grey',fontsize=15)
    savefig('../figures/paperfigs/freshcomp_symm.pdf',bbox_inches='tight')

compfresh()

srefc
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
        hexbin(daily['salinity'].values.flatten(),daily['temperature'].values.flatten(),cmap=cm.hot_r,bins='log',gridsize=200)
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



help(hexbin)

def partTS(date1,date2,textit,ax2,velvers=0):
    # if velvers==1:
    #     ax2.hexbin(daily['salinity'].values.flatten(),daily['temperature'].values.flatten(),cmap=cm.Greys,bins='log',gridsize=200)
    #     d2=500
    #     axpart=ax2.hexbin(daily['salinity'].sel(date=slice(date1,date2),depth=slice(0,d2)).values.flatten(),
    #     daily['temperature'].sel(date=slice(date1,date2),depth=slice(0,d2)).values.flatten(),
    #     C=daily['across track velocity'].sel(date=slice(date1,date2),depth=slice(0,d2)).values.flatten(),
    #     vmin=-0.6,vmax=0,cmap=univec['uacross'][2],gridsize=200)
    # else:
    ax2.hexbin(daily['salinity'].values.flatten(),daily['temperature'].values.flatten(),cmap=cm.Greys,bins='log',gridsize=200,rasterized='True')
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

def plot4TS_seas_symm():
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

def pltvelprf(dd,date1,date2,labit):
    if dd==0:
        plot(daily['across track velocity'][:6,:,:].sel(date=slice(date1,date2)).mean(dim='distance').mean(dim='date'),daily.depth,label=labit)
    else:
        skipna=10
        plot(daily['across track velocity'][:,::skipna,:].sel(date=slice(date1,date2)).sel(distance=45).mean(dim='date'),daily.depth[::skipna],label=labit)

def velprfcmp():
    fs=14
    f=figure(figsize=(4,8))
    subplot(211)
    title('Coastal current',fontsize=fs)
    pltvelprf(0,'2014-10-01','2015-1-1','Fall')
    pltvelprf(0,'2015-1-01','2015-4-1','Winter')
    pltvelprf(0,'2015-4-01','2015-9-1','Spring/Summer')
    legend()
    ylim([180,0])
    xlim([-0.7,0])
    gca().set_xticklabels('')
    subplot(212)
    title('Slope current',fontsize=fs)
    pltvelprf(1,'2014-10-01','2015-1-1','Fall')
    pltvelprf(1,'2015-1-01','2015-4-1','Winter')
    pltvelprf(1,'2015-4-01','2015-9-1','Spring/Summer')
    ylim([1100,0])
    xlim([-0.7,0])
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fs)
    xlabel('velocity [m/s]',fontsize=fs)
    savefig('../figures/paperfigs/Velprofs.pdf',bbox_inches='tight')


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
#effective degrees of freedom
Neff=140
def geterr(field):
    serr=sqrt(sum([(ss-mean(field))**2 for ss in field])/Neff)
    return serr

geterr(egic['trans'])


XXXXXXXXX
## Report max and min vals:

# ## For EGC portion transport
#
# print('first max')
# ## First max
# mdate1=datetime.datetime(2014,8,1)
# mdate2=datetime.datetime(2015,11,1)
# print(max(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)
# print('second max')
# ## Second max
# mdate1=datetime.datetime(2015,6,1)
# mdate2=datetime.datetime(2016,8,1)
# print(max(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)
# print('min')
# ## Min
# mdate1=datetime.datetime(2015,2,1)
# mdate2=datetime.datetime(2016,1,1)
# print(min(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)
#
# ### For all portions freshwater trans
#
# for field in [egic,eg]:
#     print('mean and std')
#     print(field['freshb'].mean().values)
#     print(field['freshb'].std().values)
#     print('first max')
#     ## First max
#     mdate1=datetime.datetime(2014,8,1)
#     mdate2=datetime.datetime(2015,11,1)
#     print(max(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2))))))
#     print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
#     print('second max')
#     ## Second max
#     mdate1=datetime.datetime(2015,6,1)
#     mdate2=datetime.datetime(2016,8,1)
#     print(max(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2))))))
#     print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
#     print('min')
#     ## Min
#     mdate1=datetime.datetime(2014,12,1)
#     mdate2=datetime.datetime(2016,1,1)
#     print(min(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2)))))
#     print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(sig.filtfilt(B,A,field['freshb'].sel(date=slice(mdate1,mdate2))))].values)
#
#
#
#
# print('mean and std')
# print(ic['freshb'].mean().values)
# print(ic['freshb'].std().values)
# print('first max')
# ## First max
# mdate1=datetime.datetime(2014,8,1)
# mdate2=datetime.datetime(2015,11,1)
# print(min(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
# print('second max')
# ## Second max
# mdate1=datetime.datetime(2015,6,1)
# mdate2=datetime.datetime(2016,8,1)
# print(min(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
# print('min')
# ## Min
# mdate1=datetime.datetime(2014,12,1)
# mdate2=datetime.datetime(2016,1,1)
# print(max(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,ic['freshb'].sel(date=slice(mdate1,mdate2)))))].values)
