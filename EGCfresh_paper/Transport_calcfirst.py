#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

versname='1810JHcal'
daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_'+versname+'.pickle','rb'))





daily['across track velocity']=-1*daily['across track velocity']
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
srefa=34.9187 # for lazy osnap recalc FULL
srefb=34.9
# srefc=35
srefc=34.9682 # for lazy osnap recalc EAST

##adding 10km width to CF1

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

onesxr=daily.salinity/daily.salinity

sep=9
cc={}
# cc['trans']=daily.xport[:sep,:,:].sum('depth').sum('distance')
cc['trans']=daily['xport plus'][:sep,:,:].sum('depth').sum('distance')
# cc['fresha']=(daily['xport'][:sep,:-1,:]*1e3*(daily.salinity[:sep,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
# cc['freshb']=(daily['xport'][:sep,:-1,:]*1e3*(daily.salinity[:sep,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
# cc['freshc']=(daily['xport'][:sep,:-1,:]*1e3*(daily.salinity[:sep,:-1,:]-srefc)/srefc).sum('depth').sum('distance')


cc['fresha']=-(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
cc['freshb']=-(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
cc['freshc']=-(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-srefc)/srefc).sum('depth').sum('distance')

cc['area']=(onesxr[:sep,:-1,:]*depthdiffmat[:sep,:,:]*middistmat[:sep,:,:]/1e3).sum('depth').sum('distance')
cc['sal']=(daily['xport plus'][:sep,:-1,:]*daily['salinity'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['sal'][abs(cc['sal']-35)>3]=nanmean(cc['sal'])
cc['tmp']=(daily['xport plus'][:sep,:-1,:]*daily['temperature'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['den']=(daily['xport plus'][:sep,:-1,:]*daily['potential density'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['meanvel']=cc['trans']/cc['area']

egic={}
egic['trans']=daily['xport over 27.8'][sep:,:,:].sum('depth').sum('distance')
egic['fresha']=-(daily['xport over 27.8'][sep:,:-1,:]*1e3*(daily.salinity[sep:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
egic['freshb']=-(daily['xport over 27.8'][sep:,:-1,:]*1e3*(daily.salinity[sep:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
egic['freshc']=-(daily['xport over 27.8'][sep:,:-1,:]*1e3*(daily.salinity[sep:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
egic['area']=(onesxr.where(daily['potential density']<27.8)[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport over 27.8'][sep:,:-1,:]*daily['salinity'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport over 27.8'][sep:,:-1,:]*daily['temperature'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport over 27.8'][sep:,:-1,:]*daily['potential density'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']

ic={}
ic['trans']=daily['xport over 27.8'][sep:,:,:].where(daily.salinity>srefb).sum('depth').sum('distance')
ic['fresha']=-(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
ic['freshb']=-(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
ic['freshc']=-(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
# ic['area']=(onesxr[6:,:-1,:].where(daily.salinity>srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
# ic['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['den']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
# ic['meanvel']=ic['trans']/ic['area']

eg={}
eg['trans']=daily['xport over 27.8'][sep:,:,:].where(daily.salinity<=srefb).sum('depth').sum('distance')
eg['fresha']=-(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
eg['freshb']=-(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
eg['freshc']=-(daily['xport over 27.8'][sep:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[sep:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
# eg['area']=(onesxr[6:,:-1,:].where(daily.salinity<=srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
# eg['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['den']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
# eg['meanvel']=eg['trans']/eg['area']

(cc['freshc']+egic['freshc']).sel(date=slice('2014-08-01','2016-04-30')).mean()
(cc['freshc']+egic['freshc']).sel(date=slice('2014-08-01','2016-04-30')).mean()/140

def decompose(cur,srefchoose):
    cur['abar']=nanmean(cur['area'])
    cur['vbar']=nanmean(cur['meanvel'])
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']
    cur['sbar']=-nanmean((cur['sal']-srefchoose)/srefchoose)
    cur['sprime']=(-(cur['sal']-srefchoose)/srefchoose-cur['sbar'])
    cur['vpAS']=(cur['sbar']*cur['abar']*cur['vprime']+cur['sbar']*cur['abar']*cur['vbar'])*1e3
    cur['VAsp']=(cur['sprime']*cur['abar']*cur['vbar']+cur['sbar']*cur['abar']*cur['vbar'])*1e3
    return cur

cc=decompose(cc,srefb)
egic=decompose(egic,srefb)


def savefilt(field):
    field['trans filt']=sig.filtfilt(B,A,field['trans'])
    field['fresh filt']=sig.filtfilt(B,A,field['freshb'])

    return field


egic=savefilt(egic)
eg=savefilt(eg)
ic=savefilt(ic)
cc=savefilt(cc)

# versname
#
# pickle.dump([cc,egic,eg,ic],open('../pickles/transdic_'+versname+'.pickle','wb'))


def colorstripes():
    axvspan(datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),color=egcol,alpha=0.4)


## OK, try fitting a sine func

from scipy.optimize import leastsq

def fitsin(t,data,guess_mean,guess_phase,guess_std,guess_period):

    first_guess=guess_std*np.sin(2*pi*(t+guess_phase)/guess_period) +guess_mean

    optimize_func = lambda x: x[0]*np.sin(2*pi*(t+x[1])/guess_period) + x[2] - data

    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

    # tgrid=linspace(t[0],t[-1],100)
    est_period=guess_period

    data_fit = est_std*np.sin(2*pi*(t+est_phase)/est_period) + est_mean

    # date=[]
    # for i in range(len(tgrid)):
    #     date.append(dt.date.fromtimestamp(tgrid[i]))

    return data_fit,est_std,est_period

tvec=arange(len(daily.date))

def implement_fitsin(field):
    figure(figsize=(9,9))
    # fit a sin to the filtered transport
    [trans_fit,trans_std,trans_period]=fitsin(tvec,field['trans filt'],mean(field['trans filt']),30,std(field['trans filt']),365.25)
    subplot(211)
    plot(daily.date,trans_fit)
    plot(daily.date,field['trans filt'])
    # print('trans period: ',trans_period)
    print('trans corr')
    print(corrcoef(trans_fit,field['trans filt'])[0,1])
    print('trans min')
    print(min(trans_fit))
    print(daily.date[argmin(trans_fit)].values)
    print('trans max')
    print(max(trans_fit))
    print(daily.date[argmax(trans_fit)].values)
    print('trans percent increase')
    print((max(trans_fit)-min(trans_fit))/min(trans_fit))
    # fit a sin to the filtered freshwater
    [fresh_fit,fresh_std,fresh_period]=fitsin(tvec,field['fresh filt'],mean(field['fresh filt']),30,std(field['fresh filt']),365.25)
    subplot(212)
    plot(daily.date,fresh_fit)
    plot(daily.date,field['fresh filt'])
    # print('fresh period: ',fresh_period)
    print('fresh corr')
    print(corrcoef(fresh_fit,field['fresh filt'])[0,1])
    print('fresh min')
    print(min(fresh_fit))
    print(daily.date[argmin(fresh_fit)].values)
    print('fresh max')
    print(max(fresh_fit))
    print(daily.date[argmax(fresh_fit)].values)
    print('fresh percent increase')
    print((max(fresh_fit)-min(fresh_fit))/min(fresh_fit))
    return trans_fit,fresh_fit

cc['trans seas'],cc['fresh seas']=implement_fitsin(cc)
(1.1-0.6)/1.1


egic['trans seas'],egic['fresh seas']=implement_fitsin(egic)
eg['trans seas'],eg['fresh seas']=implement_fitsin(eg)
ic['trans seas'],ic['fresh seas']=implement_fitsin(ic)

#################################################################################
######## Make 3 panel (freshwater) transport decomp for each current
####################################################################################

def eachpanel(field,colo,axit,labit='',xr=daily,pnofilt=1,letlab='',ls='-'):
    if pnofilt==1:
        axit.plot(xr.date,field,alpha=0.5,color=colo,label='',linewidth=0.75)
    axit.plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit,linestyle=ls)
    axit.set_xlabel('')
    axit.set_ylabel('')
    axit.text(0.01, 0.85,letlab,transform = axit.transAxes,fontsize=15)

years=matplotlib.dates.YearLocator()
months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')


def plottrans():
    f,((ax1), (ax2)) = plt.subplots(1, 2, sharex=True, figsize=(15,3.5))
    #1
    eachpanel(cc['trans'],ccol,ax1,labit='Coastal current')
    legend()
    ax1.set_ylim(0,2)
    ax1.set_ylabel('Transport [Sv]',fontsize=14)
    # ax1.set_title('Coastal current',fontsize=18)
    #2
    eachpanel(egic['trans'],egicol,ax2,labit='Slope current (EGIC)')
    eachpanel(eg['trans'],egcol,ax2,pnofilt=0,labit='East Greenland Current, S < 34.9')
    eachpanel(ic['trans'],icol,ax2,pnofilt=0,labit='Irminger current, S > 34.9')
    legend()
    # ax2.set_title('Slope Current',fontsize=18)
    ax2.set_ylim([0,35])
    tight_layout()
    savefig('../figures/paperfigs/voltrans.pdf',bbox_inches='tight')

plottrans()


seascol='k'

def plot3pan_coastal():
    f, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True, figsize=(8,7))
    ##2
    # eachpanel(cc['trans'],ccol,ax1)
    eachpanel(cc['trans'],ccol,ax1,letlab='a)')
    eachpanel(cc['trans seas'],seascol,ax1,pnofilt=0,ls='--')
    ax1.set_ylim(-0.1,2)

    ##2
    # eachpanel(cc['freshb'],ccol,ax2,labit='Between moorings')
    eachpanel(cc['freshb'],ccol,ax2,labit='Including shelf',letlab='b)')
    eachpanel(cc['fresh seas'],seascol,ax2,pnofilt=0,ls='--')
    # ax2.legend(loc=(0.27,0.95),framealpha=1)
    ax2.set_ylim(-10,120)
    ##3
    eachpanel(cc['freshb'],ccol,ax3,labit='vAF',pnofilt=0,letlab='c)')
    eachpanel(cc['vpAS'],'orange',ax3,pnofilt=0,labit='$\overline{v} \overline{A} \overline{F} + v\' \overline{A} \overline{F}$',ls='--')
    eachpanel(cc['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{v} \overline{A} \overline{F} + \overline{v} \overline{A}F\'$',ls='--')
    ax3.legend()
    ax3.set_ylim(20,70)
    fs=14
    xlab=-0.01
    f.text(xlab, 0.85, 'Transport [Sv]', va='center', rotation='vertical',fontsize=fs)
    f.text(xlab, 0.4, 'Freshwater transport [mSv]', va='center', rotation='vertical',fontsize=fs)
    plt.tight_layout()
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax3.xaxis.set_major_locator(years)
    ax3.xaxis.set_minor_locator(threemonth)
    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    # f.autofmt_xdate()
    ax1.set_title('Coastal current',fontsize=fs+4)
    savefig('../figures/paperfigs/coastal_tseries.pdf',bbox_inches='tight')

plot3pan_coastal()

std(egic['trans'])

def plot3pan_slope():
    f, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True, figsize=(8,7))
    ##2
    eachpanel(egic['trans'],egicol,ax1)
    eachpanel(eg['trans'],egcol,ax1,pnofilt=0)
    eachpanel(ic['trans'],icol,ax1,pnofilt=0,letlab='a)')
    eachpanel(egic['trans seas'],seascol,ax1,pnofilt=0,ls='--')
    ax1.set_ylim(0,30)
    ##2
    eachpanel(egic['freshb'],egicol,ax2)
    eachpanel(eg['freshb'],egcol,ax2,pnofilt=0,labit='East Greenland Current (S < 34.9)')
    eachpanel(ic['freshb'],icol,ax2,pnofilt=0,labit='Irminger Current (S > 34.9)',letlab='b)')
    eachpanel(egic['fresh seas'],seascol,ax2,pnofilt=0,ls='--')
    ax2.axhline(0,color='k')
    ax2.legend(loc=(0.3,0.95),framealpha=1)
    ax2.set_ylim([-35,100])
    ##3
    eachpanel(egic['freshb'],egicol,ax3,labit='vAF',pnofilt=0)
    eachpanel(egic['vpAS'],'orange',ax3,pnofilt=0,labit='$\overline{v} \overline{A} \overline{F} + v\' \overline{A} \overline{F}$',ls='--')
    eachpanel(egic['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{v} \overline{A} \overline{F} + \overline{v} \overline{A} F\'$',letlab='c)',ls='--')
    ax3.legend(loc=1)
    ax3.set_ylim(-25,80)
    fs=14
    xlab=-0.01
    f.text(xlab, 0.85, 'Transport [Sv]', va='center', rotation='vertical',fontsize=fs)
    f.text(xlab, 0.4, 'Freshwater transport [mSv]', va='center', rotation='vertical',fontsize=fs)
    plt.tight_layout()
    ax3.xaxis.set_major_locator(years)
    ax3.xaxis.set_minor_locator(threemonth)
    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax1.set_title('Slope current',fontsize=fs+4)
    savefig('../figures/paperfigs/slope_tseries.pdf',bbox_inches='tight')

plot3pan_slope()

##################################################################################
######################### PRES PLOTS ###########################################
##################################################################################

def pwf(field,colo,nofilt=0,labit='',xr=daily):
    if nofilt==0:
        plot(xr.date,field,alpha=0.5,color=colo)
    plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit)


def psf(field,colo,ylim1,ylim2,tit,nofilt=0,colcol=0,xr=daily):
    if colcol==1:
        colorstripes()
    pwf(field,colo,nofilt,xr=xr)
    ylim([ylim1,ylim2])
    xlabel('')



def OSMtrans(savename,ylab):
    ylabel(ylab)
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/'+savename+'.pdf',bbox_inches='tight')
    # savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/'+savename+'.pdf',bbox_inches='tight')



def plotccpres1():
    figure(figsize=(10,6))
    ax1=subplot(211)
    psf(cc['trans'],ccol,0,2,'cc_trans')
    gca().set_yticks(arange(0,2.5,0.5))
    ylim(0,2.1)
    text(datetime.datetime(2015,4,1),1.5,'0.86 $\pm$ 0.1 Sv',fontsize=20,color=ccol)
    plot([datetime.datetime(2016,8,1),datetime.datetime(2016,8,1)],[0.5,2],color='black',linewidth=5)
    OSMtrans('cctrans','Transport [Sv]')
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_minor_locator(threemonth)
    ax1.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    gca().set_xticklabels('')
    ax2=subplot(212)
    psf(cc['freshb'],ccol,0,100,'cc_fresh')
    plot([datetime.datetime(2016,8,1),datetime.datetime(2016,8,1)],[50,70],color='black',linewidth=5)
    text(datetime.datetime(2015,4,1),80,'42 $\pm$ 6 mSv',fontsize=20,color=ccol)
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    plt.tight_layout()
    OSMtrans('ccfresh','Freshwater transport [mSv]')

plotccpres1()

def egicpres1():
    figure(figsize=(10,6))
    ax1=subplot(211)
    psf(egic['trans'],'purple',5,30,'egic_trans')
    gca().set_yticks(arange(5,35,10))
    text(datetime.datetime(2016,1,1),7,'18 $\pm$ 1 Sv',fontsize=20,color='purple')
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_minor_locator(threemonth)
    ax1.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    gca().set_xticklabels('')
    OSMtrans('egictrans','Transport [Sv]')
    ax2=subplot(212)
    psf(egic['freshb'],'purple',-35,80,'egic_fresh')
    text(datetime.datetime(2015,1,1),-20,'14 $\pm$ 7 mSv',fontsize=20,color='purple')
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    plt.tight_layout()
    OSMtrans('egicfresh','Freshwater transport [mSv]')


egicpres1()

def plotccpres2():
    fig=figure(figsize=(10,6))
    ax1=subplot(211)
    psf(cc['trans'],ccol,0,2,'cc_trans')
    gca().set_yticks(arange(0,2.5,0.5))
    # text(datetime.datetime(2015,4,1),-1.2,'-0.6 $\pm$ 0.3 Sv',fontsize=20,color=ccol)
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_minor_locator(threemonth)
    ax1.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    gca().set_xticklabels('')
    OSMtrans('cctrans','Transport [Sv]')
    ax2=subplot(212)
    pwf(cc['freshb'],ccol,1,labit='total')
    pwf((cc['sbar']*cc['abar']*cc['vprime']+cc['sbar']*cc['abar']*cc['vbar'])*1e3,'grey',1,labit='velocity variations')
    pwf((cc['sprime']*cc['abar']*cc['vbar']+cc['sbar']*cc['abar']*cc['vbar'])*1e3,'k',1,labit='salinity variations')
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    plt.tight_layout()
    legend()
    (OSMtrans('ccfresh_decomp','Freshwater transport [mSv]'))

plotccpres2()

def egicpres2():
    fig=figure(figsize=(10,6))
    ax1=subplot(211)
    psf(egic['trans'],'purple',5,30,'egic_trans')
    gca().set_yticks(arange(5,35,10))
    # text(datetime.datetime(2016,1,1),-10,'21.6 $\pm$ 5.0 Sv',fontsize=20,color='purple')
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_minor_locator(threemonth)
    ax1.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    gca().set_xticklabels('')
    OSMtrans('egictrans','Transport [Sv]')
    ax2=subplot(212)
    pwf(egic['freshb'],egicol,1,labit='total')
    pwf((egic['sbar']*egic['abar']*egic['vprime']+egic['sbar']*egic['abar']*egic['vbar'])*1e3,'grey',1,labit='velocity variations')
    pwf((egic['sprime']*egic['abar']*egic['vbar']+egic['sbar']*egic['abar']*egic['vbar'])*1e3,'k',1,labit='salinity variations')
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    plt.tight_layout()
    legend()
    (OSMtrans('egicfresh_decomp','Freshwater transport [mSv]'))

egicpres2()

##################################################################################
######################### END PRES PLOTS ######################################
##################################################################################

def compfresh():
    fig=figure(figsize=(10,3.5))
    ax=subplot(111)
    freshtype='freshb'
    # fill_between(daily.date.values, sig.filtfilt(B,A, cc[freshtype]), sig.filtfilt(B,A,cc[freshtype+' plus']),color=ccol,label='',alpha=0.6)
    eachpanel(cc[freshtype],ccol,ax,pnofilt=0,labit='EGCC (coastal current)')
    eachpanel(eg[freshtype],egcol,ax,pnofilt=0,labit='EGC (slope current, S<34.9)')
    eachpanel(ic[freshtype],icol,ax,pnofilt=0,labit='IC (slope current, S>34.9)')
    axhline(0,color='k')
    # axvline(datetime.datetime(2015,7,1,0),color='grey',linewidth=0.8)
    # pwf(cc[freshtype]+egic[freshtype],'grey',1)
    ylim([-30,120])
    ylabel('Freshwater transport [mSv]')
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
    colorstripes()
    # fig.autofmt_xdate()
    # dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]
    # gca().set_xticks(dticks)
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_minor_locator(threemonth)
    ax.xaxis.set_minor_formatter(monthFMT)
    ax.xaxis.set_major_formatter(yearFMT)
    legend(loc=1)
    labyy=90
    text(datetime.datetime(2014,10,15),labyy,'FALL',color=ccol,fontsize=13.5)
    text(datetime.datetime(2015,1,10),labyy,'WINTER',color=egcol,fontsize=13.5)
    text(datetime.datetime(2015,4,18),labyy,'SPRING/SUMMER',color='grey',fontsize=13.5)
    # text(datetime.datetime(2015,7,5),75,'SUMMER',color='grey',fontsize=15)
    savefig('../figures/paperfigs/freshcomp.pdf',bbox_inches='tight')
    savefig('../figures/paperfigs/freshcomp.jpg',bbox_inches='tight')

compfresh()


# def comptrans():
#     fig=figure(figsize=(10,3.5))
#     ax=subplot(111)
#     freshtype='trans'
#     # fill_between(daily.date.values, sig.filtfilt(B,A, cc[freshtype]), sig.filtfilt(B,A,cc[freshtype+' plus']),color=ccol,label='',alpha=0.6)
#     eachpanel(cc[freshtype],ccol,ax,pnofilt=1,labit='Coastal current')
#     eachpanel(eg[freshtype],egcol,ax,pnofilt=1,labit='Slope current < 34.9')
#     eachpanel(ic[freshtype],icol,ax,pnofilt=1,labit='Slope current > 34.9')
#     axhline(0,color='k')
#     # axvline(datetime.datetime(2015,7,1,0),color='grey',linewidth=0.8)
#     # pwf(cc[freshtype]+egic[freshtype],'grey',1)
#     # ylim([-30,100])
#     ylabel('Transport [mSv]')
#     xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
#     colorstripes()
#     fig.autofmt_xdate()
#     dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]
#     gca().set_xticks(dticks)
#     # legend(loc=(1.05,0.5),framealpha=1)


# def comptrans():
#     fig=figure(figsize=(10,3.5))
#     ax=subplot(111)
#     freshtype='trans'
#     # fill_between(daily.date.values, sig.filtfilt(B,A, cc[freshtype]), sig.filtfilt(B,A,cc[freshtype+' plus']),color=ccol,label='',alpha=0.6)
#     eachpanel(cc[freshtype],ccol,ax,pnofilt=1,labit='Coastal current')
#     eachpanel(eg[freshtype],egcol,ax,pnofilt=1,labit='Slope current < 34.9')
#     eachpanel(ic[freshtype],icol,ax,pnofilt=1,labit='Slope current > 34.9')
#     axhline(0,color='k')
#     # axvline(datetime.datetime(2015,7,1,0),color='grey',linewidth=0.8)
#     # pwf(cc[freshtype]+egic[freshtype],'grey',1)
#     # ylim([-30,100])
#     ylabel('Transport [mSv]')
#     xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
#     colorstripes()
#     fig.autofmt_xdate()
#     dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]
#     gca().set_xticks(dticks)
#     # legend(loc=(1.05,0.5),framealpha=1)
#     legend()
#     labyy=82
#     legend()
#     labyy=82
#     # text(datetime.datetime(2014,10,15),labyy,'FALL',color=ccol,fontsize=13.5)
#     # text(datetime.datetime(2015,1,10),labyy,'WINTER',color=egcol,fontsize=13.5)
#     # text(datetime.datetime(2015,4,3),labyy,'SPRING/SUMMER',color='grey',fontsize=13.5)
#     # text(datetime.datetime(2015,7,5),75,'SUMMER',color='grey',fontsize=15)
#     savefig('../figures/paperfigs/transcomp.pdf',bbox_inches='tight')
#
# comptrans()


#
# def refcomp():
#     f, ((ax1), (ax2), (ax3), (ax4)) = plt.subplots(4, 1, sharex=True, figsize=(8,9))
#     eachpanel(cc['fresha'],'b',ax1,pnofilt=0,labit='34.8')
#     eachpanel(cc['freshb'],'r',ax1,pnofilt=0,labit='34.9')
#     eachpanel(cc['freshc'],'g',ax1,pnofilt=0,labit='35')
#     ax1.set_title('Coastal Current')
#
#     eachpanel(eg['fresha'],'b',ax2,pnofilt=0,labit='34.8')
#     eachpanel(eg['freshb'],'r',ax2,pnofilt=0,labit='34.9')
#     eachpanel(eg['freshc'],'g',ax2,pnofilt=0,labit='35')
#     ax2.legend()
#     ax2.set_title('Slope Current < 34.9')
#
#     eachpanel(ic['fresha'],'b',ax3,pnofilt=0,labit='34.8')
#     eachpanel(ic['freshb'],'r',ax3,pnofilt=0,labit='34.9')
#     eachpanel(ic['freshc'],'g',ax3,pnofilt=0,labit='35')
#     ax3.set_title('Slope Current > 34.9')
#     ax3.set_ylabel('Freshwater transport [mSv]')
#
#     eachpanel(egic['fresha']+cc['fresha'],'b',ax4,pnofilt=0,labit='34.8')
#     eachpanel(egic['freshb']+cc['freshb'],'r',ax4,pnofilt=0,labit='34.9')
#     eachpanel(egic['freshc']+cc['freshc'],'g',ax4,pnofilt=0,labit='35')
#     ax4.set_title('Coastal + Slope')
#
#
# refcomp()
# savefig('../notes/section_details/transfigs/refsalcomp.pdf',bbox_inches='tight')
#


#################################################################################
#### TS diagrams
#################################################################################
def TSplot():
        hexbin(daily['salinity'].where(daily.distance>0).values.flatten(),daily['temperature'].where(daily.distance>0).values.flatten(),cmap=cm.Blues,bins='log',gridsize=200,rasterized=True,mincnt=1,vmin=-1)
        colorbar(label='[log # measurements]')
        xlim([32,35.2])
        ylim([-2,10])
        xlabel('Salinity')
        plot([32.5,32.5,34.25,34.25,32.5],[-1.8,0,0,-1.8,-1.8],color='k')
        # den=contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(24,28,0.5))
        # manual_locations=[(32,2),(33,2),(33,2),(34,2)]
        # clabel(den,fmt='%1.1f')#,manual=manual_locations)
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
    # ax2.hexbin(daily['salinity'].where(daily.distance>0).values.flatten(),daily['temperature'].where(daily.distance>0).values.flatten(),
    # cmap=cm.Greys,bins='log',gridsize=200,rasterized='True')
    axpart=ax2.hexbin(daily['salinity'].sel(date=slice(date1,date2)).values.flatten(),daily['temperature'].sel(date=slice(date1,date2)).values.flatten(),
    cmap=cm.Blues,bins='log',mincnt=1,gridsize=200,rasterized='True',vmin=-1)
    ax2.axvline(34.9,color='k')
    ax2.plot(daily['salinity'].where(daily.distance>0).sel(date=slice(date1,date2)).sel(distance=slice(45,62)).mean(dim='date').mean(dim='distance'),
    daily['temperature'].where(daily.distance>0).sel(date=slice(date1,date2)).sel(distance=slice(45,62)).mean(dim='date').mean(dim='distance'),
    linewidth=3,color='red')
    ax2.plot(daily['salinity'].where(daily.depth<150).sel(date=slice(date1,date2)).sel(distance=slice(-15,13)).mean(dim='date').mean(dim='distance'),
    daily['temperature'].where(daily.depth<150).sel(date=slice(date1,date2)).sel(distance=slice(-15,13)).mean(dim='date').mean(dim='distance'),
    linewidth=3,color='purple')
    ax2.set_xlim([31,35.2])
    ax2.set_ylim([-2,8])
    ax2.plot([32.5,32.5,34.25,34.25,32.5],[-1.8,0,0,-1.8,-1.8],color='k')
    # ax2.contour(salvec,tmpvec,pdenmat,colors='grey',zorder=50,levels=arange(24,28,0.5),alpha=0.5)
    # ax2.text(34.9495,6,'AW',color='k',fontsize=16,zorder=60)
    # ax2.text(34.3,-1,'PW',color='k',fontsize=16,zorder=60)

    ax2.text(31.5,6,textit,fontsize=18)

    return axpart

# def plot4TS_seas(vv=0):
#     f, ((ax11, ax22), (ax33, ax44)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(9.5,7))
#     partTS('2014-10-01','2015-1-1','FALL',ax11,vv)
#     partTS('2015-01-01','2015-04-1','WINTER',ax22,vv)
#     partTS('2015-04-01','2015-7-1','SPRING',ax33,vv)
#     axts=partTS('2015-7-01','2015-10-1','SUMMER',ax44,vv)
#     f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=16)
#     f.text(-0.03, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=16)
#     plt.tight_layout()
#     cbaxes = f.add_axes([1.01, 0.15, 0.025, 0.7])
#     # if vv==1:
#     #     cbar=colorbar(axts,label='Velocity in top 500m [m/s]',cax=cbaxes)
#     #     savefig('../figures/paperfigs/TS_seasonal_vel.pdf',bbox_inches='tight')
#     if vv==0:
#         cbar=colorbar(axts,label='[log # measurements]',cax=cbaxes)
#         savefig('../figures/paperfigs/TS_seasonal.pdf',bbox_inches='tight')
#
# # plot4TS_seas(1)
#
# plot4TS_seas(0)

def plot4TS_seas_symm(vv=0):
    f, (ax11, ax22, ax33) = plt.subplots(3, 1, sharex=True, sharey=True,figsize=(4.5,9))
    partTS('2014-10-01','2015-1-1','FALL',ax11,vv)
    partTS('2015-01-01','2015-04-1','WINTER',ax22,vv)
    ax11.text(32.03,-1.4,'PW',color='k',fontsize=16)
    ax22.text(32.03,-1.4,'PW',color='k',fontsize=16)
    ax33.text(32.03,-1.4,'PW',color='k',fontsize=16)
    ax11.text(32.5,3,'shelf',color='purple',fontsize=17)
    ax11.text(34.1,6.5,'slope',color='red',fontsize=17,)
    ax22.text(34.4,5.5,'AW',color='k',fontsize=16)
    axts=partTS('2015-04-01','2015-10-1','SPRING/SUMMER',ax33,vv)
    f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=16)
    plt.tight_layout()
    cbaxes = f.add_axes([1.01, 0.3, 0.04, 0.4])
    cbar=colorbar(axts,label='[log # measurements]',cax=cbaxes, drawedges=False)
    savefig('../figures/paperfigs/TS_seasonal_symm.pdf',bbox_inches='tight',)

plot4TS_seas_symm()

def TSwinter():
    f,axi,=subplots(1,1,figsize=(4.5,3.5))
    partTS('2015-01-01','2015-04-1','',axi)
    axi.text(32.5,4.5,'WINTER',fontsize=18)
    xlim(32,35.05)
    ylim(-2,6)
    xlabel('Salinity')
    ylabel('Potential temperature [$^\circ$C]')
    text(32.75,0.5,'shelf',color='purple',fontsize=17)
    text(34.1,4.5,'slope',color='red',fontsize=17,)
    savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/TSwinter.pdf',bbox_inches='tight')
TSwinter()

#############################################################
############## Current vel profiles  ######################
#############################################################

# CTD=pickle.load(open('../pickles/Shipboard/CTD_xarray_1806_relmoor.pickle','rb'))

def pltprf(dd,date1,date2,labit,field,col,d1,d2):
    if dd==0:
        plot(daily[field].sel(date=slice(date1,date2)).sel(distance=0).mean(dim='date'),daily.depth,label=labit,linewidth=3,color=col)
    else:
        # plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').T,daily.depth,'.',alpha=0.3,color=col)
        plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='distance').mean(dim='date'),daily.depth,label=labit,linewidth=3,color=col)
        plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='distance').mean(dim='date').isel(depth=50),100,'k*',label='',markersize=10,zorder=20)
        plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='distance').mean(dim='date').isel(depth=150),300,'ko',label='',zorder=20)
    # else:
    #     skipna=10
    #     plot(daily[field].sel(date=slice(date1,date2)).sel(distance=45).mean(dim='date')[::skipna],daily.depth[::skipna],label=labit,linewidth=3,color=col)


def pltseas(dd,field,d1,d2):
    pltprf(dd,'2014-10-01','2015-1-1','Fall',field,ccol,d1,d2)
    pltprf(dd,'2015-1-01','2015-4-1','Winter',field,egcol,d1,d2)
    pltprf(dd,'2015-4-01','2015-10-1','Spring/Summer',field,'orange',d1,d2)



def pltTS(date1,date2,labit,col,d1,d2):#
    p1=int(0/2)
    p2=int(1500/2)
    plot(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).isel(depth=slice(p1,p2)).mean(dim='date').mean(dim='distance'),
    daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).isel(depth=slice(p1,p2)).mean(dim='date').mean(dim='distance'),label=labit,linewidth=3,color=col)
    plot(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=150),
    daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=150),'ko',label='',zorder=20)
    plot(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50),
    daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50),'k*',markersize=10,label='',zorder=20)
    print(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50))
    print(daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50))

def pltseas_TScomp():
    ctd1=45
    ctd2=62
    ctdp=100
    ctda=0.5
    figure(figsize=(8,3.5))
    subplot(121)
    pltseas(1,'temperature',ctd1,ctd2)
    ylim(500,0)
    xlabel('Potential temperature [$^\circ$C]')
    ylabel('depth [m]')
    alf=0.5
    axhline(100,color='k',alpha=alf)
    axhline(300,color='k',alpha=alf)
    xlim(3,6.5)
    subplot(122)
    pltseas(1,'salinity',ctd1,ctd2)
    # plot(CTD.salinity.sel(occupation='2014 (KN221)').sel(distance=45).sel(pressure=slice(0,ctdp)),
    #         CTD.pressure.sel(pressure=slice(0,ctdp)),'.-',color='orange',alpha=ctda,label='Aug 2014 profiles',markersize=10,zorder=10)
    # plot(CTD.salinity.sel(occupation='2014 (KN221)').sel(distance=slice(ctd1,ctd2)).sel(pressure=slice(0,ctdp)),
    #         CTD.pressure.sel(pressure=slice(0,ctdp)),'.-',color='orange',alpha=ctda,label='',markersize=10,zorder=10)
    xlabel('Salinity')
    axhline(100,color='k',alpha=alf)
    axhline(300,color='k',alpha=alf)
    legend(loc=(-0.3,-0.7),fontsize=12,markerscale=2)
    gca().set_yticklabels('')
    ylim([500,0])
    tight_layout()
    savefig('../figures/paperfigs/TSprofs_seasonal_slopeave.pdf',bbox_inches='tight',)
    figure()
    pltTS('2014-10-01','2015-1-1','Fall',ccol,ctd1,ctd2)
    pltTS('2015-1-01','2015-4-1','Winter',egcol,ctd1,ctd2)
    pltTS('2015-4-01','2015-10-1','Spring/Summer','orange',ctd1,ctd2)
    arrow(34.885, 6, -0.16, -1.6, head_width=0.025, head_length=0.25, fc='k', ec='k',zorder=20,linewidth=2)
    arrow(34.69, 3.99, 0.13, 0.77, head_width=0.025, head_length=0.2, fc='k', ec='k',zorder=20,linewidth=2)
    # plot(CTD.salinity.sel(occupation='2014 (KN221)').sel(distance=45).sel(pressure=slice(0,ctdp)),
    # CTD.temperature.sel(occupation='2014 (KN221)').sel(distance=45).sel(pressure=slice(0,ctdp)),'g.-',label='August 2014 CTD', alpha=0.5)
    # plot(CTD.salinity.sel(occupation='2016').sel(distance=slice(ctd1,ctd2)).sel(pressure=slice(0,ctdp)),
    # CTD.temperature.sel(occupation='2016').sel(distance=slice(ctd1,ctd2)).sel(pressure=slice(0,ctdp)),'g*-',label='July 2016 CTD', alpha=0.5)
    contour(salvec,tmpvec,pdenmat,colors='grey',zorder=10,levels=arange(27.3,27.6,0.01),alpha=0.3)
    cc=contour(salvec,tmpvec,pdenmat,colors='grey',zorder=10,levels=arange(27.3,27.6,0.1))
    xlim(34.5,35.0)
    ylim(3,6.5)
    manual_locations=[(34.55,4),(34.6,3.2),(34.72,3.2),(34.83,3.4)]
    clabel(cc,fmt='%1.1f',manual=manual_locations)
    axvline(34.9,color='k')
    xlabel('Salinity')
    ylabel('Potential temperature [$^\circ$C]')
    savefig('../figures/paperfigs/TS_seasonal_slopeave.pdf',bbox_inches='tight',)

pltseas_TScomp()

def pltseas_TScomp_yr2():
    ctd1=45
    ctd2=62
    ctdp=100
    ctda=0.5
    figure(figsize=(8,3.5))
    subplot(121)
    pltseas(1,'temperature',ctd1,ctd2)
    ylim(500,0)
    xlabel('Potential temperature [$^\circ$C]')
    ylabel('depth [m]')
    alf=0.5
    axhline(100,color='k',alpha=alf)
    axhline(300,color='k',alpha=alf)
    xlim(3,6.5)
    subplot(122)
    pltseas(1,'salinity',ctd1,ctd2)
    # plot(CTD.salinity.sel(occupation='2014 (KN221)').sel(distance=45).sel(pressure=slice(0,ctdp)),
    #         CTD.pressure.sel(pressure=slice(0,ctdp)),'.-',color='orange',alpha=ctda,label='Aug 2014 profiles',markersize=10,zorder=10)
    # plot(CTD.salinity.sel(occupation='2014 (KN221)').sel(distance=slice(ctd1,ctd2)).sel(pressure=slice(0,ctdp)),
    #         CTD.pressure.sel(pressure=slice(0,ctdp)),'.-',color='orange',alpha=ctda,label='',markersize=10,zorder=10)
    xlabel('Salinity')
    axhline(100,color='k',alpha=alf)
    axhline(300,color='k',alpha=alf)
    legend(loc=(-0.3,-0.7),fontsize=12,markerscale=2)
    gca().set_yticklabels('')
    ylim([500,0])
    tight_layout()
    savefig('../figures/paperfigs/TSprofs_seasonal_slopeave_yr2.pdf',bbox_inches='tight',)
    figure()
    pltTS('2015-9-01','2015-12-1','Fall',ccol,ctd1,ctd2)
    pltTS('2015-12-01','2016-3-1','Winter',egcol,ctd1,ctd2)
    pltTS('2016-3-01','2016-8-1','Spring/Summer','orange',ctd1,ctd2)
    arrow(34.885, 6, -0.17, -1.7, head_width=0.015, head_length=0.2, fc='k', ec='k',zorder=20)
    arrow(34.69, 3.99, 0.13, 0.7, head_width=0.02, head_length=0.15, fc='k', ec='k',zorder=20)
    # plot(CTD.salinity.sel(occupation='2014 (KN221)').sel(distance=45).sel(pressure=slice(0,ctdp)),
    # CTD.temperature.sel(occupation='2014 (KN221)').sel(distance=45).sel(pressure=slice(0,ctdp)),'g.-',label='August 2014 CTD', alpha=0.5)
    # plot(CTD.salinity.sel(occupation='2016').sel(distance=slice(ctd1,ctd2)).sel(pressure=slice(0,ctdp)),
    # CTD.temperature.sel(occupation='2016').sel(distance=slice(ctd1,ctd2)).sel(pressure=slice(0,ctdp)),'g*-',label='July 2016 CTD', alpha=0.5)
    contour(salvec,tmpvec,pdenmat,colors='grey',zorder=10,levels=arange(27.3,27.6,0.01),alpha=0.3)
    cc=contour(salvec,tmpvec,pdenmat,colors='grey',zorder=10,levels=arange(27.3,27.6,0.1))
    xlim(34.5,35.0)
    ylim(3,6.5)
    manual_locations=[(34.55,4),(34.6,3.2),(34.72,3.2),(34.83,3.4)]
    clabel(cc,fmt='%1.1f',manual=manual_locations)
    axvline(34.9,color='k')
    xlabel('Salinity')
    ylabel('Potential temperature [$^\circ$C]')
    savefig('../figures/paperfigs/TS_seasonal_slopeave_yr2.pdf',bbox_inches='tight',)

pltseas_TScomp_yr2()

def pltprf_grd(dd,date1,date2,labit,field,col,d1,d2,axx):
    if dd==0:
        axx.plot(daily[field].sel(date=slice(date1,date2)).sel(distance=0).mean(dim='date'),daily.depth,label=labit,linewidth=3,color=col)
    else:
        axx.plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='distance').mean(dim='date'),daily.depth,label=labit,linewidth=3,color=col)
        axx.plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='distance').mean(dim='date').isel(depth=50),100,'k*',label='',markersize=10,zorder=20)
        axx.plot(daily[field].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='distance').mean(dim='date').isel(depth=150),300,'ko',label='',zorder=20)


def pltseas_grd(dd,field,d1,d2,axx):
    pltprf_grd(dd,'2014-10-01','2015-1-1','Fall',field,ccol,d1,d2,axx)
    pltprf_grd(dd,'2015-1-01','2015-4-1','Winter',field,egcol,d1,d2,axx)
    pltprf_grd(dd,'2015-4-01','2015-10-1','Spring/Summer',field,'orange',d1,d2,axx)

def pltTS_grd(date1,date2,labit,col,d1,d2,axx):#
    p1=int(0/2)
    p2=int(1500/2)
    axx.plot(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).isel(depth=slice(p1,p2)).mean(dim='date').mean(dim='distance'),
    daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).isel(depth=slice(p1,p2)).mean(dim='date').mean(dim='distance'),label=labit,linewidth=3,color=col)
    axx.plot(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=150),
    daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=150),'ko',label='',zorder=20)
    axx.plot(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50),
    daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50),'k*',markersize=10,label='',zorder=20)
    # print(daily['salinity'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50))
    # print(daily['temperature'].sel(date=slice(date1,date2)).sel(distance=slice(d1,d2)).mean(dim='date').mean(dim='distance').isel(depth=50))


def pltseas_TScomp_gridspec():
    ctd1=45
    ctd2=62
    ctdp=100
    ctda=0.5
    fig=plt.figure(figsize=(7,9))
    # subplot(121)
    gs = gridspec.GridSpec(30, 12)
    ax1=fig.add_subplot(gs[:10,:6])
    pltseas_grd(1,'temperature',ctd1,ctd2,ax1)
    ylim(500,0)
    xlabel('Potential temperature [$^\circ$C]')
    ylabel('depth [m]')
    alf=0.5
    axhline(100,color='k',alpha=alf)
    axhline(300,color='k',alpha=alf)
    xlim(3,6.5)
    # subplot(122)
    ax2=fig.add_subplot(gs[:10,6:])
    pltseas_grd(1,'salinity',ctd1,ctd2,ax2)
    xlabel('Salinity')
    axhline(100,color='k',alpha=alf)
    axhline(300,color='k',alpha=alf)
    gca().set_yticklabels('')
    ylim([500,0])
    # tight_layout()
    # axleg=fig.add_subplot(gs[9:11,5:7])
    legend(fontsize=12,markerscale=2,loc=(-0.3,-0.65),)
    # savefig('../figures/paperfigs/TSprofs_seasonal_slopeave.pdf',bbox_inches='tight',)
    # figure()
    ax3=fig.add_subplot(gs[17:,1:11])
    pltTS_grd('2014-10-01','2015-1-1','Fall',ccol,ctd1,ctd2,ax3)
    pltTS_grd('2015-1-01','2015-4-1','Winter',egcol,ctd1,ctd2,ax3)
    pltTS_grd('2015-4-01','2015-10-1','Spring/Summer','orange',ctd1,ctd2,ax3)
    arrow(34.9, 5.9, -0.15, -1.5, head_width=0.025, head_length=0.25, fc='k', ec='k',zorder=20,linewidth=2)
    arrow(34.725, 4.05, 0.13, 0.77, head_width=0.025, head_length=0.2, fc='k', ec='k',zorder=20,linewidth=2)
    contour(salvec,tmpvec,pdenmat,colors='grey',zorder=10,levels=arange(27.4,27.7,0.01),alpha=0.3)
    cc=contour(salvec,tmpvec,pdenmat,colors='grey',zorder=10,levels=arange(27.4,27.72,0.1))
    xlim(34.5,35.0)
    ylim(3,6.5)
    manual_locations=[(34.55,4),(34.55,3.2),(34.7,3.4),(34.83,3.2)]
    clabel(cc,fmt='%1.1f',manual=manual_locations)
    axvline(34.9,color='k')
    xlabel('Salinity')
    ylabel('Potential temperature [$^\circ$C]')
    savefig('../figures/paperfigs/Slopeprops_seasonal.pdf',bbox_inches='tight',)


pltseas_TScomp_gridspec()
XXXXXXXXXXXX
# def pltprfs():
#     figure(figsize=(6,12))
#     subplot(311)
#     pltseas(1,'potential density')
#     ylim([1300,0])
#     subplot(312)
#     pltseas(1,'temperature')
#     ylim([1300,0])
#     legend()
#     subplot(313)
#     pltseas(1,'salinity')
#     ylim([1300,0])
#
# pltprfs()


#
# def velprfcmp():
#     fs=14
#     f=figure(figsize=(7,8))
#     subplot(221)
#     pltseas(0,'across track velocity')
#     ylim([180,0])
#     f.text(-0.03, 0.82, 'Coastal current', rotation='vertical',fontsize=fs+1)
#     text(-0.88, 15, 'a)',fontsize=fs+1)
#     xlim(-0.9,0)
#     subplot(223)
#     pltseas(1,'across track velocity')
#     ylim([1100,0])
#     f.text(0, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fs)
#     text(-0.88, 100, 'c)',fontsize=fs+1)
#     xlim(-0.9,0)
#     xlabel('velocity [m/s]',fontsize=fs)
#     f.text(-0.03, 0.35, 'Slope current', rotation='vertical',fontsize=fs+1)
#     subplot(222)
#     pltseas(0,'salinity')
#     text(32.05, 15, 'b)',fontsize=fs+1)
#     xlim(32,34.2)
#     ylim([180,0])
#     gca().set_yticklabels('')
#     subplot(224)
#     pltseas(1,'salinity')
#     text(34.37, 100, 'd)',fontsize=fs+1)
#     xlim(34.35,35)
#     gca().set_yticklabels('')
#     legend(loc=(1.05,0.9),fontsize=fs-1)
#     ylim([1100,0])
#     xlabel('salinity',fontsize=fs)
#     savefig('../figures/paperfigs/Seas_velandsalprofs.pdf',bbox_inches='tight')
#
#
# velprfcmp()




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

autocorr(cc['trans'])



len(daily.date)/10
#effective degrees of freedom
Neff=70
ttest=1.99 #t value for 90 degrees of freedome and 95% confidence
def getconf(field):
    conf=ttest*std(field)/sqrt(Neff)
    return conf

#############################################################
###################  Report means and error  ###############
#############################################################
cc['trans'].mean()
cc['freshb'].mean()
getconf(cc['trans'])
getconf(cc['freshb'])


egic['trans'].mean()
egic['freshb'].mean()

getconf(egic['trans'])
getconf(egic['freshb'])

eg['trans'].mean()
eg['freshb'].mean()

getconf(eg['trans'])
getconf(eg['freshb'])


ic['trans'].mean()
ic['freshb'].mean()

getconf(ic['trans'])
getconf(ic['freshb'])

egic.keys()

# get July/August estimates
ic['trans'].groupby('date.month').mean(dim='date').sel(month=slice(7,8)).mean()
ic['trans'].groupby('date.month').std(dim='date').sel(month=slice(7,8)).mean()

eg['trans'].groupby('date.month').mean(dim='date').sel(month=slice(7,8)).mean()
eg['trans'].groupby('date.month').std(dim='date').sel(month=slice(7,8)).mean()












XXXXXXXXXXXXXXXX
# END
#
# dlen=len(daily.date)
# int(dlen/2)
# ## Do maxima for EGCC, EGC and EGIC (should all abide by similar forms)
# def getmaxmin(field,absopt=1):
#     print('MAX 1')
#     dind=int(dlen/2)
#     if absopt==0:
#         print('freshwater no abs val')
#         print(max(field['fresh filt'][:dind]))
#         print(daily.date[argmax(field['fresh filt'][:dind])].values)
#
#     print('volume')
#     print(max(abs(field['trans filt'][:dind])))
#     print(daily.date[argmax(abs(field['trans filt'][:dind]))].values)
#     print('freshwater')
#     print(max(abs(field['fresh filt'][:dind])))
#     print(daily.date[argmax(abs(field['fresh filt'][:dind]))].values)
#
#     print('MIN')
#     mind1=int(dlen/4)
#     mind2=int(3*dlen/4)
#     if absopt==0:
#         print('freshwater no abs val')
#         print(min(field['fresh filt'][mind1:mind2]))
#         print(daily.date[mind1:mind2][argmin(field['fresh filt'][mind1:mind2])].values)
#
#     print('volume')
#     print(min(abs(field['trans filt'][mind1:mind2])))
#     print(daily.date[mind1:mind2][argmin(abs(field['trans filt'][mind1:mind2]))].values)
#     print('freshwater')
#     print(min(abs(field['fresh filt'][mind1:mind2])))
#     print(daily.date[mind1:mind2][argmin(abs(field['fresh filt'][mind1:mind2]))].values)
#
#     print('MAX 2')
#     if absopt==0:
#         print('freshwater no abs val')
#         print(max(field['fresh filt'][dind:]))
#         print(daily.date[dind:][argmax(field['fresh filt'][dind:])].values)
#
#     print('volume')
#     print(max(abs(field['trans filt'][dind:])))
#     print(daily.date[dind:][argmax(abs(field['trans filt'][dind:]))].values)
#     print('freshwater')
#     print(max(abs(field['fresh filt'][dind:])))
#     print(daily.date[dind:][argmax(abs(field['fresh filt'][dind:]))].values)
#
#
# getmaxmin(cc)
#
# getmaxmin(egic,0)
#
# getmaxmin(eg)
#
# getmaxmin(ic,0)

# #### Keeping this legacy stuff around because I might have known something. Reconcile depending on what I decide/ streamline scripts
# print('first max')
# ## First max
# mdate1=datetime.datetime(2014,8,1)
# mdate2=datetime.datetime(2015,11,1)
# print(max(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2))))))
# print(daily.date.sel(date=slice(mdate1,mdate2))[argmax(abs(sig.filtfilt(B,A,eg['trans'].sel(date=slice(mdate1,mdate2)))))].values)
# print('second max')
# ## Second ma
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
# ################################################
# #### Coastal section
# ################################################
# for key in cc:
#     if ('trans' in key) | ('fresh' in key):
#             print(key)
#             print('mean=',cc[key].mean().values)
#             print('std=',cc[key].std().values)
#
# date1=datetime.datetime(2014,8,15)
# ## Print max and min (freshwater) transports
# middate=datetime.datetime(2015,5,1,0)
# for key in cc:
#     if ('trans' in key) | ('fresh' in key):
#             print(key)
#             # print('mean=',cc[key].mean().values)
#             # print('std=',cc[key].std().values)
#             print(max(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(date1,middate))))))
#             print(daily.date[argmax(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(date1,middate)))))].values)
#
#
# lastdate=datetime.datetime(2016,7,30)
# for key in cc:
#     if ('trans' in key) | ('fresh' in key):
#         print(key)
#         print(max(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(middate,lastdate))))))
#         print(daily.date.sel(date=slice(middate,lastdate))[argmax(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(middate,lastdate)))))].values)
#
# mdate1=datetime.datetime(2014,10,1)
# mdate2=datetime.datetime(2015,9,1)
# for key in cc:
#     if ('trans' in key) | ('fresh' in key):
#         print(key)
#         print(min(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(mdate1,mdate2))))))
#         print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(mdate1,mdate2)))))].values)
#
#
#
# ## Print (freshwater) transport for first two weeks of August
# date2=datetime.datetime(2014,9,1)
# for key in cc:
#     if ('trans' in key) | ('fresh' in key):
#         print(key)
#         print(cc[key].sel(date=slice(date1,date2)).mean().values)
