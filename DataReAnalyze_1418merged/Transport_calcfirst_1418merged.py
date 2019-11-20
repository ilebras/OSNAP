#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from firstfuncs_1618 import *

daily=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid.nc')
daily['across track velocity']=-1*daily['across track velocity']
daily=daily.where(daily['across track velocity']!=0)
daily=daily.sel(date=slice('2014-8-15','2018-9-1')) #just to avoid some zeros at the beginning which mess up filtering.


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
    cur['vbar']=nanmean(~isinf(cur['meanvel']))
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']
    cur['sbar']=-nanmean((cur['sal']-srefchoose)/srefchoose)
    cur['sprime']=(-(cur['sal']-srefchoose)/srefchoose-cur['sbar'])
    cur['vpAS']=(cur['sbar']*cur['abar']*cur['vprime']+cur['sbar']*cur['abar']*cur['vbar'])*1e3
    cur['VAsp']=(cur['sprime']*cur['abar']*cur['vbar']+cur['sbar']*cur['abar']*cur['vbar'])*1e3
    return cur

cc=decompose(cc,srefb)

cc['VAsp'].plot()

egic=decompose(egic,srefb)

B,A = sig.butter(2,1./50, output='ba')

def savefilt(field):
    field['trans filt']=sig.filtfilt(B,A,field['trans'])
    field['fresh filt']=sig.filtfilt(B,A,field['freshb'])

    return field


egic=savefilt(egic)
eg=savefilt(eg)
ic=savefilt(ic)
cc=savefilt(cc)


def colorstripes():
    axvspan(datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2016,10,1),datetime.datetime(2017,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2017,10,1),datetime.datetime(2018,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2017,1,1),datetime.datetime(2017,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2018,1,1),datetime.datetime(2018,4,1),color=egcol,alpha=0.4)


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


plottrans()

seascol='k'

def plot2pan_coastal():
    f, ((ax1),  (ax3)) = plt.subplots(2, 1, sharex=True, figsize=(12,5))
    ##2
    # eachpanel(cc['trans'],ccol,ax1)
    eachpanel(cc['trans'],ccol,ax1,letlab='a)')
    eachpanel(cc['trans seas'],seascol,ax1,pnofilt=0,ls='--')
    ax1.set_ylim(-0.1,2)

    ##2
    # eachpanel(cc['freshb'],ccol,ax2,labit='Between moorings')
    eachpanel(cc['freshb'],ccol,ax3,labit='Including shelf',letlab='b)')
    eachpanel(cc['fresh seas'],seascol,ax3,pnofilt=0,ls='--')
    # ax2.legend(loc=(0.27,0.95),framealpha=1)
    ax3.set_ylim(-10,120)
    # ##3
    # eachpanel(cc['freshb'],ccol,ax3,labit='vAF',pnofilt=0,letlab='c)')
    # eachpanel(cc['vpAS'],'orange',ax3,pnofilt=1,labit='$\overline{v} \overline{A} \overline{F} + v\' \overline{A} \overline{F}$',ls='--')
    # eachpanel(cc['VAsp'],'grey',ax3,pnofilt=1,labit='$\overline{v} \overline{A} \overline{F} + \overline{v} \overline{A}F\'$',ls='--')
    # ax3.legend()
    # ax3.set_ylim(20,70)
    fs=14
    xlab=-0.01
    f.text(xlab, 0.75, 'Transport [Sv]', va='center', rotation='vertical',fontsize=fs)
    f.text(xlab, 0.25, 'Freshwater transport [mSv]', va='center', rotation='vertical',fontsize=fs)
    plt.tight_layout()
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2018,8,15)])
    ax3.xaxis.set_major_locator(years)
    ax3.xaxis.set_minor_locator(threemonth)
    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    # f.autofmt_xdate()
    ax1.set_title('Coastal current',fontsize=fs+4)
    savefig(figdir+'ReAnalyze/egc_fresh/CoastalCurrent_volfreshtrans.png',bbox_inches='tight')

plot2pan_coastal()

def plot2pan_slope():
    f, ((ax1),  (ax3)) = plt.subplots(2, 1, sharex=True, figsize=(12,5))
    ##2
    eachpanel(egic['trans'],egicol,ax1)
    eachpanel(eg['trans'],egcol,ax1,pnofilt=0)
    eachpanel(ic['trans'],icol,ax1,pnofilt=0,letlab='a)')
    eachpanel(egic['trans seas'],seascol,ax1,pnofilt=0,ls='--')
    ax1.set_ylim(0,30)
    ##2
    eachpanel(egic['freshb'],egicol,ax3)
    eachpanel(eg['freshb'],egcol,ax3,pnofilt=0,labit='East Greenland Current (S < 34.9)')
    eachpanel(ic['freshb'],icol,ax3,pnofilt=0,labit='Irminger Current (S > 34.9)',letlab='b)')
    eachpanel(egic['fresh seas'],seascol,ax3,pnofilt=0,ls='--')
    ax3.axhline(0,color='k')
    ax3.legend(loc=(0.3,1.05),framealpha=1)
    ax3.set_ylim([-35,100])
    # ##3
    # eachpanel(egic['freshb'],egicol,ax3,labit='vAF',pnofilt=0)
    # eachpanel(egic['vpAS'],'orange',ax3,pnofilt=0,labit='$\overline{v} \overline{A} \overline{F} + v\' \overline{A} \overline{F}$',ls='--')
    # eachpanel(egic['VAsp'],'grey',ax3,pnofilt=0,labit='$\overline{v} \overline{A} \overline{F} + \overline{v} \overline{A} F\'$',letlab='c)',ls='--')
    # ax3.legend(loc=1)
    # ax3.set_ylim(-25,80)
    fs=14
    xlab=-0.01
    f.text(xlab, 0.75, 'Transport [Sv]', va='center', rotation='vertical',fontsize=fs)
    f.text(xlab, 0.25, 'Freshwater transport [mSv]', va='center', rotation='vertical',fontsize=fs)
    plt.tight_layout()
    ax3.xaxis.set_major_locator(years)
    ax3.xaxis.set_minor_locator(threemonth)
    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2018,8,15)])
    ax1.set_title('Slope current',fontsize=fs+4)
    savefig(figdir+'ReAnalyze/egc_fresh/SlopeCurrent_volfreshtrans.png',bbox_inches='tight')

plot2pan_slope()

def compfresh():
    fig=figure(figsize=(14,3.5))
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
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2018,8,15)])
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_minor_locator(threemonth)
    ax.xaxis.set_minor_formatter(monthFMT)
    ax.xaxis.set_major_formatter(yearFMT)
    legend(loc=(0.01,0.7))
    labyy=90
    # text(datetime.datetime(2014,10,15),labyy,'FALL',color=ccol,fontsize=13.5)
    # text(datetime.datetime(2015,1,10),labyy,'WINTER',color=egcol,fontsize=13.5)
    # text(datetime.datetime(2015,4,18),labyy,'SPRING/SUMMER',color='grey',fontsize=13.5)
    # text(datetime.datetime(2015,7,5),75,'SUMMER',color='grey',fontsize=15)
    savefig(figdir+'ReAnalyze/egc_fresh/Freshcomp.png',bbox_inches='tight')


compfresh()


#################################################################################
#### TS diagrams
#################################################################################
def TSplot():
        hexbin(daily['salinity'].where(daily.distance>0).values.flatten(),daily['temperature'].where(daily.distance>0).values.flatten(),cmap=cm.Blues,bins='log',gridsize=200,rasterized=True,mincnt=1)
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
        savefig(figdir+'ReAnalyze/egc_fresh/TSall.png',bbox_inches='tight')


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
    cmap=cm.Blues,bins='log',mincnt=1,gridsize=200,rasterized='True')
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


def plot4TS_seas_symm(vv=0):
    for ii in range(1,5):
        f, (ax11, ax22, ax33) = plt.subplots(3, 1, sharex=True, sharey=True,figsize=(4.5,9))
        title('Year '+str(ii))
        partTS(str(2013+ii)+'-10-01',str(2014+ii)+'-1-1','FALL',ax11,vv)
        partTS(str(2014+ii)+'-01-01',str(2014+ii)+'-04-1','WINTER',ax22,vv)
        ax11.text(32.03,-1.4,'PW',color='k',fontsize=16)
        ax22.text(32.03,-1.4,'PW',color='k',fontsize=16)
        ax33.text(32.03,-1.4,'PW',color='k',fontsize=16)
        ax11.text(32.5,3,'shelf',color='purple',fontsize=17)
        ax11.text(34.1,6.5,'slope',color='red',fontsize=17,)
        ax22.text(34.4,5.5,'AW',color='k',fontsize=16)
        axts=partTS(str(2014+ii)+'-04-01',str(2014+ii)+'-10-1','SPRING/SUMMER',ax33,vv)
        f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=16)
        f.text(-0.03, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=16)
        plt.tight_layout()
        cbaxes = f.add_axes([1.01, 0.3, 0.04, 0.4])
        cbar=colorbar(axts,label='[log # measurements]',cax=cbaxes, drawedges=False)
        savefig(figdir+'ReAnalyze/egc_fresh/TSseas_year'+str(ii)+'.png',bbox_inches='tight')


plot4TS_seas_symm()
