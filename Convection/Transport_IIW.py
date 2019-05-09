#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

versname='1810JHIL'
daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_'+versname+'.pickle','rb'))

daily['across track velocity']=-1*daily['across track velocity']

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

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

onesxr=daily.salinity/daily.salinity


srefb=34.9
sep=9

egic={}
egic['trans']=daily['xport over 27.8'][sep:,:,:].sum('depth').sum('distance')
egic['area']=(onesxr.where(daily['potential density']<27.8)[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport over 27.8'][sep:,:-1,:]*daily['salinity'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport over 27.8'][sep:,:-1,:]*daily['temperature'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport over 27.8'][sep:,:-1,:]*daily['potential density'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']

for ii in range(0,200,20):
    figure()
    (daily['potential density'][sep:,:-1,:].where((daily['potential density']<27.74)&(daily['potential density']>=27.68)).isel(date=ii).T).plot()
    (daily['potential density'][sep:,:-1,:].where((daily['potential density']<27.8)&(daily['potential density']>=27.74)).isel(date=ii).T).plot()


uIIW={}
uIIW['trans']=daily.xport.where((daily['potential density']<27.74)&(daily['potential density']>=27.68)).sum('distance').sum('depth')

dIIW={}
dIIW['trans']=daily.xport.where((daily['potential density']<27.8)&(daily['potential density']>=27.74)).sum('distance').sum('depth')

IIW={}
IIW['trans']=daily.xport.where((daily['potential density']<27.8)&(daily['potential density']>=27.68)).sum('distance').sum('depth')


def savefilt(field):
    field['trans filt']=sig.filtfilt(B,A,field['trans'])

    return field


egic=savefilt(egic)
uIIW=savefilt(uIIW)
dIIW=savefilt(dIIW)
IIW=savefilt(IIW)

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
    f,ax3=subplots(1,1,figsize=(8,2.5),sharex=True)
    # eachpanel(egic['trans'],egicol,ax1,labit='Slope current transport')
    # eachpanel(IIW['trans'],'k',ax1,labit='Total IIW transport')
    eachpanel(uIIW['trans'],'k',ax3,labit='upper IIW transport')
    eachpanel(dIIW['trans'],'C7',ax3,labit='deep IIW transport')
    # ax1.legend(loc=2)
    ax3.legend(loc=2)
    ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax3.xaxis.set_major_locator(years)
    ax3.xaxis.set_minor_locator(threemonth)
    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    ax3.set_yticks(arange(0,15,4))
    # ax1.set_ylabel('[Sv]')
    ax3.set_ylabel('[Sv]')
    # ax1.set_ylim(5,35)
    savefig(figdir+'MixedLayer/pres/IIW_trans.pdf',bbox_inches='tight')


plottrans()

How much is vel, how much is layer thickness!?

XXXXXXXXXXXXXXX
# def implement_fitsin(field):
#     figure(figsize=(9,9))
#     # fit a sin to the filtered transport
#     [trans_fit,trans_std,trans_period]=fitsin(tvec,field['trans filt'],mean(field['trans filt']),30,std(field['trans filt']),365.25)
#     subplot(211)
#     plot(daily.date,trans_fit)
#     plot(daily.date,field['trans filt'])
#     # print('trans period: ',trans_period)
#     print('trans corr')
#     print(corrcoef(trans_fit,field['trans filt'])[0,1])
#     print('trans min')
#     print(min(trans_fit))
#     print(daily.date[argmin(trans_fit)].values)
#     print('trans max')
#     print(max(trans_fit))
#     print(daily.date[argmax(trans_fit)].values)
#     print('trans percent increase')
#     print((max(trans_fit)-min(trans_fit))/min(trans_fit))
#     # fit a sin to the filtered freshwater
#     [fresh_fit,fresh_std,fresh_period]=fitsin(tvec,field['fresh filt'],mean(field['fresh filt']),30,std(field['fresh filt']),365.25)
#     subplot(212)
#     plot(daily.date,fresh_fit)
#     plot(daily.date,field['fresh filt'])
#     # print('fresh period: ',fresh_period)
#     print('fresh corr')
#     print(corrcoef(fresh_fit,field['fresh filt'])[0,1])
#     print('fresh min')
#     print(min(fresh_fit))
#     print(daily.date[argmin(fresh_fit)].values)
#     print('fresh max')
#     print(max(fresh_fit))
#     print(daily.date[argmax(fresh_fit)].values)
#     print('fresh percent increase')
#     print((max(fresh_fit)-min(fresh_fit))/min(fresh_fit))
#     return trans_fit,fresh_fit
#
# cc['trans seas'],cc['fresh seas']=implement_fitsin(cc)
# egic['trans seas'],egic['fresh seas']=implement_fitsin(egic)
# eg['trans seas'],eg['fresh seas']=implement_fitsin(eg)
# ic['trans seas'],ic['fresh seas']=implement_fitsin(ic)
