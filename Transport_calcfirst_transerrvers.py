#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1803bathy.pickle','rb'))

# mask the fields based on bathymetry file (this version does not have extra fields on either side, thats just for plotting)


[bathatmoor,newdpth]=pickle.load(open('../pickles/moordpths.pickle','rb'))

daily=newgrid.copy()
for vv in daily:
    if vv[0]!='d':
        print(vv)
        for dd,adist in enumerate(newgrid.distance):
            daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=newdpth[dd])

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

# figure()
# plot(bathdist,bathbath,'-')
# plot(newgrid.distance,bathonmygrid,'o')

# hexbin(daily.salinity.values.flatten(),daily.temperature.values.flatten(),C=daily['across track velocity'].values.flatten(),cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4)
# axvline(34.8,color='k')
# axvline(34.9,color='k')
# axvline(35,color='k')
# colorbar()
# savefig('../figures/xport/TSall_velcont.png')

# hexbin(daily.salinity.mean(dim='date').values.flatten(),daily.temperature.mean(dim='date').values.flatten(),C=daily['across track velocity'].mean('date').values.flatten(),cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4)
# axvline(34.8,color='k')
# axvline(34.9,color='k')
# axvline(35,color='k')
# colorbar()
# savefig('../figures/xport/TSmean_velcont.png')

#################################################################################
###################### (Freshwater) Transport #################################
#################################################################################

mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))



# Adding back minind to document the difference
## surface def
minind=zeros(len(daily.date))
for tt,na in enumerate(daily.date):
    minind[tt]=where(max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt])[0][0]+2

minind=[int(mm) for mm in minind]

mindxr=daily.date.copy()
mindxr.values=daily.distance[minind]




ccvel=daily['across track velocity'].copy()
for tt,mm in enumerate(minind):
        ccvel[mm:,:,tt]=NaN

srefa=34.8
srefb=34.9
srefc=35

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3

onesxr=daily.salinity/daily.salinity


cc={}
cc['trans']=daily.xport[:6,:,:].sum('depth').sum('distance')
cc['trans sal def']=(daily['across track velocity'].where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cc['trans plus sal']=(daily['across track velocity'].where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3)[:6,:,:].sum('depth').sum('distance')
cc['trans min vel']=(ccvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cc['trans min vel plus sal']=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

cc['fresha']=(daily['xport'][:6,:-1,:]*1e3*(daily.salinity[:6,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
cc['freshb']=(daily['xport'][:6,:-1,:]*1e3*(daily.salinity[:6,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
cc['freshc']=(daily['xport'][:6,:-1,:]*1e3*(daily.salinity[:6,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
cc['area']=(onesxr[:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('depth').sum('distance')
cc['sal']=(daily['xport'][:6,:-1,:]*daily['salinity'][:6,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['tmp']=(daily['xport'][:6,:-1,:]*daily['temperature'][:6,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['den']=(daily['xport'][:6,:-1,:]*daily['potential density'][:6,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['meanvel']=cc['trans']/cc['area']

egic={}
egic['trans']=daily.xport[6:,:,:].sum('depth').sum('distance')
egic['fresha']=(daily['xport'][6:,:-1,:]*1e3*(daily.salinity[6:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
egic['freshb']=(daily['xport'][6:,:-1,:]*1e3*(daily.salinity[6:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
egic['freshc']=(daily['xport'][6:,:-1,:]*1e3*(daily.salinity[6:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
egic['area']=(onesxr[6:,:-1,:]*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport'][6:,:-1,:]*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport'][6:,:-1,:]*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport'][6:,:-1,:]*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']


ic={}
ic['trans']=daily.xport[6:,:,:].where(daily.salinity>srefb).sum('depth').sum('distance')
ic['fresha']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[6:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
ic['freshb']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[6:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
ic['freshc']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*1e3*(daily.salinity[6:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
ic['area']=(onesxr[6:,:-1,:].where(daily.salinity>srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
ic['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
ic['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
ic['den']=(daily['xport'][6:,:-1,:].where(daily.salinity>srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/ic['trans']
ic['meanvel']=ic['trans']/ic['area']


eg={}
eg['trans']=daily.xport[6:,:,:].where(daily.salinity<=srefb).sum('depth').sum('distance')
eg['fresha']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[6:,:-1,:]-srefa)/srefa).sum('depth').sum('distance')
eg['freshb']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[6:,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
eg['freshc']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*1e3*(daily.salinity[6:,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
eg['area']=(onesxr[6:,:-1,:].where(daily.salinity<=srefb)*depthdiffmat[6:,:,:]*middistmat[6:,:,:]/1e3).sum('depth').sum('distance')
eg['sal']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['salinity'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
eg['tmp']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['temperature'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
eg['den']=(daily['xport'][6:,:-1,:].where(daily.salinity<=srefb)*daily['potential density'][6:,:-1,:]).sum('distance').sum('depth')/eg['trans']
eg['meanvel']=eg['trans']/eg['area']

mean(eg['freshb'])+mean(cc['freshb'])

def colorstripes():
    axvspan(datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,9,1),datetime.datetime(2015,12,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2015,12,1),datetime.datetime(2016,3,1),color=egcol,alpha=0.4)
    # axvspan(datetime.datetime(2015,6,1),datetime.datetime(2015,9,1),color=icol,alpha=0.4)
    # axvspan(datetime.datetime(2016,5,1),datetime.datetime(2016,8,1),color=icol,alpha=0.4)

import scipy.signal as sig
# Design the Buterworth filter
N  = 2    # Filter order
Wn = 0.02 # Cutoff frequency
B, A = sig.butter(N, Wn, output='ba')



def pwf(field,colo,nofilt=0,labit='',xr=daily):
    if nofilt==0:
        field.plot(alpha=0.5,color=colo,label='')
    plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit)


def psf(field,colo,ylim1,ylim2,tit,nofilt=0,colcol=0,xr=daily,labit=''):
    if colcol==1:
        colorstripes()
    pwf(field,colo,nofilt,xr=xr,labit=labit)
    ylim([ylim1,ylim2])
    xlabel('')
    # if nofilt==0:
    #     savefig('../figures/xport/'+tit+'_nocol.png')



def OSMtrans(savename,ylab):
    ylabel(ylab)
    savefig('../figures/OSM_post/'+savename+'.pdf',bbox_inches='tight')

1/0.02

figure(figsize=(10,3))
psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
psf(cc['trans plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
legend()
gca().set_yticks(arange(-1.5,0.3,0.5))
ylabel('Transport [Sv]')
savefig('../figures/xport/cc_trans_fixedpos_salcomp.pdf',bbox_inches='tight')

figure(figsize=(10,3))
psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
psf(cc['trans min vel plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
legend()
gca().set_yticks(arange(-1.5,0.3,0.5))
ylabel('Transport [Sv]')
savefig('../figures/xport/cc_trans_minvel_salcomp.pdf',bbox_inches='tight')


figure(figsize=(10,9))
subplot(311)
psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
psf(cc['trans sal def'],'orange',-1.5,0,'cc_trans',labit='sal def')
legend()
gca().set_yticks(arange(-1.5,0.3,0.5))
gca().set_xticklabels('')
ylabel('Transport [Sv]')
subplot(312)
psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
psf(cc['trans min vel plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
legend()
gca().set_yticks(arange(-1.5,0.3,0.5))
ylabel('Transport [Sv]')
gca().set_xticklabels('')
subplot(313)
psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
psf(cc['trans plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
legend()
gca().set_yticks(arange(-1.5,0.3,0.5))
ylabel('Transport [Sv]')
savefig('../figures/xport/cc_trans_allcomp.pdf',bbox_inches='tight')


figure(figsize=(10,3))
psf(mindxr,'purple',0,40,'position of surface velocity minimum')
ylabel('')
title('position of surface velocity minimum between currents [km]')
savefig('../figures/xport/minpos_smooth.pdf',bbox_inches='tight')

figure(figsize=(10,6))
subplot(211)
psf(cc['trans'],ccol,-1.5,0,'cc_trans')
gca().set_yticks(arange(-1.5,0.3,0.5))
text(datetime.datetime(2015,4,1),-1.2,'-0.6 $\pm$ 0.3 Sv',fontsize=20,color=ccol)
gca().set_xticklabels('')
OSMtrans('cctrans','Transport [Sv]')
subplot(212)
psf(cc['freshb'],ccol,0,70,'cc_fresh')
text(datetime.datetime(2015,4,1),50,'26 $\pm$ 14 mSv',fontsize=20,color=ccol)
OSMtrans('ccfresh','Freshwater flux [mSv]')


figure(figsize=(10,6))
subplot(211)
psf(egic['trans'],'purple',-35,-5,'egic_trans')
gca().set_yticks(arange(-30,-5,10))
text(datetime.datetime(2016,1,1),-10,'-21.6 $\pm$ 5.0 Sv',fontsize=20,color='purple')
gca().set_xticklabels('')
OSMtrans('egictrans','Transport [Sv]')
subplot(212)
psf(egic['freshb'],'purple',-35,80,'egic_fresh')
text(datetime.datetime(2015,1,1),-20,'13 $\pm$ 27 mSv',fontsize=20,color='purple')
OSMtrans('egicfresh','Freshwater flux [mSv]')

def decompose(cur,srefchoose):
    cur['abar']=nanmean(cur['area'])
    cur['vbar']=nanmean(cur['meanvel'])
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']

    cur['sbar']=nanmean((cur['sal']-srefchoose)/srefchoose)
    cur['sprime']=((cur['sal']-srefchoose)/srefchoose-cur['sbar'])

    return cur


cc=decompose(cc,srefb)
egic=decompose(egic,srefb)

egic['trans'].mean()
fig=figure(figsize=(10,6))
subplot(211)
psf(cc['trans'],ccol,-1.5,0,'cc_trans')
gca().set_yticks(arange(-1.5,0.3,0.5))
text(datetime.datetime(2015,4,1),-1.2,'-0.6 $\pm$ 0.3 Sv',fontsize=20,color=ccol)
gca().set_xticklabels('')
OSMtrans('cctrans','Transport [Sv]')
subplot(212)
pwf(cc['freshb'],ccol,1,labit='$vAS$')
pwf((cc['sbar']*cc['abar']*cc['vprime']+cc['sbar']*cc['abar']*cc['vbar'])*1e3,'grey',1,labit='$v\'\overline{A}\overline{S}$')
pwf((cc['sprime']*cc['abar']*cc['vbar']+cc['sbar']*cc['abar']*cc['vbar'])*1e3,'k',1,labit='$\overline{v}\overline{A}S\'$')
fig.autofmt_xdate()
legend()
(OSMtrans('ccfresh_decomp','Freshwater flux [mSv]'))

fig=figure(figsize=(10,6))
subplot(211)
psf(egic['trans'],'purple',-35,-5,'egic_trans')
gca().set_yticks(arange(-30,-5,10))
text(datetime.datetime(2016,1,1),-10,'-21.6 $\pm$ 5.0 Sv',fontsize=20,color='purple')
gca().set_xticklabels('')
OSMtrans('egictrans','Transport [Sv]')
subplot(212)
pwf(egic['freshb'],egicol,1,labit='$vAS$')
pwf((egic['sbar']*egic['abar']*egic['vprime']+egic['sbar']*egic['abar']*egic['vbar'])*1e3,'grey',1,labit='$v\'\overline{A}\overline{S}$')
pwf((egic['sprime']*egic['abar']*egic['vbar']+egic['sbar']*egic['abar']*egic['vbar'])*1e3,'k',1,labit='$\overline{v}\overline{A}S\'$')
fig.autofmt_xdate()
legend()
(OSMtrans('egicfresh_decomp','Freshwater flux [mSv]'))



def compfresh(freshtype,ylim1,ylim2,fresh=1):
    fig=figure(figsize=(10,3.5))
    pwf(cc[freshtype],ccol,1,'Coastal Current')
    pwf(eg[freshtype],egcol,1,'Slope Current < 34.9')
    pwf(ic[freshtype],icol,1,'Slope Current > 34.9')
    axhline(0,color='k')
    # pwf(cc[freshtype]+egic[freshtype],'grey',1)
    ylim([ylim1,ylim2])
    colorstripes()
    fig.autofmt_xdate()

    if fresh==1:
        ylabel('Freshwater flux [mSv]')
    else:
        ylabel('Transport [Sv]')
    savefig('../figures/xport/'+freshtype+'_all_comps.png')

dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]

compfresh('freshb',-35,80)
gca().set_xticks(dticks)
text(datetime.datetime(2014,10,15),68,'FALL',color=ccol,fontsize=15)
legend()
text(datetime.datetime(2015,1,7),68,'WINTER',color=egcol,fontsize=15)
savefig('../figures/OSM_post/FreshB_all_comps.pdf',bbox_inches='tight')

#### TS diagram
def TSplot():
        hexbin(daily['salinity'].values.flatten(),daily['temperature'].values.flatten(),cmap=cm.hot_r,bins='log')
        colorbar(label='[log # measurements]')
        xlim([32,35.2])
        ylim([-2,10])
        xlabel('salinity')
        plot([32.5,32.5,34.25,34.25,32.5],[-1.8,0,0,-1.8,-1.8],color='k')
        den=contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(24,28,0.5))
        manual_locations=[(32,2),(33,2),(33,2),(34,2)]
        clabel(den,fmt='%1.1f')#,manual=manual_locations)
        text(34.95,7,'AW',color='k',fontsize=20)
        text(34.1,0.1,'PW',color='k',fontsize=20)
        ylabel('pot. temperature [$^\circ$C]')
        axvline(34.9,color='k')
        savefig('../figures/OSM_post/TSall.pdf',bbox_inches='tight')

TSplot()



psf(egic['trans'],'purple',-35,-5,'')
colorstripes()

psf(winddat['across flow wind'].mean(dim='distance').resample('D',dim='date'),'k',-10,5,'',xr=winddat.resample('D',dim='date'))
colorstripes()


psf(eg['trans'],egcol,-10,0,'eg_trans')
psf(ic['trans'],icol,-30,-10,'ic_trans')
psf(egic['trans'],egicol,-30,-10,'egic_trans')

psf(cc['fresh'],ccol,0,50,'cc_fresh')
psf(eg['freshb'],egcol,0,150,'eg_fresh')
psf(ic['fresh'],icol,-5,2,'ic_fresh')

psf(egic['fresh'],egicol,-5,150,'egic_fresh')

psf(egic['fresh'],egicol,-5,150,'egic_fresh',nofilt=1)
pwf(ic['fresh'],icol,nofilt=1)
pwf(eg['fresh'],egcol,nofilt=1)
pwf(cc['fresh'],ccol,nofilt=1)
pwf(cc['fresh']+egic['fresh'],'k',nofilt=1)
axhline(0,color='grey')
savefig('../figures/xport/fresh_all.png')










compfresh('freshc',0,80)

pwf(cc['trans'],ccol)
ylim([-1,0.2])
axvspan(datetime.datetime(2014,9,1),datetime.datetime(2015,2,1),color=ccol,alpha=0.4)
axvspan(datetime.datetime(2015,9,1),datetime.datetime(2016,2,1),color=ccol,alpha=0.4)

pwf(eg['trans'],egcol)
ylim([-10,0])
axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,5,1),color=egcol,alpha=0.4)
axvspan(datetime.datetime(2016,1,1),datetime.datetime(2016,5,1),color=egcol,alpha=0.4)


#################################################################################
###################### Decomposition #################################
#################################################################################


def decompose(cur,srefchoose):
    cur['abar']=nanmean(cur['area'])
    cur['vbar']=nanmean(cur['meanvel'])
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']

    cur['sbar']=nanmean((cur['sal']-srefchoose)/srefchoose)
    cur['sprime']=((cur['sal']-srefchoose)/srefchoose-cur['sbar'])

    return cur


def plot_decompfresh(cur,colcol,y1,y2,tit,freshchoose):
    figure(figsize=(18,10))
    subplot(321)
    plot(cur['vprime'].date,cur['sbar']*cur['abar']*cur['vprime']*1e3,color=colcol)
    title('$\overline{s} \overline{a} v\'$: '+str(int(corrcoef(cur['sbar']*cur['abar']*cur['vprime'],cur[freshchoose])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    ylabel('[mSv]')
    subplot(322)
    plot(cur['vprime'].date,cur['sbar']*cur['aprime']*cur['vbar']*1e3,color=colcol)
    title('$\overline{s} a\' \overline{v} $: '+str(int(corrcoef(cur['sbar']*cur['aprime']*cur['vbar'],cur[freshchoose])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    subplot(323)
    plot(cur['vprime'].date,cur['sprime']*cur['abar']*cur['vbar']*1e3,color=colcol)
    title('$s\' \overline{a} \overline{v} $: '+str(int(corrcoef(cur['sprime']*cur['abar']*cur['vbar'],cur[freshchoose])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    ylabel('[mSv]')
    subplot(324)
    plot(cur['vprime'].date,cur['sprime']*cur['abar']*cur['vprime']*1e3,color=colcol)
    title('$s\' \overline{a} v\'} $: '+str(int(corrcoef(cur['sprime']*cur['abar']*cur['vprime'],cur[freshchoose])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    subplot(325)
    remterms=cur['sprime']*cur['aprime']*cur['vbar']+cur['sprime']*cur['aprime']*cur['vprime']+cur['sbar']*cur['aprime']*cur['vprime']
    plot(cur['vprime'].date,(remterms)*1e3,color=colcol)
    ylim([y1,y2])
    title('remaining terms: '+str(int(corrcoef(remterms,cur[freshchoose])[1,0]*100)/100))
    ylabel('[mSv]')
    subplot(326)
    plot(cur['vprime'].date,cur[freshchoose]-cur['abar']*cur['sbar']*cur['vbar']*1e3,color=colcol)
    title('Total variability')
    ylim([y1,y2])
    savefig('../figures/xport/Decomp_fresh_'+tit+'.png')


def decomplot(current,curtit,colo,ylim1=-100,ylim2=100):
    current=decompose(current,srefa)
    plot_decompfresh(current,colo,ylim1,ylim2,curtit+'_fresha','fresha')
    current=decompose(current,srefb)
    plot_decompfresh(current,colo,ylim1,ylim2,curtit+'_freshb','freshb')
    current=decompose(current,srefc)
    plot_decompfresh(current,colo,ylim1,ylim2,curtit+'_freshc','freshc')

decomplot(egic,'egic',egicol)


decomplot(ic,'ic',icol)


decomplot(eg,'eg',egcol)


def plot_decomptrans(cur,y1,y2,colcol,tit):
    figure(figsize=(10,12))
    subplot(411)
    plot(cur['vprime'].date,cur['abar']*cur['vprime'],color=colcol)
    title('$\overline{a} v\'$: '+str(int(corrcoef(cur['abar']*cur['vprime'],cur['trans'])[1,0]*100)/100))
    gca().set_xticklabels('')
    ylabel('[Sv]')
    ylim([y1,y2])
    subplot(412)
    plot(cur['vprime'].date,cur['aprime']*cur['vbar'],color=colcol)
    title('$a\' \overline{v}$: '+str(int(corrcoef(cur['aprime']*cur['vbar'],cur['trans'])[1,0]*100)/100))
    gca().set_xticklabels('')
    ylabel('[Sv]')
    ylim([y1,y2])
    subplot(413)
    plot(cur['vprime'].date,cur['aprime']*cur['vprime'],color=colcol)
    title('$a\' v\'$: '+str(int(corrcoef(cur['aprime']*cur['vprime'],cur['trans'])[1,0]*1000)/1000))
    gca().set_xticklabels('')
    ylim([y1,y2])
    ylabel('[Sv]')
    subplot(414)
    plot(cur['vprime'].date,cur['abar']*cur['vprime']+cur['aprime']*cur['vbar']+cur['aprime']*cur['vprime'],color=colcol)
    title('Total transport variability')
    ylabel('[Sv]')
    ylim([y1,y2])
    savefig('../figures/newtrans/Decomp_trans_'+tit+'.png')



def curwmprop(cur,coll,tit):
    # figure(figsize=(10,8))
    subplot(311)
    plot(daily.date,cur['sal'],color=coll)
    axhline(mean(cur['sal']),color=coll)
    ylabel('salinity')
    title('Water mass properties in the '+tit)
    subplot(312)
    plot(daily.date,cur['tmp'],color=coll)
    axhline(mean(cur['tmp']),color=coll)
    ylabel('potential temperature [$^\circ$ C]')
    subplot(313)
    plot(daily.date,cur['den'],color=coll)
    axhline(mean(cur['den']),color=coll)
    ylabel('potential density [kg/$m^3$]')

figure(figsize=(10,8))
curwmprop(egrem,egcol,'East Greenland Current Remnant')
savefig('../figures/newtrans/WMprops_EGCrem.png',bbox_inches='tight')

figure(figsize=(10,8))
curwmprop(egic,egicol,'Slope Current')
savefig('../figures/newtrans/WMprops_EGIC.png',bbox_inches='tight')

figure(figsize=(10,8))
curwmprop(cc,ccol,'Coastal Current')
savefig('../figures/newtrans/WMprops_EGCC.png',bbox_inches='tight')

figure(figsize=(10,8))
curwmprop(egrem,egcol,'East Greenland Current Remnant')
curwmprop(egic,egicol,'Slope Current')
curwmprop(cc,ccol,'Coastal Current')
subplot(311)
title('Water mass properties in all current components')
savefig('../figures/newtrans/WMprops_all.png',bbox_inches='tight')

figure(figsize=(10,4))
(cctrans.resample('W',dim='date',how='std')/std(cctrans)).plot(color=ccol,label='coastal current')
(egictrans.resample('W',dim='date',how='std')/std(egictrans)).plot(color=egicol,label='slope current')
title('Normalized standard deviation of transport')
legend()
savefig('../figures/newtrans/Transvar_both.png',bbox_inches='tight')


fig, ax1 = plt.subplots(figsize=(12,3),)
egictrans.plot(alpha=0.5,ax=ax1,color=egicol)
egictrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=egicol,ax=ax1)
ax1.set_ylabel('East Greenland-Irminger Current',color=egicol)
ax1.tick_params('y', colors=egicol)
ax1.set_ylim([-25,-12])
ax2 = ax1.twinx()
cctrans.plot(alpha=0.5,ax=ax2,color=ccol)
cctrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax2)
ax2.set_title('Transports in the EGC system [Sv]')
ax2.set_ylabel('Coastal Current',color=ccol)
ax2.tick_params('y', colors=ccol)
ax2.set_ylim([-2,0])
savefig('../figures/newtrans/EGCC-EGIC_trans.png')
