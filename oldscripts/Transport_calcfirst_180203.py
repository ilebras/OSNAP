#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1801.pickle','rb'))

# mask the fields based on bathymetry file (this version does not have extra fields on either side, thats just for plotting)


bathf=interpolate.interp1d(bathdist,bathbath)
bathonmygrid=bathf(newgrid.distance)

# figure()
# plot(bathdist,bathbath,'-')
# plot(newgrid.distance,bathonmygrid,'o')

daily=newgrid.copy()
for vv in newgrid:
    if vv[0]!='d':
        print(vv)
        for dd,adist in enumerate(daily.distance):
            daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=bathonmygrid[dd])


# Establish min in abs vel to differentiate EGCC
minind=zeros(len(daily.date))
for tt,na in enumerate(daily.date):
    minind[tt]=where(max(daily['across track velocity'][2:15,:,tt].mean(dim='depth'))==daily['across track velocity'][2:15,:,tt].mean(dim='depth'))[0][0]+2

minind=[int(mm) for mm in minind]

mindxr=daily.date.copy()
mindxr.values=daily.distance[minind]



mindxr.plot()

figure(figsize=(10,4))
plot(daily.date,run_ave(mindxr.values,20))
title('distance at which minimum between currents occurs [km]')
ylabel('')
savefig('../figures/xport/minpos_smooth_dpthave.png')

ccvel=daily['across track velocity'].copy()
ccarea=daily['across track velocity'].copy()
ccarea[:]=1


for tt,mm in enumerate(minind):
        ccvel[mm:,:,tt]=NaN
        ccarea[mm:,:,tt]=NaN

ccsal=daily['salinity'].copy()
for tt,mm in enumerate(minind):
        ccsal[mm:,:,tt]=NaN


egvel=daily['across track velocity'].where(daily['potential density']<27.8)
egvel_fulld=daily['across track velocity'].copy()
egarea=daily['across track velocity'].copy()
egarea[:]=1
egarea=egarea.where(daily['potential density']<27.8)
for tt,mm in enumerate(minind):
        egvel[:mm,:,tt]=NaN
        egvel_fulld[:mm,:,tt]=NaN
        egarea[:mm,:,tt]=NaN

#################################################################################
###################### (Freshwater) Transport #################################
#################################################################################

mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

srefa=34
srefb=34.8
srefc=35


ccf={}
ccf['trans']=(ccvel.where(daily.salinity<33)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
ccf['fresh']=(ccvel.where(daily.salinity<33)[:,:-1,:]*(daily.salinity[:,:-1,:]-srefb)/srefb*depthdiffmat*middistmat).sum('depth').sum('distance')
onesxr=daily.salinity/daily.salinity

hexbin(daily.salinity,daily.temperature)

daily.plot(kind='hexbin',x='salinity',y='temperature')

cc={}
# cc['trans to5 wsal']=(daily['across track velocity'].where(daily.salinity<34)[:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('depth').sum('distance')
# cc['trans just sal']=(daily['across track velocity'].where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
# cc['trans moving minvel']=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance') #decided this was too noisy for unclear reasons, lets leave it out for now
# cc['fresh']=(ccvel.where(daily.salinity<34)[:,:-1,:]*(daily.salinity[:,:-1,:]-srefb)/srefb*depthdiffmat*middistmat).sum('depth').sum('distance')
# cc['area']=(ccarea.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
# cc['sal']=(ccvel.where(daily.salinity<34)[:,:-1,:]*daily['salinity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')/cc['trans']
cc['trans']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('depth').sum('distance')
cc['fresh']=(daily['across track velocity'][:6,:-1,:]*(daily.salinity[:6,:-1,:]-srefb)/srefb*depthdiffmat[:6,:,:]*middistmat[:6,:,:]).sum('depth').sum('distance')
cc['area']=(onesxr[:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('depth').sum('distance')
cc['sal']=(daily['across track velocity'][:6,:-1,:]*daily['salinity'][:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('distance').sum('depth')/cc['trans']
# cc['tmp']=daily.temperature.where(~isnan(ccvel)).where(daily.salinity<34).mean(dim='depth').mean(dim='distance')
# cc['den']=daily['potential density'].where(~isnan(ccvel)).where(daily.salinity<34).mean(dim='depth').mean(dim='distance')
cc['meanvel']=((daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:,:]*middistmat[:6,:,:]/1e3).sum('depth').sum('distance'))/cc['area']
# cc['meanvel'][isnan(cc['meanvel'])]=0
# cc['sal'][isnan(cc['sal'])]=nanmean(cc['sal'])


egrem={}
egrem['trans']=(egvel.where(daily.salinity<34.8)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
egrem['sal']=daily.salinity.where(~isnan(egvel)).where(daily.salinity<34.8).mean(dim='depth').mean(dim='distance')
egrem['tmp']=daily.temperature.where(~isnan(egvel)).where(daily.salinity<34.8).mean(dim='depth').mean(dim='distance')
egrem['den']=daily['potential density'].where(~isnan(egvel)).where(daily.salinity<34.8).mean(dim='depth').mean(dim='distance')
egrem['fresh']=(egvel*(daily.where(daily.salinity<34.8)['salinity'][:,:-1,:]-srefb)/srefb*depthdiffmat*middistmat).sum('distance').sum('depth')


egic={}
egic['trans from5']=(egvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
egic['trans']=(egvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
egic['fresh']=(egvel*(daily['salinity'][:,:-1,:]-srefb)/srefb*depthdiffmat*middistmat).sum('distance').sum('depth')
egic['fresh 35']=(egvel*(daily['salinity'][:,:-1,:]-srefc)/srefc*depthdiffmat*middistmat).sum('distance').sum('depth')
egic['area']=(egarea[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
egic['meanvel']=((egvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance'))/egic['area']
egic['sal']=(egvel*daily['salinity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')/egic['trans']
egic['tmp']=((daily.temperature.where(~isnan(egvel))[:,:-1,:]*depthdiffmat*middistmat/1e3).sum(dim='depth').sum(dim='distance'))/egic['area']
egic['den']=((daily['potential density'].where(~isnan(egvel))[:,:-1,:]*depthdiffmat*middistmat/1e3).sum(dim='depth').sum(dim='distance'))/egic['area']
egic['meanvel'][isnan(egic['meanvel'])]=0


def run_ave(vec,rundiv):
    runave=vec.copy()
    for ii in range(int(rundiv/2)):
        runave[ii]=nanmean(vec[:int(ii+rundiv/2)])
    for ii in range(int(rundiv/2),len(vec)-int(rundiv/2)):
        runave[ii]=nanmean(vec[int(ii-rundiv/2):int(ii+rundiv/2)])
    for ii in range(len(vec)-int(rundiv/2),len(vec)):
        runave[ii]=nanmean(vec[int(ii-rundiv/2):])
    # runave=array(runave)
    return runave

import scipy.signal as sig
# First, design the Buterworth filter
N  = 2    # Filter order
Wn = 0.02 # Cutoff frequency
B, A = sig.butter(N, Wn, output='ba')
# Second, apply the filter

def psf(field,colo,ylim1,ylim2):
    figure(figsize=(16,6))
    axvspan(datetime.datetime(2014,m1,1),datetime.datetime(2015,m2,1),color='grey',alpha=0.5)
    axvspan(datetime.datetime(2015,m1,1),datetime.datetime(2016,m2,1),color='grey',alpha=0.5)
    field.plot(alpha=0.5,color=colo)
    plot(daily.date,signal.filtfilt(B,A, field),linewidth=2,color=colo)
    ylim([ylim1,ylim2])
    title(tit)


psf(cc['fresh'],ccol,0,50)

psf(cc['trans'],ccol,-1.5,0)
savefig('../figures/xport/cctrans_movingmin.png')


psf(cc['trans to5 nosal'],ccol,-2,0)
savefig('../figures/xport/cctrans_to5_nosal.png')

psf(cc['trans to5'],ccol,-2,0)
savefig('../figures/xport/cctrans_to5_wsal.png')
egicol='purple'

psf(egic['trans from5'],egicol,-30,-5)

psf(egic['trans'],egicol,-30,-5)

#################################################################################
###################### Decomposition #################################
#################################################################################

def decompose(cur,srefchoose):
    cur['abar']=mean(cur['area'])
    cur['vbar']=mean(cur['meanvel'])
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']
    cur['sbar']=mean((cur['sal']-srefchoose)/srefchoose)
    cur['sprime']=((cur['sal']-srefchoose)/srefchoose-cur['sbar'])

    return cur

egic=decompose(egic)
cc=decompose(cc)

def plot_decompfresh(cur,colcol,y1,y2,tit):
    figure(figsize=(18,10))
    subplot(321)
    plot(cur['vprime'].date,cur['sbar']*cur['abar']*cur['vprime']*1e3,color=colcol)
    title('$\overline{s} \overline{a} v\'$: '+str(int(corrcoef(cur['sbar']*cur['abar']*cur['vprime'],cur['fresh'])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    ylabel('[mSv]')
    subplot(322)
    plot(cur['vprime'].date,cur['sbar']*cur['aprime']*cur['vbar']*1e3,color=colcol)
    title('$\overline{s} a\' \overline{v} $: '+str(int(corrcoef(cur['sbar']*cur['aprime']*cur['vbar'],cur['fresh'])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    subplot(323)
    plot(cur['vprime'].date,cur['sprime']*cur['abar']*cur['vbar']*1e3,color=colcol)
    title('$s\' \overline{a} \overline{v} $: '+str(int(corrcoef(cur['sprime']*cur['abar']*cur['vbar'],cur['fresh'])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    ylabel('[mSv]')
    subplot(324)
    plot(cur['vprime'].date,cur['sprime']*cur['abar']*cur['vprime']*1e3,color=colcol)
    title('$s\' \overline{a} v\'} $: '+str(int(corrcoef(cur['sprime']*cur['abar']*cur['vprime'],cur['fresh'])[1,0]*100)/100))
    ylim([y1,y2])
    gca().set_xticklabels('')
    subplot(325)
    remterms=cur['sprime']*cur['aprime']*cur['vbar']+cur['sprime']*cur['aprime']*cur['vprime']+cur['sbar']*cur['aprime']*cur['vprime']
    plot(cur['vprime'].date,(remterms)*1e3,color=colcol)
    ylim([y1,y2])
    title('remaining terms: '+str(int(corrcoef(remterms,cur['fresh'])[1,0]*100)/100))
    ylabel('[mSv]')
    subplot(326)
    plot(cur['vprime'].date,cur['fresh']-cur['abar']*cur['sbar']*cur['vbar']*1e3,color=colcol)
    title('Total variability')
    ylim([y1,y2])
    savefig('../figures/newtrans/Decomp_fresh_'+tit+'.png')


plot_decompfresh(egic,egicol,-100,150,'egic')

plot_decompfresh(cc,ccol,-50,100,'cc')

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

plot_decomptrans(cc,-1,1,ccol,'cc')

plot_decomptrans(egic,-15,15,egicol,'egic')

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
