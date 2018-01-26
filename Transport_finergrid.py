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


#################################################################################
#################################################################################
############# Get EGCC and EGC transports ####################################
#################################################################################
#################################################################################

#quick monthly plot

def monthplot(afield):
    figure()
    afield.resample('M',dim='date',how='mean').plot(x='distance', y='depth', col='date', col_wrap=4)

monthplot(daily['across track velocity'])
ylim([1000,0])


#################################################################################
################ Find and examine isohalines ###################################
#################################################################################
#
# def find_isohaline(which):
#
#     maxdepth=pd.DataFrame(index=daily.date, columns=daily.distance)
#
#     for j, m in enumerate(daily.distance):
#         for i, d in enumerate(daily.date):
#             thissal=daily.salinity[j,:,i]
#             nanind=~isnan(thissal)
#             if sum(nanind)==0:
#                 maxdepth.iloc[i,j]=nan
#             elif sum((thissal[nanind]>which))==0:
#                 maxdepth.iloc[i,j]=max(daily.depth[nanind])
#             else:
#                 maxdepth.iloc[i,j]=float(daily.depth[nanind][(thissal[nanind]>which)].min())
#
#     maxdepth=maxdepth.astype('float')
#     return maxdepth
#
#
# max34depth=find_isohaline(34)
# max348depth=find_isohaline(34.8)
#
# colors=pal.cubehelix.perceptual_rainbow_16.get_mpl_colormap()
#
# fig, ax = plt.subplots(1)
# fig.set_size_inches(12,4)
# max34depth.plot(ax=ax, cmap=colors, alpha=0.5,label=False)
# g=max34depth.resample('M',closed='right').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2)
# legend(loc=(1.05,0))
# gca().invert_yaxis()
# title('Depth of 34 isohaline along CF array')
# savefig('../figures/isohalines/34tseries.png')
#
# fig, ax = plt.subplots(1)
# fig.set_size_inches(12,4)
# max348depth.plot(ax=ax, cmap=colors, alpha=0.5,label=False)
# num=max348depth.resample('M').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2)
# num.legend(loc=(1.05,0))
# gca().invert_yaxis()
# title('Depth of 34.8 isohaline along CF array')
# savefig('../figures/isohalines/348tseries.png')
#
# fig, ax = plt.subplots(1)
# fig.set_size_inches(12,4)
# num=max34depth.resample('M').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2,linestyle='--')
# max348depth.resample('M').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2)
# num.legend(loc=(1.05,0))
# title('Depths of 34 and 34.8 isohalines along CF array')
# gca().invert_yaxis()
# savefig('../figures/isohalines/34and348tseries.png')


#################################################################################
# Transport
#################################################################################
mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))




# First, get transport across full array
#
# fulltrans=(daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
# fulltrans_south=(daily.where(daily['across track velocity']<0)['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
#
# fulltrans.plot(figsize=(14,3),label='sum all',color='k')
# fulltrans_south.plot(label='only sum southward transport',color='grey')
# legend()
# savefig('../figures/newtrans/fulltrans_zerotest.png')
#
#
#
# figure(figsize=(12,3))
# plot(fulltrans.date,fulltrans,color='k')
# axhline(mean(fulltrans),color='k')
# # fulltrans.resample('M',how='mean',dim='date').plot(linewidth=2,color='red')
# ylabel('Sv')
# title('Transport across the full array')
# text(datetime.datetime(2016,3,1),-40,'-21.5 $\pm$ 4.8 Sv',fontsize=18,color='k')
# savefig('../figures/newtrans/fulltrans.png')
# savefig('../figures/newtrans/fulltrans.pdf')

# Establish min in abs vel to differentiate EGCC
minind=zeros(len(daily.date))
for tt,na in enumerate(daily.date):
    minind[tt]=where(max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt])[0][0]+2

minind=[int(mm) for mm in minind]

mindxr=daily.date.copy()
mindxr.values=daily.distance[minind]

mindxr.resample('2W',dim='date').plot()

figure(figsize=(10,4))
plot(daily.date,daily.distance[minind],color='grey')
mindxr.resample('W',dim='date').plot(color='k')
title('distance at which minimum between currents occurs [km]')
ylabel('')
savefig('../figures/newtrans/minpos.png')

ccvel=daily['across track velocity'].copy()
ccarea=daily['across track velocity'].copy()
ccarea[:]=1
for tt,mm in enumerate(minind):
        ccvel[mm:,:,tt]=NaN
        ccarea[mm:,:,tt]=NaN

ccsal=daily['salinity'].copy()
for tt,mm in enumerate(minind):
        ccsal[mm:,:,tt]=NaN


#################################################################################
### Try using minpos as the barrier
#################################################################################

cctrans_nosal=(ccvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cctrans=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

cctrans_nosal.plot(figsize=(12,3),label='All transport west of vel min',color=ccol)
axhline(0,color='grey')
cctrans_nosal.resample('M',how='mean',dim='date').plot(linewidth=2,label='',color=ccol)
cctrans.plot(label='Also fresher than 34',color='blue')
cctrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',color='blue')
legend()
ylabel('Transport (Sv)')
title('EGCC transport')
savefig('../figures/newtrans/CC_trans_wnosal.png')


egvel=daily['across track velocity'].where(daily['potential density']<27.8)
egvel_fulld=daily['across track velocity'].copy()
egarea=daily['across track velocity'].copy()
egarea[:]=1
egarea=egarea.where(daily['potential density']<27.8)
for tt,mm in enumerate(minind):
        egvel[:mm,:,tt]=NaN
        egvel_fulld[:mm,:,tt]=NaN
        egarea[:mm,:,tt]=NaN

egtrans=(egvel.where(daily.salinity<34.85)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

#note: similar as for cctrans, doesn't make much of a difference whether moving minvel is used or not, but using it for thoroughness here.

egtrans.plot(figsize=(12,3),label='East Greenlandic Current Waters',color=egcol)
egtrans.resample('M',how='mean',dim='date').plot(linewidth=2,color=egcol)
ylabel('Transport (Sv)')
title('East Greenland Current transport')
savefig('../figures/newtrans/EGC_trans.png')




#################################################################################
###################### Freshwater transport #####################################
#################################################################################


srefa=34
srefb=34.8
srefc=35

ccfresh=(ccvel*(daily.salinity[:,:-1,:]-srefa)/srefa*depthdiffmat[:,:,:]*middistmat[:,:,:]).sum('depth').sum('distance')
ccfresh_refb=(ccvel*(daily.salinity[:,:-1,:]-srefb)/srefb*depthdiffmat[:,:,:]*middistmat[:,:,:]).sum('depth').sum('distance')
ccfresh_refc=(ccvel*(daily.salinity[:,:-1,:]-srefc)/srefc*depthdiffmat[:,:,:]*middistmat[:,:,:]).sum('depth').sum('distance')

egfresh=(egvel*(daily.where(daily.salinity<34.85)['salinity'][:,:-1,:]-srefb)/srefb*depthdiffmat*middistmat).sum('distance').sum('depth')

egicfresh=(egvel*(daily['salinity'][:,:-1,:]-srefc)/srefc*depthdiffmat*middistmat).sum('distance').sum('depth')


figure()
ccfresh.plot(figsize=(12,3),color=ccol)
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol)
title('Freshwater transport in the EGCC referenced to 34')
ylabel('mSv')
savefig('../figures/newtrans/CC_fresh.png')

figure()
ccfresh_refb.plot(figsize=(12,3),color=ccol)
ccfresh_refb.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol)
title('Freshwater transport in the EGCC referenced to 34.8')
ylabel('mSv')
savefig('../figures/newtrans/CC_fresh_refb.png')


figure()
egfresh.plot(figsize=(12,3),color=egcol)
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=egcol)
title('Freshwater transport in the EGC referenced to 34.8')
ylabel('mSv')
savefig('../figures/newtrans/EGC_fresh.png')



fig, ax1 = plt.subplots(figsize=(12,3),)
ccfresh.plot(alpha=0.5,ax=ax1,color=ccol)
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax1,label='Coastal Current')
egfresh.plot(alpha=0.5,ax=ax1,color=egcol)
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax1,label='East Greenland Current')
ax1.set_ylabel('Freshwater transport ref to 34.8 [mSv]')
ax1.set_ylim([0,100])
legend()
savefig('../figures/newtrans/EGCboth_fresh_sameaxis.png')

fig, ax1 = plt.subplots(figsize=(12,3),)
egtrans.plot(alpha=0.5,ax=ax1,color=egcol)
egtrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax1)
ax1.set_ylabel('East Greenland Current',color=egcol)
ax1.tick_params('y', colors=egcol)
ax1.set_ylim([-8,0])
ax2 = ax1.twinx()
cctrans.plot(alpha=0.5,ax=ax2,color=ccol)
cctrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax2)
ax2.set_title('Transports in the EGC system [Sv]')
ax2.set_ylabel('Coastal Current',color=ccol)
ax2.tick_params('y', colors=ccol)
ax2.set_ylim([-2.5,0])
savefig('../figures/newtrans/EGCboth_trans.png')




fig, ax1 = plt.subplots(figsize=(12,3),)
egfresh.plot(alpha=0.5,ax=ax1,color=egcol)
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax1)
ax1.set_ylabel('East Greenland Current',color=egcol)
ax1.tick_params('y', colors=egcol)
ax1.set_ylim([0,100])
ax2 = ax1.twinx()
ccfresh.plot(alpha=0.5,ax=ax2,color=ccol)
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax2)
ax2.set_title('Freshwater transport in the EGC system [mSv]')
ax2.set_ylabel('Coastal Current',color=ccol)
ax2.tick_params('y', colors=ccol)
ax2.set_ylim([0,60])
savefig('../figures/newtrans/EGCboth_fresh.png')


fig, ax1 = plt.subplots(figsize=(12,3),)
ccfresh_refc.plot(alpha=0.5,ax=ax1,color=ccol)
ccfresh_refc.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax1,label='Coastal Current')
egicfresh.plot(alpha=0.5,ax=ax1,color=egicol)
egicfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=egicol,ax=ax1,label='Slope Current')
ax1.set_ylabel('Freshwater transport ref to 35 [mSv]')
ax1.set_ylim([0,120])
legend()
savefig('../figures/newtrans/EGCC-EGIC_fresh_sameaxis.png')


fig, ax1 = plt.subplots(figsize=(12,3),)
egicfresh.plot(alpha=0.5,ax=ax1,color=egicol)
egicfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color=egicol,ax=ax1)
ax1.set_ylabel('East Greenland - Irminger Current',color=egicol)
ax1.tick_params('y', colors=egicol)
ax1.set_ylim([0,120])
ax2 = ax1.twinx()
ccfresh_refc.plot(alpha=0.5,ax=ax2,color=ccol)
ccfresh_refc.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax2)
ax2.set_title('Freshwater transport in the EGC system [mSv]')
ax2.set_ylabel('Coastal Current',color=ccol)
ax2.tick_params('y', colors=ccol)
ax2.set_ylim([0,60])
savefig('../figures/newtrans/EGCC-EGIC_fresh.png')

#########################################################################
## Look at rel between transports and wind
#########################################################################
winddat=pickle.load(open('../pickles/wind/NARR_linet_newtheta.pickle','rb'))

wlcol='y'
wccol='purple'

fig, ax1 = plt.subplots(figsize=(12,3),)
winddat['across track wind speed'][:,2].resample('D',dim='date',how='mean').plot(ax=ax1,color=wlcol)
winddat['across track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1,color=wlcol)
# winddat['along track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1)
ax1.tick_params('y', colors=wlcol)
ax1.axhline(-10,color=wlcol)
ax2 = ax1.twinx()
ax2.tick_params('y', colors=ccol)
ax2.set_ylabel('EGCC transport [Sv]')
ax1.set_title('Along flow wind - Coastal current transport')
cctrans.plot(alpha=0.5,color=ccol,ax=ax2)
cctrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=ccol,ax=ax2)
savefig('../figures/newtrans/cctrans_acrosswind.png')

fig, ax1 = plt.subplots(figsize=(12,3),)
# winddat['across track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1)
winddat['along track wind speed'][:,2].resample('3D',dim='date',how='mean').plot(ax=ax1,color='purple',alpha=0.5)
winddat['along track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1,color='purple')
ax1.set_ylim([-15,0])
ax2 = ax1.twinx()
fulltrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color='k',ax=ax2,alpha=0.5)
fulltrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='k',ax=ax2)
ax2.set_ylim([-30,-15])
ax1.set_title('Full transport and across flow wind speed')
ax2.set_ylabel('Full transport [Sv]')
savefig('../figures/newtrans/fulltrans_alongwind.png')

fig, ax1 = plt.subplots(figsize=(12,3),)
# winddat['across track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1)
winddat['along track wind speed'][:,2].resample('3D',dim='date',how='mean').plot(ax=ax1,color='purple',alpha=0.5)
winddat['along track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1,color='purple')
ax1.set_ylim([-15,0])
ax2 = ax1.twinx()
# fulltrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2,alpha=0.5)
# fulltrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2)
# ax2.set_ylim([-30,-15])
egtrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax2,alpha=0.5)
egtrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax2)

ax1.set_title('EGC transport and across flow wind speed')
ax2.set_ylabel('EGC transport [Sv]')
savefig('../figures/newtrans/egc_alongwind.png')

fig, ax1 = plt.subplots(figsize=(12,3),)
# winddat['across track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1)
winddat['along track wind speed'][:,2].resample('3D',dim='date',how='mean').plot(ax=ax1,color='purple',alpha=0.5)
winddat['along track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1,color='purple')
ax1.set_ylim([-15,0])
ax2 = ax1.twinx()
# fulltrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2,alpha=0.5)
# fulltrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2)
# ax2.set_ylim([-30,-15])
egictrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax2,alpha=0.5)
egictrans.resample('M',dim='date',how='mean').plot(linewidth=2,color=egcol,ax=ax2)

ax1.set_title('EGIC transport and across flow wind speed')
ax2.set_ylabel('EGIC transport [Sv]')
savefig('../figures/newtrans/egic_alongwind.png')

#########################################################################
## Freshwater trans for full/rest of array
#########################################################################

#
# fullfresh_34=(daily['across track velocity'][:,:-1,:]*(daily['salinity'][:,:-1,:]-srefa)/srefa*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')
# fullfresh_34.plot(figsize=(12,3),color='k')
# fullfresh_34.resample('M',dim='date',how='mean').plot(linewidth=2,color='k')
# title('Freshwater transport for full array (ref to 34)')
# ylabel('mSv')
# savefig('../figures/newtrans/fullfresh_34.png')
#
# fullfresh_35=(daily['across track velocity'][:,:-1,:]*(daily['salinity'][:,:-1,:]-srefc)/srefc*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')
# fullfresh_35.plot(figsize=(12,3),color='k')
# fullfresh_35.resample('M',dim='date',how='mean').plot(linewidth=2,color='k')
# title('Freshwater transport for full array (ref to 35)')
# ylabel('mSv')
# savefig('../figures/newtrans/fullfresh_35.png')
# srefmean=mean(daily['salinity'][:,:-1,:])
#
# fullfresh_mean=(daily['across track velocity'][:,:-1,:]*(daily['salinity'][:,:-1,:]-srefmean)/srefmean*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')
#
# fullfresh_mean.plot(figsize=(12,3),color='k')
# fullfresh_mean.resample('M',dim='date',how='mean').plot(linewidth=2,color='k')
# axhline(mean(fullfresh_mean),color='k')
# title('Freshwater transport for full array (ref to mean salinity)')
# ylabel('mSv')
# savefig('../figures/newtrans/fullfresh_mean.png')

#########################################################################
################# Water mass properties of each current ##################
#########################################################################

cc={}
cc['trans']=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cc['fresh']=(ccvel.where(daily.salinity<34)[:,:-1,:]*(daily.salinity[:,:-1,:]-srefc)/srefc*depthdiffmat*middistmat).sum('depth').sum('distance')
cc['area']=(ccarea.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cc['sal']=(ccvel.where(daily.salinity<34)[:,:-1,:]*daily['salinity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')/cc['trans']

cc['tmp']=daily.temperature.where(~isnan(ccvel)).where(daily.salinity<34).mean(dim='depth').mean(dim='distance')
cc['den']=daily['potential density'].where(~isnan(ccvel)).where(daily.salinity<34).mean(dim='depth').mean(dim='distance')
cc['meanvel']=((ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance'))/cc['area']
cc['meanvel'][isnan(cc['meanvel'])]=0
cc['sal'][isnan(cc['sal'])]=nanmean(cc['sal'])

cc['sal'].plot()

cc['fresh'].plot()

egrem={}
egrem['trans']=egtrans
egrem['sal']=daily.salinity.where(~isnan(egvel)).where(daily.salinity<34.85).mean(dim='depth').mean(dim='distance')
egrem['tmp']=daily.temperature.where(~isnan(egvel)).where(daily.salinity<34.85).mean(dim='depth').mean(dim='distance')
egrem['den']=daily['potential density'].where(~isnan(egvel)).where(daily.salinity<34.85).mean(dim='depth').mean(dim='distance')
egrem['fresh']=(egvel*(daily.where(daily.salinity<34.85)['salinity'][:,:-1,:]-srefb)/srefb*depthdiffmat*middistmat).sum('distance').sum('depth')


# egictrans_fulld=(egvel_fulld[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
# egictrans.plot()
# egictrans_fulld.plot()
# (egictrans_fulld-egictrans).plot()


egic={}
egic['trans']=(egvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
egic['fresh']=(egvel*(daily['salinity'][:,:-1,:]-srefc)/srefc*depthdiffmat*middistmat).sum('distance').sum('depth')
egic['area']=(egarea[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
egic['meanvel']=((egvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance'))/egic['area']
egic['sal']=(egvel*daily['salinity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')/egic['trans']
egic['tmp']=((daily.temperature.where(~isnan(egvel))[:,:-1,:]*depthdiffmat*middistmat/1e3).sum(dim='depth').sum(dim='distance'))/egic['area']
egic['den']=((daily['potential density'].where(~isnan(egvel))[:,:-1,:]*depthdiffmat*middistmat/1e3).sum(dim='depth').sum(dim='distance'))/egic['area']
egic['meanvel'][isnan(egic['meanvel'])]=0


egic['area'].plot()

def decompose(cur):
    cur['abar']=mean(cur['area'])
    cur['vbar']=mean(cur['meanvel'])
    cur['aprime']=cur['area']-cur['abar']
    cur['vprime']=cur['meanvel']-cur['vbar']
    cur['sbar']=mean((cur['sal']-srefc)/srefc)
    cur['sprime']=((cur['sal']-srefc)/srefc-cur['sbar'])

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

#########################################################################
## Decomposition of transports into components
#########################################################################
cc['trans'].plot()


#########################################################################
## Position of the currents? (velocity maxima and width of current)
#########################################################################


# egictrans_jneg=(egvel.where(egvel<0)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance') #Very little difference if I only addup the negative vels...
# Tried defining EGIC based on a velocity criterion, but this ended up not really making sense to do -- may as well just use it all.
# egvel_cond=egvel.copy()
# for dd,kk in enumerate(daily.date.values):
#     egvel_cond[:,:,dd]=egvel[:,:,dd].where(egvel[:,:,dd]<0.5*egvel[:,:,dd].min())
#
#
# egvel_cond[:,:,20].plot()
# egictrans_velcond=(egvel_cond[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
#
# egictrans_velcond.plot(color='grey',label='only vels greater than 0.5 the max each day')
# egictrans.plot(color=egicol,alpha=0.85,label='all vels off-shore of min')
# ylabel('Transport (Sv)')
# legend()
# savefig('../figures/newtrans/EGIC_transcomp050.png')


figure(figsize=(10,3))
egictrans.resample('M',dim='date').plot(color=egicol,label='EG-IC')
egictrans.plot(color=egicol,alpha=0.5)
axhline(mean(egictrans),color=egicol)
ylabel('EGIC transport (Sv)')
savefig('../figures/newtrans/EGIC_trans.png')
egtrans.resample('M',dim='date').plot(color=egcol,label='EGC')
ylabel('Transport (Sv)')
egtrans.plot(color=egcol,alpha=0.5)
axhline(mean(egtrans),color=egcol)
legend()
savefig('../figures/newtrans/EGCIC_bothtrans.png')



figure(figsize=(10,3))
(egictrans-egtrans).resample('M',dim='date').plot(color=icol,label='IC')
(egictrans-egtrans).plot(color=icol,alpha=0.5)
axhline(mean(egictrans-egtrans),color=icol)
ylabel('IC transport (Sv)')
savefig('../figures/newtrans/IC_trans.png')

egicmaxvel=egvel.min('depth').min(dim='distance')

figure(figsize=(10,3))
egicmaxvel.plot(color=egicol,alpha=0.8)
axhline(mean(egicmaxvel),color=egicol)
title('Maximum velocity in EGIC')
savefig('../figures/newtrans/EGIC_maxvel.png')

ccmaxvel=ccvel.min('depth').min(dim='distance')

figure(figsize=(10,3))
ccmaxvel.plot(color=ccol,label='coastal current',alpha=0.8)
axhline(mean(ccmaxvel),color=ccol)
title('Maximum velocity in EGCC')
savefig('../figures/newtrans/EGCC_maxvel.png')
egicmaxvel.plot(color=egicol,label='slope current',alpha=0.8)
axhline(mean(egicmaxvel),color=egicol)
title('Maximum velocities')
legend()
savefig('../figures/newtrans/EGC_maxvelcomp.png')

egind=egvel.min(dim='depth').argmin(dim='distance')
ccind=ccvel.min(dim='depth').argmin(dim='distance')


figure(figsize=(14,3))
plot(daily.date,egind,color=egicol,label='pos of max slope current vel')
plot(daily.date,minind,color='grey',label='mi vel between currents')
plot(daily.date,ccind,color=ccol,label='pos of max coastal current vel')
plot(daily.date,ccind.where(ccind>5),'bo')
plot(daily.date,egind.where(egind<13),'ro')

ccind>5
dat=daily.copy()

def onecont(field,tit,vrange,coloor,hlevs,nomoorlines=0):
    ax1=contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    colorbar()
    ax2=contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
    clabel(ax2)
    fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance (km)')
    ylabel('depth (m)')
    xlim([-5,100])
    ylim([2200,0])
    title(tit)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return ax1,ax2


#
# for cc in where(array(ccind>5))[0]:
#         figure(figsize=(12,3))
#         subplot(121)
#         onecont(dat['across track velocity'][:,:,cc].T,'Across track velocity on '+str(daily.date[cc].values)[:10],arange(-0.8,0.81,0.01),cm.RdBu_r,arange(-0.8,0.9,0.2))
#         subplot(122)
#         onecont(dat['salinity'][:,:,cc].T,'Salinity on '+str(daily.date[cc].values)[:10],univec['sal'][1],sal_cmap,univec['sal'][3])
#         savefig('../figures/curpos/cc_close/salvelsec_'+str(cc)+'.png')
# #
#
# for cc in where(array(egind<13))[0]:
#         figure(figsize=(12,3))
#         subplot(121)
#         onecont(dat['across track velocity'][:,:,cc].T,'Across track velocity on '+str(daily.date[cc].values)[:10],arange(-0.8,0.81,0.01),cm.RdBu_r,arange(-0.8,0.9,0.2))
#         subplot(122)
#         onecont(dat['salinity'][:,:,cc].T,'Salinity on '+str(daily.date[cc].values)[:10],univec['sal'][1],sal_cmap,univec['sal'][3])
#         savefig('../figures/curpos/eg_close/salvelsec_'+str(cc)+'.png')
