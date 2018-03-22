#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


newgrid=pickle.load(open('../pickles/CF_xarray_gridplot_notid_newtheta.pickle','rb'))

# mask the fields based on bathymetry file (this version does not have extra fields on either side, thats just for plotting)


bathf=interpolate.interp1d(bathdist,bathbath)
bathonmygrid=bathf(newgrid.distance)

figure()
plot(bathdist,bathbath,'-')
plot(newgrid.distance,bathonmygrid,'o')

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

fulltrans=(daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
fulltrans_south=(daily.where(daily['across track velocity']<0)['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')


fulltrans.plot(figsize=(14,3),label='sum all')
fulltrans_south.plot(label='only sum southward transport')
legend()
savefig('../figures/newtrans/fulltrans_zerotest.png')



figure(figsize=(12,3))
plot(fulltrans.date,fulltrans,color='red')
axhline(mean(fulltrans),color='red')
# fulltrans.resample('M',how='mean',dim='date').plot(linewidth=2,color='red')
ylabel('Sv')
title('Transport across the full array')
text(datetime.datetime(2016,3,1),-40,'-21.5 $\pm$ 4.8 Sv',fontsize=18,color='red')
savefig('../figures/newtrans/fulltrans.png')
savefig('../figures/newtrans/fulltrans.pdf')

# Then, establish whether a min can be found to differentiate EGCC

minpos=zeros(len(daily.date))
minind=zeros(len(daily.date))
for tt,na in enumerate(daily.date):
    minpos[tt]=float(daily.distance[2:15][max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt]])
    minind[tt]=int(where(max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt])[0][0]+2)

minind=[int(mm) for mm in minind]

figure(figsize=(12,3))
plot(minind)


figure(figsize=(12,3))
plot(minpos)
axhline(mean(minpos))
savefig('../figures/newtrans/minvel_btwn_cf1_cf5_new.png')

#################################################################################
### Try using minpos as the barrier
#################################################################################

ccvel=daily['across track velocity'].copy()

for tt,mm in enumerate(minind):
        ccvel[mm:,:,tt]=NaN


cctrans=(ccvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cctrans_sal=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')


cctrans.plot(figsize=(12,3),label='Full CF1 water column')
axhline(0)
cctrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
cctrans_sal.plot(label='Fresher than 34 at CF1')
cctrans_sal.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
legend()
ylabel('Transport (Sv)')
title('Transport at CF1 (EGCC)')
savefig('../figures/newtrans/CF1trans_minpos.png')
min(cctrans_sal.resample('M',dim='date'))
max(cctrans_sal.resample('M',dim='date'))
mean(cctrans)

newcctrans=cctrans_sal.copy()

#################################################################################
### Coastal current transport
#################################################################################

# Use cf2 position as the division point
curdiv=5

cf1vel=daily['across track velocity'][:curdiv,:-1,:]

cctrans=(cf1vel*depthdiffmat[:curdiv,:,:]*middistmat[:curdiv,:,:]/1e3).sum('depth').sum('distance')
cctrans_sal=(daily.where(daily.salinity<34)['across track velocity'][:curdiv,:-1,:]*depthdiffmat[:curdiv,:,:]*middistmat[:curdiv,:,:]/1e3).sum('depth').sum('distance')


cctrans.plot(figsize=(12,3),label='Full CF1 water column')
axhline(0)
cctrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
cctrans_sal.plot(label='Fresher than 34 at CF1')
legend()
ylabel('Transport (Sv)')
title('Transport at CF1 (EGCC)')
savefig('../figures/newtrans/CF1trans.png')

cctrans.plot(figsize=(20,5),label='Fixed integration point')
newcctrans.plot(label='Integrate to minimum in abs vel')
legend()

newcctrans[0]

newcctrans[-1]

mean(cctrans)
cctrans.resample('W',how='mean',dim='date').plot(figsize=(12,3))


EGtottrans=(daily['across track velocity'][curdiv:,:-1,:]*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]/1e3).sum('distance').sum('depth')
EGtottrans_vel=(daily.where(daily['across track velocity']<0)['across track velocity'][curdiv:,:-1,:]*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]/1e3).sum('distance').sum('depth')

EGtottrans.plot(figsize=(12,3),label='Full water columns')
# axhline(0)
EGtottrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
EGtottrans_vel.plot(label='Only negative velocities')
fulltrans.plot(label='full array transport')
ylabel('Transport (Sv)')
legend()
title('Transport at CF2-M1 (EGC system)')
savefig('../figures/newtrans/CF2-8trans.png')

egtrans=(daily.where(daily.salinity<34.85)['across track velocity'][curdiv:,:-1,:]*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]/1e3).sum('distance').sum('depth')

newegvel=daily['across track velocity'].copy()
for tt,mm in enumerate(minind):
        newegvel[:mm,:,tt]=NaN
newegtrans=(newegvel.where(daily.salinity<34.85)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

egtrans.plot(figsize=(15,5),label='Fixed boundary')
newegtrans.plot(label='Moving boundary')


ictrans=(daily.where(daily.salinity>=34.85)['across track velocity'][curdiv:,:-1,:]*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]/1e3).sum('distance').sum('depth')
cctrans.plot(figsize=(12,3),label='East Greenland COASTAL Current')
egtrans.plot(label='East Greenlandic Current Waters')
# axhline(0)
# egtrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
ictrans.plot(label='Irminger Current')
ylabel('Transport (Sv)')
legend()
title('EGC system transports')
savefig('../figures/newtrans/EGsystem_trans.png')


egtrans.plot(figsize=(12,3),label='East Greenlandic Current Waters')
egtrans.resample('M',how='mean',dim='date').plot(linewidth=2)
ylabel('Transport (Sv)')
title('East Greenlandic Current transport')
savefig('../figures/newtrans/EGC_trans.png')




#################################################################################
###################### Freshwater transport #####################################
#################################################################################


srefa=34
srefb=34.8

ccfresh=(cf1vel*(daily.salinity[:curdiv,:-1,:]-srefa)/srefa*depthdiffmat[:curdiv,:,:]*middistmat[:curdiv,:,:]).sum('depth').sum('distance')
ccfresh_refb=(cf1vel*(daily.salinity[:curdiv,:-1,:]-srefb)/srefb*depthdiffmat[:curdiv,:,:]*middistmat[:curdiv,:,:]).sum('depth').sum('distance')


figure()
ccfresh.plot(figsize=(12,3),color='orange')
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='orange')
title('Freshwater transport in the EGCC referenced to 34')
ylabel('mSv')
savefig('../figures/newtrans/CC_fresh.png')

figure()
ccfresh_refb.plot(figsize=(12,3),color='orange')
ccfresh_refb.resample('M',dim='date',how='mean').plot(linewidth=2,color='orange')
title('Freshwater transport in the referenced to 35')
ylabel('mSv')
savefig('../figures/newtrans/CC_fresh_refb.png')

egfresh=(daily.where(daily.salinity<34.85)['across track velocity'][curdiv:,:-1,:]*(daily.where(daily.salinity<34.85)['salinity'][curdiv:,:-1,:]-srefb)/srefb*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]).sum('distance').sum('depth')

figure()
egfresh.plot(figsize=(12,3))
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport in the EGC')
ylabel('mSv')
savefig('../figures/newtrans/EGC_fresh.png')



fig, ax1 = plt.subplots(figsize=(12,3),)
ccfresh.plot(alpha=0.5,ax=ax1,color='red')
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='red',ax=ax1,label='Coastal Current')
egfresh.plot(alpha=0.5,ax=ax1)
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b',ax=ax1,label='East Greenland Current')
ax1.set_ylabel('Freshwater transport [mSv]')
ax1.set_ylim([0,100])
legend()
savefig('../figures/newtrans/EGCboth_fresh_sameaxis.png')
fig, ax1 = plt.subplots(figsize=(12,3),)
egtrans.plot(alpha=0.5,ax=ax1)
egtrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='b',ax=ax1)
ax1.set_ylabel('East Greenland Current',color='blue')
ax1.tick_params('y', colors='b')
ax1.set_ylim([-8,0])
ax2 = ax1.twinx()
cctrans.plot(alpha=0.5,ax=ax2,color='red')
cctrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='red',ax=ax2)
ax2.set_title('Transports in the EGC system [Sv]')
ax2.set_ylabel('Coastal Current',color='red')
ax2.tick_params('y', colors='r')
ax2.set_ylim([-2.5,0])
savefig('../figures/newtrans/EGCboth_trans.png')





fig, ax1 = plt.subplots(figsize=(12,3),)
egfresh.plot(alpha=0.5,ax=ax1)
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b',ax=ax1)
ax1.set_ylabel('East Greenland Current',color='blue')
ax1.tick_params('y', colors='b')
ax1.set_ylim([0,100])
ax2 = ax1.twinx()
ccfresh.plot(alpha=0.5,ax=ax2,color='red')
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='red',ax=ax2)
ax2.set_title('Freshwater transport in the EGC system [mSv]')
ax2.set_ylabel('Coastal Current',color='red')
ax2.tick_params('y', colors='r')
ax2.set_ylim([0,60])
savefig('../figures/newtrans/EGCboth_fresh.png')

#########################################################################
## Look at rel between transports and wind
#########################################################################
winddat=pickle.load(open('../pickles/wind/NARR_linet_newtheta.pickle','rb'))

fig, ax1 = plt.subplots(figsize=(12,3),)
winddat['across track wind speed'][:,2].resample('D',dim='date',how='mean').plot(ax=ax1,color='green')
winddat['across track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1,color='green')
# winddat['along track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1)
ax1.tick_params('y', colors='green')
ax1.axhline(-10,color='green')
ax2 = ax1.twinx()
ax2.tick_params('y', colors='r')
ax2.set_ylabel('EGCC transport [Sv]')
ax1.set_title('Across wind - Coastal current transport')
cctrans.plot(alpha=0.5,color='red',ax=ax2)
cctrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='red',ax=ax2)
savefig('../figures/newtrans/cctrans_acrosswind.png')

fig, ax1 = plt.subplots(figsize=(12,3),)
# winddat['across track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1)
winddat['along track wind speed'][:,2].resample('3D',dim='date',how='mean').plot(ax=ax1,color='purple',alpha=0.5)
winddat['along track wind speed'][:,2].resample('M',dim='date',how='mean').plot(ax=ax1,color='purple')
ax1.set_ylim([-15,0])
ax2 = ax1.twinx()
fulltrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2,alpha=0.5)
fulltrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2)
ax2.set_ylim([-30,-15])
ax1.set_title('Full transport and along track wind speed')
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
egtrans.resample('3D',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2,alpha=0.5)
egtrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='blue',ax=ax2)
# ax1.set_title('')
ax1.set_title('EGC transport and along track wind speed')
ax2.set_ylabel('EGC transport [Sv]')
savefig('../figures/newtrans/egc_alongwind.png')

#########################################################################
## Freshwater trans for full/rest of array
#########################################################################
ax2_srefc=35

icfresh=(daily.where(daily.salinity>=34.85)['across track velocity'][curdiv:,:-1,:]*(daily.where(daily.salinity>=34.8)['salinity'][curdiv:,:-1,:]-srefc)/srefc*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]/1e3).sum('distance').sum('depth')
icfresh.plot(figsize=(12,3))
icfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport in the IC (ref to 35)')
ylabel('mSv')
savefig('../figures/newtrans/icfresh.png')

egicfresh=(daily['across track velocity'][curdiv:,:-1,:]*(daily['salinity'][curdiv:,:-1,:]-srefc)/srefc*depthdiffmat[curdiv:,:,:]*middistmat[curdiv:,:,:]/1e3).sum('distance').sum('depth')
egicfresh.plot(figsize=(12,3))
egicfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport for full EGC/IC current (ref to 35)')
ylabel('mSv')
savefig('../figures/newtrans/egicfresh.png')

fullfresh_34=(daily['across track velocity'][:,:-1,:]*(daily['salinity'][:,:-1,:]-srefa)/srefa*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')
fullfresh_34.plot(figsize=(12,3))
fullfresh_34.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport for full array (ref to 34)')
ylabel('mSv')
savefig('../figures/newtrans/fullfresh_34.png')

fullfresh_35=(daily['across track velocity'][:,:-1,:]*(daily['salinity'][:,:-1,:]-srefc)/srefc*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')
fullfresh_35.plot(figsize=(12,3))
fullfresh_35.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport for full array (ref to 35)')
ylabel('mSv')
savefig('../figures/newtrans/fullfresh_35.png')
srefmean=mean(daily['salinity'][:,:-1,:])

fullfresh_mean=(daily['across track velocity'][:,:-1,:]*(daily['salinity'][:,:-1,:]-srefmean)/srefmean*depthdiffmat*middistmat/1e3).sum('distance').sum('depth')

fullfresh_mean.plot(figsize=(12,3))
fullfresh_mean.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
axhline(mean(fullfresh_mean))
title('Freshwater transport for full array (ref to mean salinity)')
ylabel('mSv')
savefig('../figures/newtrans/fullfresh_mean.png')

#########################################################################
## Get transport time series for top 50% of vels east of shelf
#########################################################################
