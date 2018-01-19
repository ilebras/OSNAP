#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT ###################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


daily=pickle.load(open('../pickles/xarray/CF_xarray_notid_SAtheta.pickle','rb'))

#################################################################################
# Have a quick look at CF1 evolution in time
#################################################################################

def plotmoortime(moornum):
    figure(figsize=(12,3))
    ax=contourf(daily.date.data,daily.depth,daily['across track velocity'][moornum-1,:,:],cmap=cm.RdBu_r,vmin=-1.25,vmax=1.25)
    colorbar(ticks=[-1.5,-1,-0.5,0,0.5])
    contour(daily.date.data,daily.depth,daily['across track velocity'][moornum-1,:,:],[-0.75],colors='k')
    ylim([170,0])
    ylabel('depth (m)')
    xlabel('date')
    title('CF'+str(moornum)+' across track velocity')
    savefig('../figures/hovmueller/cf'+str(moornum)+'_vel.png',bbox_inches='tight')
    savefig('../figures/hovmueller/cf'+str(moornum)+'_vel.pdf',bbox_inches='tight')

def plotmoortime(moornum):
    figure(figsize=(12,3))
    ax=contourf(daily.date.data,daily.depth,daily['across track velocity'][moornum-1,:,:],cmap=cm.RdBu_r,vmin=-1.25,vmax=1.25)
    colorbar(ticks=[-1.5,-1,-0.5,0,0.5])
    contour(daily.date.data,daily.depth,daily['across track velocity'][moornum-1,:,:],[-0.75],colors='k')
    ylim([170,0])
    ylabel('depth (m)')
    xlabel('date')
    title('CF'+str(moornum)+' across track velocity')
    savefig('../figures/hovmueller/cf'+str(moornum)+'_vel.png',bbox_inches='tight')
    savefig('../figures/hovmueller/cf'+str(moornum)+'_vel.pdf',bbox_inches='tight')

for rr in range(1,9):
    plotmoortime(rr)

#################################################################################
#################################################################################
############# Get EGCC and EGC transports ####################################
#################################################################################
#################################################################################
#
# #################################################################################
# # Quick code for looking at monthly averages
# #################################################################################
#
# def monthplot(afield):
#     figure()
#     afield.resample('M',dim='date',how='mean')[:12,:,:].plot(x='distance', y='depth', col='date', col_wrap=4)
#
# monthplot(daily['across track velocity'])
# ylim([1000,0])


#################################################################################
################ Find and examine isohalines ###################################
#################################################################################

def find_isohaline(which):

    maxdepth=pd.DataFrame(index=daily.date, columns=daily.distance)

    for j, m in enumerate(daily.distance):
        for i, d in enumerate(daily.date):
            thissal=daily.salinity[j,:,i]
            nanind=~isnan(thissal)
            if sum(nanind)==0:
                maxdepth.iloc[i,j]=nan
            elif sum((thissal[nanind]>which))==0:
                maxdepth.iloc[i,j]=max(daily.depth[nanind])
            else:
                maxdepth.iloc[i,j]=float(daily.depth[nanind][(thissal[nanind]>which)].min())

    maxdepth=maxdepth.astype('float')
    return maxdepth


max34depth=find_isohaline(34)
max348depth=find_isohaline(34.8)

colors=pal.cubehelix.perceptual_rainbow_16.get_mpl_colormap()

fig, ax = plt.subplots(1)
fig.set_size_inches(12,4)
max34depth.plot(ax=ax, cmap=colors, alpha=0.5,label=False)
g=max34depth.resample('M',closed='right').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2)
legend(loc=(1.05,0))
gca().invert_yaxis()
title('Depth of 34 isohaline along CF array')
savefig('../figures/isohalines/34tseries.png')

fig, ax = plt.subplots(1)
fig.set_size_inches(12,4)
max348depth.plot(ax=ax, cmap=colors, alpha=0.5,label=False)
num=max348depth.resample('M').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2)
num.legend(loc=(1.05,0))
gca().invert_yaxis()
title('Depth of 34.8 isohaline along CF array')
savefig('../figures/isohalines/348tseries.png')

fig, ax = plt.subplots(1)
fig.set_size_inches(12,4)
num=max34depth.resample('M').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2,linestyle='--')
max348depth.resample('M').mean().plot(ax=ax, cmap=colors, alpha=1, lw=2)
num.legend(loc=(1.05,0))
title('Depths of 34 and 34.8 isohalines along CF array')
gca().invert_yaxis()
savefig('../figures/isohalines/34and348tseries.png')

#################################################################################
###         Look at velocity magnitudes at different moorings
#################################################################################

figure(figsize=(14,3))
for rr in range(3):
    plot(daily.date,daily['across track velocity'].min(dim='depth')[rr],alpha=0.5,label='CF'+str(rr+1))
    plot(daily.resample('M',dim='date',how='mean').date,daily['across track velocity'].resample('M',dim='date',how='mean').min(dim='depth')[rr])
legend(loc=(1.05,0))
plot(daily.date,0.15*daily['across track velocity'].min(dim='depth')[0],'k')
savefig('../figures/minvels/CF1-2.png')


figure(figsize=(14,3))
for rr in range(1,3):
    plot(daily.date,daily['across track velocity'].min(dim='depth')[rr],alpha=0.75,label='CF'+str(rr+1))
    # plot(daily.resample('M',dim='date',how='mean').date,daily['across track velocity'].resample('M',dim='date',how='mean').min(dim='depth')[rr])
legend(loc=(1.05,0.2))
title('CF2 and 3 track each other closely')
savefig('../figures/minvels/CF2-3.png')


for rr in range(8):
    figure(figsize=(14,3))
    # plot(daily.date,daily['across track velocity'].min(dim='depth')[rr],alpha=0.5,label='CF'+str(rr+1))
    plot(daily.resample('M',dim='date',how='mean').date,daily['across track velocity'].resample('M',dim='date',how='mean').min(dim='depth')[rr],label='min vel')
    title('CF'+str(rr+1))
    plot(daily.resample('M',dim='date',how='mean').date,daily['across track velocity'].resample('M',dim='date',how='mean')[rr,0,:],label='surface vel')
    legend(loc=(1.05,0.2))
    ylabel('velocity (m/s)')
    savefig('../figures/velstats/CF'+str(rr+1)+'_minvelcomp_monthly.png',bbox_inches='tight')


for rr in range(8):
    figure(figsize=(14,3))
    plot(daily.date,daily['across track velocity'].min(dim='depth')[rr],label='min vel')
    axhline(0)
    title('CF'+str(rr+1))
    plot(daily.date,daily['across track velocity'][rr,0,:],label='surface vel')
    legend(loc=(1.05,0.2))
    ylabel('velocity (m/s)')
    savefig('../figures/velstats/CF'+str(rr+1)+'_minvelcomp_daily.png',bbox_inches='tight')

daily.dims

figure(figsize=(14,3))
for rr in range(8):
    plot(daily.resample('M',dim='date',how='mean').date,daily['across track velocity'].resample('M',dim='date',how='mean')[rr,0,:],label='CF'+str(rr+1))
legend(loc=(1.05,0.2))
savefig('../figures/velstats/Monthlyave_surf_all.png')



#################################################################################
# Transport -- define as solely at CF1 for now
#################################################################################
mid_dist=hstack((12,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,17))
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

shape(middistmat[:,:,:])

cf1vel=daily['across track velocity'][0,:-1,:]

cctrans=(cf1vel*depthdiffmat[0,:,:]*middistmat[0,:,:]/1e3).sum('depth')
cctrans_sal=(daily.where(daily.salinity<34)['across track velocity'][0,:-1,:]*depthdiffmat[0,:,:]*middistmat[0,:,:]/1e3).sum('depth')


cctrans.plot(figsize=(12,3),label='Full CF1 water column')
axhline(0)
cctrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
cctrans_sal.plot(label='Fresher than 34 at CF1')
legend()
ylabel('Transport (Sv)')
title('Transport at CF1 (EGCC)')
savefig('../figures/trans/CF1trans.png')

cctrans_scaled=cctrans*3

cctrans.plot(figsize=(12,3),label='')
axhline(0)
cctrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
# cctrans_sal.plot(label='Fresher than 34 at CF1')
legend()
ylabel('[Sv]')
title('EG Coastal Current transport')
savefig('../figures/trans/EGCC_trans.pdf')


cctrans.resample('W',how='mean',dim='date').plot(figsize=(12,3))


EGtottrans=(daily['across track velocity'][1:,:-1,:]*depthdiffmat[1:,:,:]*middistmat[1:,:,:]/1e3).sum('distance').sum('depth')
EGtottrans_vel=(daily.where(daily['across track velocity']<0)['across track velocity'][1:,:-1,:]*depthdiffmat[1:,:,:]*middistmat[1:,:,:]/1e3).sum('distance').sum('depth')

EGtottrans.plot(figsize=(12,3),label='Full water columns')
# axhline(0)
EGtottrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
EGtottrans_vel.plot(label='Only negative velocities')
ylabel('Transport (Sv)')
legend()
title('Transport at CF2-M1 (EGC system)')
savefig('../figures/trans/CF2-8trans.png')

egtrans=(daily.where(daily.salinity<34.8)['across track velocity'][1:,:-1,:]*depthdiffmat[1:,:,:]*middistmat[1:,:,:]/1e3).sum('distance').sum('depth')


ictrans=(daily.where(daily.salinity>=34.85)['across track velocity'][1:,:-1,:]*depthdiffmat[1:,:,:]*middistmat[1:,:,:]/1e3).sum('distance').sum('depth')
cctrans.plot(figsize=(12,3),label='East Greenland COASTAL Current')
egtrans.plot(label='East Greenlandic Current Waters')
# axhline(0)
# egtrans.resample('M',how='mean',dim='date').plot(linewidth=2,label='',)
ictrans.plot(label='Irminger Current')
ylabel('Transport (Sv)')
legend()
title('EGC system transports')
savefig('../figures/trans/EGsystem_trans.png')


egtrans.plot(figsize=(12,3),label='East Greenlandic Current Waters')
axhline(0)
egtrans.resample('M',how='mean',dim='date').plot(linewidth=2)
ylabel('[Sv]')
title('East Greenlandic Current transport')
savefig('../figures/trans/EGC_trans.png')
savefig('../figures/trans/EGC_trans.pdf')

figure()
egtrans.plot(figsize=(12,3),alpha=0.5,label='')
egtrans.resample('M',dim='date',how='mean').plot(linewidth=2,color='b',label='East Greenland Current')
cctrans_scaled.plot(alpha=0.5,label='')
cctrans_scaled.resample('M',dim='date',how='mean').plot(linewidth=2,color='orange',label='Coastal Current (x 3)')
title('Transport in the EGC system')
ylabel('[Sv]')
legend()
savefig('../figures/trans/EGCboth_trans.png')
savefig('../figures/trans/EGCboth_trans.pdf',bbox_inches='tight')

ictrans.plot(figsize=(12,3))
ictrans.resample('M',how='mean',dim='date').plot(linewidth=2)
ylabel('Transport (Sv)')
title('Irminger Current transport')
savefig('../figures/trans/IC_trans.png')


hexbin(daily.salinity.values.flatten(),daily.temperature.values.flatten(),bins='log',cmap=cm.hot_r)
axvline(34.8,color='k')
colorbar(label='[log of number of measurements]')

ylabel('potential temperature [$^\circ$ C]')
xlabel('salinity')
title('Separation of Polar and Atlantic Water at 34.8')
savefig('../figures/trans/TS_separation.png')
savefig('../figures/trans/TS_separation.pdf',bbox_inches='tight')


#################################################################################
###################### Freshwater transport #####################################
#################################################################################


srefa=34
srefb=34.8

ccfresh=(cf1vel*(daily.salinity[0,:-1,:]-srefa)/srefa*depthdiffmat[0,:,:]*middistmat[0,:,:]).sum('depth')
ccfresh_refb=(cf1vel*(daily.salinity[0,:-1,:]-srefb)/srefb*depthdiffmat[0,:,:]*middistmat[0,:,:]).sum('depth')

ccfresh_scaled=ccfresh*2

figure()
ccfresh.plot(figsize=(12,3),color='orange')
ccfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='orange')
title('Freshwater transport in the EGCC referenced to 34')
ylabel('mSv')
savefig('../figures/trans/CC_fresh.png')

figure()
ccfresh_refb.plot(figsize=(12,3),color='orange')
ccfresh_refb.resample('M',dim='date',how='mean').plot(linewidth=2,color='orange')
title('Freshwater transport in the referenced to 35')
ylabel('mSv')
savefig('../figures/trans/CC_fresh_refb.png')

egfresh=(daily.where(daily.salinity<34.85)['across track velocity'][1:,:-1,:]*(daily.where(daily.salinity<34.85)['salinity'][1:,:-1,:]-srefb)/srefb*depthdiffmat[1:,:,:]*middistmat[1:,:,:]).sum('distance').sum('depth')

figure()
egfresh.plot(figsize=(12,3))
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport in the EGC')
ylabel('mSv')
savefig('../figures/trans/EGC_fresh.png')

figure()
egfresh.plot(figsize=(12,3),alpha=0.5)
egfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b',label='East Greenland Current')
ccfresh_scaled.plot(alpha=0.5)
ccfresh_scaled.resample('M',dim='date',how='mean').plot(linewidth=2,color='orange',label='Coastal Current (x 2)')
title('Freshwater transport in the EGC system')
ylabel('mSv')
legend()
savefig('../figures/trans/EGCboth_fresh.png')
savefig('../figures/trans/EGCboth_fresh.pdf',bbox_inches='tight')

icfresh=(daily.where(daily.salinity>=34.85)['across track velocity'][1:,:-1,:]*(daily.where(daily.salinity>=34.85)['salinity'][1:,:-1,:]-srefb)/srefb*depthdiffmat[1:,:,:]*middistmat[1:,:,:]/1e3).sum('distance').sum('depth')
icfresh.plot(figsize=(12,3))
icfresh.resample('M',dim='date',how='mean').plot(linewidth=2,color='b')
title('Freshwater transport in the IC')
ylabel('mSv')
