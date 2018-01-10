#################################################################################
#################################################################################
#################################################################################
###### CREATE PLOTS FOR PRESENTATION AT IRMINGER SEA WORKSHOP NOV 2017  #########
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

#################################################################################
#################### TRANSPORT PLOTS  ###########################################
#################################################################################


dat=pickle.load(open('../pickles/CF_xarray_gridplot_notid_newtheta.pickle','rb'))

bathf=interpolate.interp1d(bathdist,bathbath)
bathonmygrid=bathf(dat.distance)

plot(bathbath,'o')

daily=dat.copy()
for vv in dat:
    if vv[0]!='d':
        print(vv)
        for dd,adist in enumerate(daily.distance):
            daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=bathonmygrid[dd])

mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

# First, get transport across full array

fulltrans=(daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')


fulltransmean=(daily['across track velocity'].mean(dim='date')[:,:-1]*depthdiffmat[:,:,0]*middistmat[:,:,0]/1e3).sum('depth').sum('distance')

mean(fulltrans)

std(fulltrans)

mean(fulltrans)
fulltransmean

fulltrans.plot(figsize=(12,3),color='red')
axhline(mean(fulltrans),color='red')
# axhline(float(fulltransmean),linestyle='--',color='purple')
ylabel('Sv')
xlabel('')
title('Transport across the full array',fontsize=24)
text(datetime.datetime(2016,2,1),-38,'-22.2 $\pm$ 5.1 Sv',fontsize=24,color='red')

savefig('../figures/irminger17/fulltrans.pdf',bbox_inches='tight')

# sidenote: freshwater transport across full line doesn't have large eddy comp either... (~ 3 mSv) mean ~ 24 mSv
# sref=34.8
# fullfresh=(daily['across track velocity'][:,:-1,:]*(sref-daily.salinity[:,:-1,:])/sref*depthdiffmat*middistmat).sum('depth').sum('distance')
#
# fullfreshmean=(daily['across track velocity'].mean(dim='date')[:,:-1]*(sref-daily.salinity.mean(dim='date')[:,:-1])/sref*depthdiffmat[:,:,0]*middistmat[:,:,0]).sum('depth').sum('distance')
#
# fullfresh.plot()
# axhline(mean(fullfresh))
# axhline(fullfreshmean,color='r')
# mean(fullfresh)
#
# mean(fullfresh)-fullfreshmean



# Then, establish min in abs vel to differentiate EGCC
minind=zeros(len(daily.date))
for tt,na in enumerate(daily.date):
    minind[tt]=where(max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt])[0][0]+2

minind=[int(mm) for mm in minind]

ccvel=daily['across track velocity'].copy()
for tt,mm in enumerate(minind):
        ccvel[mm:,:,tt]=NaN

ccsal=daily['salinity'].copy()
for tt,mm in enumerate(minind):
        ccsal[mm:,:,tt]=NaN

figure(figsize=(12,4))
ax0=subplot(111)
ccsal.where(daily.salinity<34).mean(dim='depth').mean(dim='distance').plot(ax=ax0)
ccsal.mean(dim='depth').mean(dim='distance').plot(ax=ax0)
ccsal.where(daily.salinity<34).mean(dim='depth').mean(dim='distance').resample('M',dim='date').plot(ax=ax0,color='blue')
ccsal.mean(dim='depth').mean(dim='distance').resample('M',dim='date').plot(ax=ax0,color='red')

plot(minind)

cctrans=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

cctrans.plot(figsize=(12,3),color='blue')
cctrans.resample('M',dim='date').plot(color='blue')
ylabel('Sv',fontsize=18)
xlabel('')
title('Transport of East Greenland COASTAL Current',fontsize=24)
text(datetime.datetime(2015,5,1),-2,'-0.5 $\pm$ 0.4 Sv',fontsize=24,color='blue')
savefig('../figures/irminger17/cctrans.pdf',bbox_inches='tight')

egvel=daily['across track velocity'].where(daily.salinity<34.85)
for tt,mm in enumerate(minind):
        egvel[:mm,:,tt]=NaN
egtrans=(egvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')


mean(egtrans)

std(egtrans)

egtrans.plot(figsize=(12,3),color='orange')
egtrans.resample('M',dim='date').plot(color='orange')
ylabel('Sv',fontsize=18)
xlabel('')
title('Transport of East Greenland Current Remnant',fontsize=24)
text(datetime.datetime(2015,5,1),-8,'-2.3 $\pm$ 1.6 Sv',fontsize=18,color='orange')
savefig('../figures/irminger17/egtrans.pdf',bbox_inches='tight')


fig, ax1 = plt.subplots(figsize=(12,3),)
cctrans.plot(alpha=0.6,ax=ax1,color='blue')
cctrans.resample('M',dim='date',how='mean').plot(linewidth=3,color='b',ax=ax1)
ax1.tick_params('y')
ax1.set_ylabel('Coastal Current',color='blue',fontsize=18)
ax1.set_yticks(arange(-1.5,0.1,0.5))
ax1.set_ylim([-1.5,0])
ax2 = ax1.twinx()
egtrans.plot(alpha=0.6,ax=ax2,color='orange')
ax2.plot(datetime.datetime(2014,8,1),0,'b',linewidth=3,label='Coastal Current')
egtrans.resample('M',dim='date',how='mean').plot(linewidth=3,color='orange',ax=ax2,label='EGC remnant')
ax2.set_title('Transports in the EGC system [x 10$^6$ m$^3$/s]',fontsize=24)
ax2.set_ylabel('EGC remnant',color='orange',fontsize=20)
ax2.tick_params('y')#, colors='orange')
ax2.set_ylim([-8,0])
ax1.set_xlabel('')
ax1.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,9,1)])
# legend()
savefig('../figures/irminger17/EGCboth_trans_forposter.pdf',bbox_inches='tight')
# savefig('../figures/irminger17/EGCboth_trans.pdf',bbox_inches='tight')



sref=34.8

ccfresh=(ccvel[:,:-1,:]*(sref-daily.salinity[:,:-1,:])/sref*depthdiffmat*middistmat).sum('depth').sum('distance')
egfresh=(egvel[:,:-1,:]*(sref-daily.salinity[:,:-1,:])/sref*depthdiffmat*middistmat).sum('depth').sum('distance')

ccfresh.plot(figsize=(12,3),color='blue',alpha=0.6)
ccfresh.resample('M',dim='date').plot(color='blue',label='Coastal Current',linewidth=3)
egfresh.plot(color='orange',alpha=0.6)
egfresh.resample('M',dim='date').plot(color='orange',label='EGC remnant',linewidth=3)
xlabel('')
ylim([-120,10])
# legend()
text(datetime.datetime(2014,8,1),-100,'Coastal Current',color='blue',fontsize=18)
text(datetime.datetime(2014,8,1),-115,'EGC remnant',color='orange',fontsize=18)
title('Freshwater transports referenced to 34.8 [mSv]',fontsize=24)
savefig('../figures/irminger17/freshtrans.pdf',bbox_inches='tight')



#################################################################################
############################## TS plot  ########################################
#################################################################################

hexbin(daily.salinity.values.flatten(),daily.temperature.values.flatten(),bins='log',cmap=cm.hot_r)
axvline(34.8,color='k')
axvline(34,color='k')

fslab=16
text(32.5,-2,'EGCC',color='blue',fontsize=fslab)
text(34.1,-2,'EGC',color='orange',fontsize=fslab)
text(34.85,-2,'AW',color='black',fontsize=fslab)
ylim([-3,10])
colorbar(label='[log # measurements]')

ylabel('potential temperature  [$^\circ$ C]')
xlabel('salinity')
savefig('../figures/irminger17/TS_separation.pdf',bbox_inches='tight')


#################################################################################
############################## SECTIONS  ########################################
#################################################################################

def mkbox(axc,x1,x2,y1,y2):
    axc.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],color='purple',linewidth=3)

# needing to reload because its blanking out, not sure why...
dat=pickle.load(open('../pickles/CF_xarray_gridplot_notid_pdenfix.pickle','rb'))

def overcont(field1,field2,tit,vrange,coloor,hlevs,nomoorlines=0):
    ax1=contourf(dat.distance,dat.depth,field1,vrange,cmap=coloor)
    ax2=contour(dat.distance,dat.depth,field2,levels=hlevs,colors='k')
    fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance (km)')
    ylabel('depth (m)')
    xlim([-5,100])
    ylim([2200,0])
    title(tit,fontsize=18)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return ax1


univec['pden']=['potential density',linspace(26,28,41),cm.BuPu,arange(26,28.1,0.2),'[kg/m$^3$]']


def meanplot2(field1,field2,noylab=0,hq=0):
    if hq==0:
        hlevs=univec[field2][3]
    else:
        hlevs=hq
    ax1=overcont(dat[univec[field1][0]].mean(dim='date').T,dat[univec[field2][0]].mean(dim='date').T,'Mean '+univec[field1][0],univec[field1][1],univec[field1][2],hlevs)
    colorbar(ax1,label=univec[field1][4],ticks=univec[field1][3][::2])
    if noylab==1:
        gca().set_yticklabels('')
        ylabel('')



dat.salinity.mean(dim='date')

figure(figsize=(8,5))
meanplot2('uacross','uacross')
ylim([400, 0])
axvline(12.5,zorder=2,color='r',linewidth=3)
text(15,25,'min vel',color='red',fontsize=18)
text(25,70,'34',color='orange',fontsize=20)
contour(dat.distance,dat.depth,dat.salinity.mean(dim='date').T,[34],colors='orange',linewidths=3)
title('')
savefig('../figures/irminger17/fullmean_uacross_cczoom.pdf',bbox_inches='tight')


figure(figsize=(8,5))
meanplot2('uacross','uacross')
ylim([400, 0])
axvline(12.5,zorder=2,color='r',linewidth=3)
text(15,25,'min vel',color='red',fontsize=18)
text(48,125,'34.8',color='orange',fontsize=20)
contour(dat.distance,dat.depth,dat.salinity.mean(dim='date').T,[34.8],colors='orange',linewidths=3)
title('')
savefig('../figures/irminger17/fullmean_uacross_egzoom.pdf',bbox_inches='tight')

figure(figsize=(8,5))
meanplot2('uacross','uacross')
text(-5,350,'EGCC',color='white',fontsize=24,zorder=1000)
text(40,350,'EGC/IC',fontsize=24,zorder=1000)
title('Mean across track velocity \n (2014-2016)',fontsize=25)
savefig('../figures/irminger17/fullmean_uacross_fulldpth.pdf',bbox_inches='tight')


figure(figsize=(8,5))
# ax1=subplot(121)
meanplot2('uacross','uacross')
text(-5,600,'Coastal\n Current',color='white',fontsize=24,zorder=1000)
text(40,600,'East Greenland \n Current',fontsize=24,zorder=1000)
# plotinstpos(ax1,'sal')
savefig('../figures/irminger17/fullmean_uacross_forposter.pdf',bbox_inches='tight')

figure(figsize=(8,5))
# ax1=subplot(121)
meanplot2('sal','sal',hq=[34,34.8,34.95])
# plotinstpos(ax1,'sal')
savefig('../figures/irminger17/fullmean_sal_forposter.pdf',bbox_inches='tight')

figure(figsize=(8,5))
# ax2=subplot(122)
meanplot2('tmp','tmp')
# plotinstpos(ax2,'sal')
title('Mean potential temperature',fontsize=18)
savefig('../figures/irminger17/fullmean_tmp_forposter.pdf',bbox_inches='tight')

figure(figsize=(13,4))
ax1=subplot(121)
meanplot2('sal','pden')
plotinstpos(ax1,'sal')
ax2=subplot(122)
meanplot2('tmp','pden',1)
plotinstpos(ax2,'sal')
title('Mean potential temperature',fontsize=18)
savefig('../figures/irminger17/fullmean_saltmp_fulldpth.pdf',bbox_inches='tight')

####################################################################
#################### SHIPBOARD ADCP COMPARISON  #################################
#################################################################################

adcp_dist=array(pd.DataFrame.from_csv(datadir+'Shipboard/adcp_distances.dat').index)

with open(datadir+'Shipboard/adcp_distances.dat', 'r') as f:
    reader = csv.reader(f)
    adcp_dist = list(reader)

adcp_dist=[float(dd[0]) for dd in adcp_dist]

adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')
adcp_16=io.loadmat(datadir+'Shipboard/ar07_2016/ar07_2016_vm_adcp_ilebras.mat')

def acrosstrackname(adic):
    adic['across track velocity']=adic['v_dt']*sin(theta)+adic['u_dt']*cos(theta)

    return adic

adcp_14=acrosstrackname(adcp_14)
adcp_16=acrosstrackname(adcp_16)

len(adcp_14['depth'][:,0])
len(adcp_16['depth'][:,0])

mid_dist_adcp=hstack((diff(adcp_dist)[0]/2,(diff(adcp_dist)[:-1]+diff(adcp_dist)[1:])/2,diff(adcp_dist)[-1]/2))
middistmat_a14=tile(mid_dist_adcp,[len(adcp_14['depth'][:,0]),1])
middistmat_a16=tile(mid_dist_adcp[:-3],[len(adcp_16['depth'][:,0]),1])
shape(middistmat_a16)
depthdiffmat_a14=tile(diff(hstack((0,adcp_14['depth'][:,0]))),[len(adcp_dist),1]).T
depthdiffmat_a16=tile(diff(hstack((0,adcp_16['depth'][:,0]))),[len(adcp_dist[:-3]),1]).T
shape(depthdiffmat_a16)

# minind[tt]=where(max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt])[0][0]+2

plot(adcp_dist[-20:],adcp_14['across track velocity'][0,-20:]);
plot(adcp_dist[-23:-3],adcp_16['across track velocity'][0,-20:]);
plot(adcp_dist[-13:-10],adcp_16['across track velocity'][0,-10:-7]);

minind14=where(adcp_14['across track velocity'][0,-20:]==max(adcp_14['across track velocity'][0,-20:]))[0][0]+len(adcp_14['across track velocity'][0,:])-20
minind16=where(adcp_16['across track velocity'][0,-10:-7]==max(adcp_16['across track velocity'][0,-10:-7]))[0][0]+len(adcp_16['across track velocity'][0,:])-10

minind14
minind16
trans14_less=nansum(adcp_14['across track velocity'][:,minind14:-3]*depthdiffmat_a14[:,minind14:-3]*middistmat_a14[:,minind14:-3])/1e3
trans14=nansum(adcp_14['across track velocity'][:,minind14:]*depthdiffmat_a14[:,minind14:]*middistmat_a14[:,minind14:])/1e3
trans16=nansum(adcp_16['across track velocity'][:,minind16:]*depthdiffmat_a16[:,minind16:]*middistmat_a16[:,minind16:])/1e3

trans14
trans14_less
trans16
#new function just has minimal framework
def contadcp_min(axspec,distvec,adic):
    axspec.contourf(distvec,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(distvec,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)

#same for mooring version
def contmoor_min(axspec,datechoice):
    im=axspec.contourf(dat.distance,dat.depth,dat['across track velocity'][:,:,datechoice].T,arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(dat.distance,dat.depth,dat['across track velocity'][:,:,datechoice].T,[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='w',linewidth=2) for mm in distvec]
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return im



def plot4pan():
    fs=18
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(10,7.5))
    contadcp_min(ax1,adcp_dist,adcp_14)
    ax1.text(-5,600,'-0.97 Sv',color='white',zorder=1000,fontsize=18)
    ax1.text(-5,725,'-1.35 Sv',color='white',zorder=1000,fontsize=18)
    ax2.text(-5,600,'-0.46 Sv',color='white',zorder=1000,fontsize=18)
    ax1.set_xlim([-15,110])
    ax1.set_ylim([800,0])
    im=contmoor_min(ax2,0)
    plotinstpos(ax2,'v')
    mkbox(ax1,-12,adcp_dist[minind14],150,5)
    mkbox(ax2,-5,dat.distance[minind[0]],150,5)
    contadcp_min(ax3,adcp_dist[:-3],adcp_16)
    im=contmoor_min(ax4,-1)
    mkbox(ax3,-5,adcp_dist[minind16],150,5)
    mkbox(ax4,-5,daily.distance[minind[-1]],150,5)
    plotinstpos(ax4,'v')
    ax3.text(-5,600,'-0.76 Sv',color='white',zorder=1000,fontsize=18)
    ax4.text(-5,600,'-0.80 Sv',color='white',zorder=1000,fontsize=18)

    ax1.set_title('Vessel mounted ADCP',fontsize=fs)
    ax2.set_title('Moored array measurements',fontsize=fs)
    ax1.set_ylabel('August, 2014 \n \n depth (m)',fontsize=fs)
    ax3.set_ylabel('July, 2016 \n \n depth (m)',fontsize=fs)
    # ax1.set_ylabel('depth (m)')
    # ax3.set_ylabel('depth (m)')
    ax3.set_xlabel('distance (km)',fontsize=fs)
    ax4.set_xlabel('distance (km)',fontsize=fs)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax,label='across-track velocity [m/s]',ticks=univec['uacross'][3][::2])

    savefig('../figures/irminger17/ADCPmoorcomp_4panel.pdf',bbox_inches='tight')


plot4pan()


####################################################################
#################### SEASONAL SECTION  #################################
#################################################################################
titvec=['August - October','November - January','February - April','May - July']

def seasplot(field,field2):
    fig=figure(figsize=(10,7.5))
    for num in range(4):
        subplot(2,2,num+1)
        elfieldo=(dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,num]+dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,num+4]).T/2
        elfieldo2=(dat[univec[field2][0]].resample('3M',how='mean',dim='date')[:,:,num]+dat[univec[field2][0]].resample('3M',how='mean',dim='date')[:,:,num+4]).T/2
        ax1=overcont(elfieldo,elfieldo,titvec[num],univec[field][1],univec[field][2],
        univec[field][3],nomoorlines=1)
        contour(dat.distance,dat.depth,elfieldo2,[34,34.8],colors='orange',linewidths=3)
        axvline(12,zorder=2,color='r',linewidth=3)
        ylim([400, 0])

    fig.subplots_adjust(hspace=0.15,wspace=0.1)
    subplot(221)
    gca().set_xticklabels('')
    xlabel('')
    subplot(222)
    gca().set_xticklabels('')
    gca().set_yticklabels('')
    xlabel('')
    ylabel('')
    subplot(224)
    gca().set_yticklabels('')
    ylabel('')

    fig.subplots_adjust(right=0.8)
    cbaxes = fig.add_axes([0.85, 0.15, 0.025, 0.6])
    cbar=colorbar(ax1, cax = cbaxes,ticks=univec[field][3][::2],label=univec[field][4])

    savefig('../figures/irminger17/seasonal_'+field+'.pdf',bbox_inches='tight')



seasplot('uacross','sal')
