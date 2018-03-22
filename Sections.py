#################################################################################
#################################################################################
#################################################################################
######################## Make sections ###################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


dat=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1803bathy.pickle','rb'))

# bathf=interpolate.interp1d(bathdist,bathbath)
# bathonmygrid=bathf(dat.distance)
# justsal=dat['salinity'].copy()
# justtmp=dat['temperature'].copy()
# for dd,adist in enumerate(dat.distance):
#             justsal[dd,:,:]=justsal[dd,:,:].where(justsal[dd,:,:].depth<=bathonmygrid[dd])
#             justtmp[dd,:,:]=justtmp[dd,:,:].where(justtmp[dd,:,:].depth<=bathonmygrid[dd])

#################################################################################
# Make mean and var sections of all fields
#################################################################################


def onecont(field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,addcont=0):
    if axx==0:
        axx=subplot(111)
    ax1=axx.contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    # if 'sal' in tit:
    #     axx.contour(dat.distance,dat.depth,dat['salinity'].mean(dim='date').T,levels=[34.9],colors='w',linewidths=12)
    if addcont==1:
        ax2=contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
        clabel(ax2,fmt='%1.1f')
    axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance [km]',fontsize=18)
    ylabel('depth [m]',fontsize=18)
    xlim([-5,100])
    ylim([2200,0])
    title(tit,fontsize=22)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axx.set_yticks(range(500,2100,500))
    return ax1


def partmean(date1,date2,textit,axx=0,cbarito=1,cwidth=0.025,nomoorlines=1):
    field='uacross'
    if axx==0:
        axx=subplot(111)
    axvel=onecont(dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T,'',univec[field][1],univec[field][2],univec[field][3],axx,nomoorlines=nomoorlines)

    axsal=contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[33.5,34,34.4,34.95],cmap=cm.YlOrRd,linewidths=2)
    manual_locations = [(10,100), (30, 100), (80,500)]
    clabel(axsal,fmt='%1.2f',manual=manual_locations)
    ax349=contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.9],colors='w',linewidths=12)
    ax349=contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.9],colors='k',linewidths=4)
    clabel(ax349,fmt='%1.1f',manual=[(50,250)])
    if cbarito==1:
        cbaxes = fig.add_axes([0.95, 0.15, cwidth, 0.7])
        cbar=colorbar(axvel,ticks=univec[field][3],label='[m/s]',cax=cbaxes)
    axx.text(0,1250,textit,color='white',fontsize=20,zorder=200)

    return axx

def plotbox(ax,x1,x2,y1,y2,colo):
    ax.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=colo,linewidth=4,zorder=300)

def plotkinkybox(ax,x1,x2,x3,y1,y2,y3,colo):
    ax.plot([x1,x2,x3,x3,x2,x2,x1,x1],[y1,y1,y1,y3,y3,y2,y2,y1],color=colo,linewidth=4,zorder=300)


fig=figure(figsize=(16,12))
ax1=subplot(221)
partmean('2014-12-01','2015-1-01','Dec 2014 \nCC offshore?',ax1,nomoorlines=0)
ax2=subplot(222)
axx=partmean('2015-08-01','2015-10-01','Aug/Sep 2015 \nCC offshore?',ax2,nomoorlines=0)
ax3=subplot(223)
partmean('2015-03-01','2015-06-01','Mar-May 2015 \nCC inshore?',ax3,nomoorlines=0)
ax4=subplot(224)
axx=partmean('2016-04-01','2016-05-01','Apr 2016 \nCC inshore?',ax4,nomoorlines=0)
savefig('../figures/xport/CurrentMotionTest.pdf',bbox_inches='tight')

fig=figure(figsize=(16,12))
ax1=subplot(221)
partmean('2015-05-01','2015-06-01','May 2015 \nCC offshore?',ax1,nomoorlines=0)
ax2=subplot(222)
axx=partmean('2016-04-01','2016-05-01','Apr 2016 \nCC offshore?',ax2,nomoorlines=0)
ax3=subplot(223)
partmean('2014-09-01','2014-10-01','Sep 2014 \nCC inshore?',ax3,nomoorlines=0)
ax4=subplot(224)
axx=partmean('2015-10-01','2015-11-01','Oct 2015 \nCC inshore?',ax4,nomoorlines=0)
savefig('../figures/xport/CurrentMotionTest2.pdf',bbox_inches='tight')


# fig=figure(figsize=(8,6))
# axx=partmean('2014-1-01','2016-12-01','Mean 2014-2016')
# axx.text(-20,-50,'Coastal current',fontsize=20)
# axx.text(30,-50,'Slope current',fontsize=20)
# axx.text(40,150,'EGC',fontsize=20)
# axx.text(60,400,'IC',fontsize=20)
# axx.text(85,1900,'DWBC',fontsize=20)
# savefig('../../confschools/1802_oceansciences/presentation/figures/velpropmean.pdf',bbox_inches='tight')
# # plotbox(axx,-4.5,15,20,200,'r')
# # savefig('../../confschools/1802_oceansciences/presentation/figures/velpropmean_wcc.pdf',bbox_inches='tight')
# plotkinkybox(axx,15,35.5,99.5,20,200,1750,'red')
# savefig('../../confschools/1802_oceansciences/presentation/figures/velpropmean_wegic.pdf',bbox_inches='tight')
#
#
#
# fig=figure(figsize=(16,6))
# ax1=subplot(121)
# partmean('2014-10-01','2015-1-1','FALL',ax1,0)#,'Coastal current is fastest\n\nTS suggests local origin\n\nStill high salinity off-shore')
# ax1.text(35,150,'EGC',fontsize=20)
# ax1.text(60,400,'IC',fontsize=20)
# ax2=subplot(122)
# partmean('2015-1-1','2015-4-1','WINTER',ax2,1,0.015)#,'Shelf freshens\n\nLow salinity off-shore\n\nTS shows mixed PW')
# ax2.text(45,200,'EGC',fontsize=20)
# ax2.text(60,500,'IC',fontsize=20)
# ax2.set_yticklabels('')
# ax2.set_ylabel('')
# savefig('../../confschools/1802_oceansciences/presentation/figures/fallvswintersections.pdf',bbox_inches='tight')


# partmean('2015-6-1','2015-9-1','')
# partmean('2015-9-1','2015-12-1','')
# partmean('2015-12-1','2016-3-1','')
# partmean('2016-5-1','2016-8-1','')


def plotinstpos(axchoice,savename):
        if 'uac' in savename:
            for rr in range(7):
                if rr==0:
                    axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,zorder=40,label='ADCP')
                else:
                    axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,zorder=40)
        mm=0
        for key in depths:
            for dd in range(len(depths[key])):
                if ('sal' in savename) | ('tmp' in savename) | ('pden' in savename):
                    if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
                        axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=4)
                    # if 'tmp' in savename:
                    #     if ('tidbit' in inst[key][dd]) | ('olo' in inst[key][dd]):
                    #          axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35)

                elif 'uac' in savename:
                    if ('AQ' in inst[key][dd]):
                        if dd==0:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='Current meter')
                        else:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6)
            mm+=1

univec['tmp'][-1]
def meanOSM(field,textit,ax1=0,ypos=1250,cwidth=0.025,cpos=0.95):
    if ax1==0:
        ax1=subplot(111)
    axvel=onecont(dat[univec[field][0]].mean(dim='date').T,field,univec[field][1],univec[field][2],univec[field][3],ax1,addcont=1,nomoorlines=0)
    cbaxes = fig.add_axes([cpos, 0.15,cwidth , 0.7])
    colorbar(axvel,extend='both',label=univec[field][-1],cax=cbaxes)
    ax1.text(0,ypos,textit,color='white',fontsize=20,zorder=100)
    ax1.set_title('')
    plotinstpos(ax1,field)
    if field=='uacross':
        ax1.legend(loc=(0.05,0.08),fontsize=18).set_zorder(100)
    savefig('../figures/xport/'+field+'_mean_winst.png',bbox_inches='tight')
    savefig('../../confschools/1802_oceansciences/presentation/figures/'+field+'_mean_winst.pdf',bbox_inches='tight')

    return ax1

fig=figure(figsize=(8,6))
meanOSM('uacross','Mean velocity\n\n2014-2016')
fig=figure(figsize=(16,6))
ax1=subplot(121)
meanOSM('sal','Salinity',ax1,ypos=1750,cwidth=0.015,cpos=0.49)
ax2=subplot(122)
ax2=meanOSM('tmp','Pot. Temperature',ax2,ypos=1750,cwidth=0.015,cpos=0.915)
ax2.set_yticks([])
ax2.set_yticklabels('')
ax2.set_ylabel('')
savefig('../../confschools/1802_oceansciences/presentation/figures/saltmp_mean_winst.pdf',bbox_inches='tight')


hexbin(justsal.values.flatten(),justtmp.values.flatten(),cmap=cm.Greys,bins='log')

def partTS(date1,date2):
    figure(figsize=(5,4))
    ax2 = subplot(111)
    # ax2.plot(justsal.values.flatten(),justtmp.values.flatten(),'.',color='grey')
    # ax2.plot(justsal.sel(date=slice(date1,date2)).values.flatten(),justtmp.sel(date=slice(date1,date2)).values.flatten(),'r.')
    # ax2.hexbin(justsal.values.flatten(),justtmp.values.flatten(),cmap=cm.Greys,bins='log')
    ax2.hexbin(justsal.sel(date=slice(date1,date2)).values.flatten(),justtmp.sel(date=slice(date1,date2)).values.flatten(),
    cmap=cm.hot_r,bins='log',mincnt=1)
    ax2.set_xlim([32,35.2])
    ax2.set_ylim([-2,8])
    ax2.set_xlabel('salinity')
    ax2.plot([32.5,32.5,34.25,34.25,32.5],[-1.8,0,0,-1.8,-1.8],color='k')
    ax2.contour(salvec,tmpvec,pdenmat,colors='darkgrey',zorder=50)
    ax2.text(34.95,6,'AW',color='k',fontsize=16)
    ax2.text(34.35,-1,'PW',color='k',fontsize=16)
    ax2.set_ylabel('pot. temperature [$^\circ$C]')
    ax2.axvline(34.9,color='k')

partTS('2014-10-01','2015-1-1')
savefig('../../confschools/1802_oceansciences/presentation/figures/TS_fall_hex.pdf',bbox_inches='tight')

partTS('2015-01-01','2015-4-1')
savefig('../../confschools/1802_oceansciences/presentation/figures/TS_winter_hex.pdf',bbox_inches='tight')



def meanplot(field):
    # figure()
    ax1,ax2=onecont(dat[univec[field][0]].mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3])
    colorbar(ax1,ticks=univec[field][3])
    # xlabel('')
    savefig('../figures/sections/fullmean_'+field+'_fulldpth.png',bbox_inches='tight')
    savefig('../figures/sections/fullmean_'+field+'_fulldpth.pdf',bbox_inches='tight')


def stdplot(field,limoverride):

    figure()
    ax1,ax2=onecont(dat[univec[field][0]].std(dim='date').T,'Standard deviation of '+univec[field][0],limoverride,cm.Purples,limoverride[::5])
    colorbar(ax1,label=univec[field][4])
    savefig('../figures/sections/fullstd_'+field+'_fulldpth.png',bbox_inches='tight')
    savefig('../figures/sections/fullstd_'+field+'_fulldpth.pdf',bbox_inches='tight')

titvec=['August - October','November - January','February - April','May - July']

def fourplot(field,secondyear=0):
    fig=figure(figsize=(8,6))
    for num in range(4):
        if secondyear==1:
            datenum=num+4
            datename='second'
        elif secondyear==0:
            datenum=num
            datename='first'
        else:
            datename='both'

        subplot(2,2,num+1)
        if datename=='both':
            elfieldo=(dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,num]+dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,num+4]).T/2
        else:
            elfieldo=dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,datenum].T
        ax1,ax2=onecont(elfieldo,titvec[num],univec[field][1],univec[field][2],
        univec[field][3],nomoorlines=1)

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
    cbar=colorbar(ax1, cax = cbaxes,ticks=univec[field][3],label=univec[field][4])
    suptitle(univec[field][0],fontsize=20)
    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year.png')
    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year.pdf')
    # figure(figsize=(12,6))
    for nn in range(1,5):
        subplot(2,2,nn)
        ylim([400,0])
        clabel(ax2)

    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year_zoom.png')
    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year_zoom.pdf')

def plot12all(field):
    fourplot(field)
    fourplot(field,1)
    fourplot(field,4)


def seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='both'):
    if ychoice=='both':
        sind=0
        fin=8
    elif ychoice=='first':
        sind=0
        fin=4
    elif ychoice=='second':
        sind=4
        fin=8
    figure()
    for ii in range(sind,fin):
        if atype=='1depth':
            dathere=dat[field].resample('3M',dim='date',how='mean')[:,int(dpchoice/2),ii]
            disthere=dat.distance
        elif atype=='diff':
            dathere=diff(dat[field].resample('3M',dim='date',how='mean')[:,int(dpchoice/2),ii],axis=0)[1:-1]
            disthere=dat.distance[1:-2]+diff(dat.distance[1:-1])/2
        elif atype=='mean':
            dathere=mean(dat[field].resample('3M',dim='date',how='mean')[:,:int(dpchoice/2),ii],axis=1)
            disthere=dat.distance
        if ii>=4:
            lsty='--'
            lab=''
        else:
            lsty='-'
            lab=titvec[ii]
        if ychoice=='second':
            lab=titvec[ii-4]
        plot(disthere,dathere,label=lab,linestyle=lsty,color=colvec[mod(ii,4)])
    plot(-5,mean(dathere),'k-',label='First year, 2014-2015')
    plot(-5,mean(dathere),'k--',label='Second year, 2015-2016')
    ylabel(ylab)
    xlabel('distance (km)')
    title(field+tit)
    legend(loc=(1.05,0.2))
    xlim([-5,100])
    savefig('../figures/seasonal/seasonline_'+savename+'_'+ychoice+'.png',bbox_inches='tight')


def seasonline_all(field,dpchoice,ylab,tit,atype,savename,ychoice='both'):
    seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='both')
    seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='first')
    seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='second')

dat['cross-slope heat flux']=dat['along track velocity']*dat['temperature']
dat['cross-slope salt flux']=dat['along track velocity']*dat['salinity']

plot(dat['cross-slope salt flux'].mean(dim='date').T);

univec['uT']=['cross-slope heat flux',arange(-0.4,0.4,0.01),cm.RdBu_r,arange(-0.4,0.4,0.2),'[m$^\circ$C/s]']

univec['uS']=['cross-slope salt flux',arange(-4,4,0.1),cm.RdBu_r,arange(-4,4,1),'[m/s]']

meanplot('uT')

meanplot('uacross')

meanplot('tmp')


stdplot('ualong',arange(0,0.3,0.01))


stdplot('uT',arange(0,0.7,0.05))

univec['tmp'][3]

figure(figsize=(20,7))
ax1=subplot(121)
meanplot('uacross')
plotinstpos(ax1,'v')
title('Mean velocity [m/s]',fontsize=22)
text(0,1750,'Shelf of \nGreenland',color='white',zorder=120,fontsize=26)
ax2=subplot(122)
meanplot('tmp')
title('Mean temperature [$^\circ$ C]',fontsize=22)
plotinstpos(ax2,'sal')
savefig('../figures/sections/tmpvel.pdf',bbox_inches='tight')

stdplot('uacross',arange(0,0.3,0.01))

meanplot('pden')
meanplot('sal')


meanplot('tmp')

stdplot('pden',arange(0,0.4,0.01))
stdplot('sal',arange(0,0.5,0.01))
stdplot('tmp',arange(0,2,0.1))

#################################################################################
# Make seasonal sections of all fields
#################################################################################

plot12all('uacross')

plot12all('ualong')

plot12all('sal')
plot12all('tmp')

plot12all('pden')



#################################################################################
# Line plots of seasonal evolution
#################################################################################

colvec=['#e41a1c','#ff7f00','#377eb8','#4daf4a']

dat


seasonline_all('potential density',0,'[kg/m$^3$]',' difference at surface','diff','diffpdensurf')
seasonline('potential density',0,'[kg/m$^3$]',' difference at surface','diff','diffpdensurf','first')
seasonline('potential density',0,'[kg/m$^3$]',' difference at surface','diff','diffpdensurf','second')

seasonline_all('potential density',100,'[kg/m$^3$]',' difference at 100m','diff','diffpden100')
seasonline_all('across track velocity',0,'[m/s]',' at surface','1depth','surfvel')

seasonline_all('across track velocity',100,'[m/s]',' at 100m','1depth','surfvel')
seasonline_all('across track velocity',100,'[m/s]',' averaged over top 100m','mean','mveltop100')


#################################################################################
# Save snapshots to make a movie of velocity
#################################################################################

# def moviedaily(field):
#     dat.where(dat[univec[field][0]]<-1.5)[univec[field][0]]=-1.5
#     for ii,na in enumerate(dat.date):
#         ax1,ax2=onecont(dat[univec[field][0]][:,:,ii].T,univec[field][0],arange(-1.5,1.5,0.05),univec[field][2],arange(-1.8,1.8,0.2),nomoorlines=1)
#         colorbar(ax1,label=univec[field][4],ticks=arange(-1.5,1.5,0.5))
#         t= pd.to_datetime(str(na.data))
#         title(t.strftime('%m - %d - %Y'))
#         savefig('../figures/vmovie/fulldepth/vel_'+str("%03d"%ii)+'_fulldpth.png',bbox_inches='tight',dpi=300)
#         clf()
#
#
#
# def moviedailypden(field):
#     for ii,na in enumerate(dat.date):
#         ax1=onecont(dat[univec[field][0]][:,:,ii].T,univec[field][0],univec[field][0],univec[field][1],univec[field][3],nomoorlines=1)
#         colorbar(ax1,label=univec[field][4],ticks=arange(26,28,0.5))
#         t= pd.to_datetime(str(na.data))
#         title(t.strftime('%m - %d - %Y'))
#         savefig('../figures/vmovie/fulldepth/pden_'+str("%03d"%ii)+'_fulldpth.png',bbox_inches='tight',dpi=300)
#         clf()
#
#
# def moviesnaps(field):
#     dat.where(dat[univec[field][0]]<-1.5)[univec[field][0]]=-1.5
#     for ii,na in enumerate(dat.date[::7]):
#         ax1=onecont(dat[univec[field][0]].resample('W',dim='date',how='mean')[:,:,ii].T,univec[field][0],arange(-1.5,1.5,0.05),univec[field][2],arange(-1.8,1.8,0.2),nomoorlines=1)
#         colorbar(ax1,label=univec[field][4],ticks=arange(-1.5,1.5,0.5))
#         t= pd.to_datetime(str(na.data))
#         title(t.strftime('%m - %d - %Y'))
#         ylim([800,0])
#         savefig('../figures/vmovie/weekly/vel_weekly_'+str("%03d"%ii)+'.png',bbox_inches='tight',dpi=300)
#         clf()


#################################################################################
# Make sections of velocity for first and last time point to compare with shipboard
#################################################################################
# Note: making a four-panel figure in Shipboard_ADCP.py for comparison sake.
#
#
# def contmoorvel(datechoice,tit):
#     contourf(dat.distance,dat.depth,dat['across track velocity'][:,:,datechoice].T,arange(-0.8,0.8,0.05),cmap=cm.RdBu_r);
#     colorbar(label='m/s')
#     contour(dat.distance,dat.depth,dat['across track velocity'][:,:,datechoice].T,[-0.4,-0.2],colors='k');
#     fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
#     gca().invert_yaxis()
#     title('Mooring measured across-track velocity, '+tit)
#     ylim([800,0])
#     xlabel('distance (km)')
#     ylabel('depth (m)')
#     xlim([-20,90])
#     savefig('../figures/shipboard/Mooringsection_forADCPcomp_'+tit+'.pdf')
#
# contadcp(0,'2014')
#
# contadcp(-1,'2016')
