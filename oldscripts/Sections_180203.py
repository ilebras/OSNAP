#################################################################################
#################################################################################
#################################################################################
######################## Make sections ###################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


dat=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1801.pickle','rb'))

#################################################################################
# Make mean and var sections of all fields
#################################################################################


def onecont(field,tit,vrange,coloor,hlevs,nomoorlines=1):
    ax1=contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    ax2=contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
    # clabel(ax2)
    fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance (km)',fontsize=18)
    ylabel('depth (m)',fontsize=18)
    xlim([-5,100])
    ylim([2200,0])
    title(tit)
    #
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return ax1,ax2

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

seasonline_all('across track velocity',200,'[m/s]',' averaged over top 200m','mean','mveltop200')

seasonline_all('temperature',100,'',' difference at 100m','diff','difftmp')
seasonline_all('salinity',100,'',' difference at 100m','diff','diffsal')
seasonline_all('temperature',100,'',' at 100m','1depth','100tmp')
seasonline_all('salinity',100,'',' at 100m','1depth','100sal')



xxxxx

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
