#################################################################################
#################################################################################
#################################################################################
######################## Make sections ###################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

23.54+3.36
1.80+0.12

newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1804shelf.pickle','rb'))

savename='_1804shelf'

dat=newgrid.copy()

# bathf=interpolate.interp1d(bathdist,bathbath)
# bathonmygrid=bathf(newgrid.distance)

# dat=newgrid.copy()
# for vv in newgrid:
#     if vv[0]!='d':
#         print(vv)
#         for dd,adist in enumerate(dat.distance):
#             dat[vv][dd,:,:]=dat[vv][dd,:,:].where(dat[vv][dd,:,:].depth<=(bathonmygrid[dd]-2)) # limit so that labels are more likely to be in the right place, but still covering gaps

bathdistnew=hstack((-15,-11,-10,bathdist))
bathbathnew=hstack((0,0,100,bathbath))


#################################################################################
# Plotting tools
#################################################################################

def plotinstpos(axchoice,savename):
        if ('uac' in savename) | ('v' in savename):
            for rr in range(7):
                if rr==0:
                    axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,label='ADCP')
                else:
                    axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12)
        mm=0
        for key in depths:
            for dd in range(len(depths[key])):
                if ('sal' in savename) | ('tmp' in savename) | ('pden' in savename):
                    if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=4)

                elif 'uac' in savename:
                    if ('AQ' in inst[key][dd]):
                        if dd==0:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=6,label='Current meter')
                        else:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=6)
            mm+=1


def onecont(field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,addcont=1,alone=1):
    if axx==0:
        axx=subplot(111)
    ax1=axx.contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    # if 'sal' in tit:
    #     axx.contour(dat.distance,dat.depth,dat['salinity'].mean(dim='date').T,levels=[34.9],colors='w',linewidths=12)
    if addcont==1:
        ax2=axx.contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
        if ('den' in tit) | ('sal' in tit):
            fstr='%1.2f'
        elif 'temp' in tit:
            fstr='%1.0f'
        if 'vel' in tit:
            manual_locations = [(60, 1200), (60, 50),]
            clabels=clabel(ax2, fmt='%1.1f', fontsize=10, manual=manual_locations)
        else:
            clabels=clabel(ax2,fmt=fstr)
        [txt.set_backgroundcolor('white') for txt in clabels]
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0,alpha=0.8)) for txt in clabels]
        # [txt.set_alpha(0.5) for txt in clabels]
    else:
        ax2=0
    axx.fill_between(bathdistnew,bathbathnew,2500*ones(len(bathbathnew)),color='k',zorder=22)
    if alone==1:
        xlabel('distance [km]',fontsize=18)
        ylabel('depth [m]',fontsize=18)
        title(tit,fontsize=22)
    xlim([-15,100])
    ylim([2200,-20])
    if nomoorlines==0:
        [axx.axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axx.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axx.set_yticks(range(500,2100,500))
    return ax1,ax2


#################################################################################
# Make 4 panel plot of all means for paper with instrument positions
#################################################################################


def eachpan(field,axit,atit):
    ax1,ax2=onecont(dat[univec[field][0]].mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3],axx=axit,alone=0,nomoorlines=0)
    plotinstpos(axit,field)
    if ('sal' in field):
        cticks=[34,  34.8, 34.92,34.96, 35]
    elif ('tmp' in field):
        cticks=range(0,10,2)
    else:
        cticks=univec[field][3]
    colorbar(ax1,ticks=cticks,ax=axit)
    axit.text(-10,2100,atit,color='white',zorder=1e3,fontsize=14)


univec['sal']=['salinity',array([33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35]),sal_cmap,array([33,34, 34.8,34.92,34.94,34.96,34.98, 35]),'']
univec['tmp']=['temperature',linspace(-1,8,31),cm.RdYlBu_r,[2,4,5,6,7,8],'[$^\circ$C]']
def plot4pan_means():
    f, ((ax11, ax22), (ax33, ax44)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(11,7))
    eachpan('uacross',ax11,'a) Velocity [m/s]')
    eachpan('pden',ax22,'b) Potential density [kg/m$^3$]')
    eachpan('sal',ax33,'c) Salinity')
    eachpan('tmp',ax44,'d) Potential temperature [$^\circ$C]')
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
    ax11.text(-20,-100,'Coastal current',fontsize=14)
    ax11.text(35,-100,'Slope current',fontsize=14)
    [ax22.text(distvec[ii]-4,-100,'CF'+str(ii+1),fontsize=13) for ii in range(7)]
    ax22.text(distvec[-1]-4,-100,'M1',fontsize=13)
    plt.tight_layout()
    savefig('../figures/paperfigs/meansections.pdf',bbox_inches='tight')


plot4pan_means()

def partmean(date1,date2,textit,axx):
    field='uacross'
    axvel,na=onecont(dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T,'vel',univec[field][1],
    univec[field][2],univec[field][3],axx,addcont=0,alone=0)
    ax34=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.2],colors='orange',linewidths=2)
    ax349=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.9],colors='k',linewidths=4)

    clab34=clabel(ax34,fmt='%1.1f',manual=[(25,100)],colors='k')
    clab349=clabel(ax349,fmt='%1.1f',manual=[(50,250)])
    if 'WINTER' not in textit:
        ax3495=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.95],colors='darkred',linewidths=2)
        clab3495=clabel(ax3495,fmt='%1.2f',manual=[(80,500)],colors='k')
        clablist=[clab34,clab349,clab3495]
    else:
        clablist=[clab34,clab349]
    for clabels in clablist:
        [txt.set_backgroundcolor('white') for txt in clabels]
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0,alpha=0.8)) for txt in clabels]
    axx.text(-10,2100,textit,color='white',fontsize=20,zorder=200)

    return axvel


def plot4vel_seas():
    f, ((ax11, ax22), (ax33, ax44)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(9.5,7))
    partmean('2014-10-01','2015-1-1','FALL',ax11)
    partmean('2015-01-01','2015-04-1','WINTER',ax22)
    partmean('2015-04-01','2015-7-1','SPRING',ax33)
    axvel=partmean('2015-7-01','2015-10-1','SUMMER',ax44)
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
    plt.tight_layout()
    cbaxes = f.add_axes([1.01, 0.15, 0.025, 0.7])
    cbar=colorbar(axvel,ticks=univec['uacross'][3],label='Velocity [m/s]',cax=cbaxes)
    savefig('../figures/paperfigs/vel_seasonal.pdf',bbox_inches='tight')


plot4vel_seas()


def plot4vel_seas_symm():
    f, (ax11, ax22, ax33) = plt.subplots(3,1, sharex=True, sharey=True,figsize=(5,10))
    partmean('2014-10-01','2015-1-1','FALL',ax11)
    partmean('2015-01-01','2015-04-1','WINTER',ax22)
    axvel=partmean('2015-04-01','2015-9-1','SPRING/SUMMER',ax33)
    # axvel=partmean('2015-7-01','2015-10-1','SUMMER',ax44)
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
    plt.tight_layout()
    cbaxes = f.add_axes([1.0, 0.3, 0.06, 0.45])
    cbar=colorbar(axvel,ticks=univec['uacross'][3],label='Velocity [m/s]',cax=cbaxes)
    savefig('../figures/paperfigs/vel_seasonal_symm.pdf',bbox_inches='tight')


plot4vel_seas_symm()

XXXXXXXXXXXXXX
# def meanOSM(field,textit,ax1=0,ypos=1250,cwidth=0.025,cpos=0.95):
#     if ax1==0:
#         ax1=subplot(111)
#     axvel=onecont(dat[univec[field][0]].mean(dim='date').T,field,univec[field][1],univec[field][2],univec[field][3],ax1,addcont=1,nomoorlines=0)
#     cbaxes = fig.add_axes([cpos, 0.15,cwidth , 0.7])
#     colorbar(axvel,extend='both',label=univec[field][-1],cax=cbaxes)
#     ax1.text(0,ypos,textit,color='white',fontsize=20,zorder=100)
#     ax1.set_title('')
#     plotinstpos(ax1,field)
#     if field=='uacross':
#         ax1.legend(loc=(0.05,0.08),fontsize=18).set_zorder(100)
#
#
#
#     return ax1
#
#
#
# fig=figure(figsize=(8,6))
# meanOSM('uacross','Mean velocity\n\n2014-2016')
# fig=figure(figsize=(16,6))
# ax1=subplot(121)
# meanOSM('sal','Salinity',ax1,ypos=1750,cwidth=0.015,cpos=0.49)
# ax2=subplot(122)
# ax2=meanOSM('tmp','Pot. Temperature',ax2,ypos=1750,cwidth=0.015,cpos=0.915)
# ax2.set_yticks([])
# ax2.set_yticklabels('')
# ax2.set_ylabel('')
#

#################################################################################
# Make simple plots for reference
#################################################################################

def meanplot(field):
    figure()
    ax1,ax2=onecont(dat[univec[field][0]].mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3])
    colorbar(ax1,ticks=univec[field][3])
    # xlabel('')
    savefig('../figures/sections/fullmean_'+field+'_fulldpth'+savename+'.png',bbox_inches='tight')
    # savefig('../figures/sections/fullmean_'+field+'_fulldpth'+savename+'.pdf',bbox_inches='tight')
#
#
# def stdplot(field,limoverride):
#
#     figure()
#     ax1,ax2=onecont(dat[univec[field][0]].std(dim='date').T,'Standard deviation of '+univec[field][0],limoverride,cm.Purples,limoverride[::5])
#     colorbar(ax1,label=univec[field][4])
#     savefig('../figures/sections/fullstd_'+field+'_fulldpth'+savename+'.png',bbox_inches='tight')
#     # savefig('../figures/sections/fullstd_'+field+'_fulldpth'+savename+'.pdf',bbox_inches='tight')
#
#
#
# meanplot('uacross')
# meanplot('ualong')
# meanplot('pden')
meanplot('sal')

# meanplot('tmp')
#
#
# stdplot('ualong',arange(0,0.3,0.01))
# stdplot('uacross',arange(0,0.3,0.01))
# stdplot('pden',arange(0,0.4,0.01))
# stdplot('sal',arange(0,0.5,0.01))
# stdplot('tmp',arange(0,2,0.1))
#
