#################################################################################
#################################################################################
#################################################################################
######################## Make sections ###################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *


dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))
dat['across track velocity']=-1*dat['across track velocity']

#################################################################################
# Make mean and var sections of all fields
#################################################################################
bathdistnew=hstack((-15,-11,-10,bathdist))
bathbathnew=hstack((0,0,100,bathbath))

def onecont(field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,addcont=0):
    if axx==0:
        axx=subplot(111)
    ax1=axx.contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    # if 'sal' in tit:
    #     axx.contour(dat.distance,dat.depth,dat['salinity'].mean(dim='date').T,levels=[34.9],colors='w',linewidths=12)
    if addcont==1:
        ax2=contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
        clabel(ax2,fmt='%1.1f')
    axx.fill_between(bathdistnew,bathbathnew,2500*ones(len(bathbathnew)),color='k',zorder=22)
    xlabel('distance [km]',fontsize=18)
    ylabel('depth [m]',fontsize=18)
    xlim([-15,95])
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
    axx.text(-10,1250,textit,color='white',fontsize=20,zorder=200)

    return axx

def plotbox(ax,x1,x2,y1,y2,colo):
    ax.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=colo,linewidth=4,zorder=300)

def plotkinkybox(ax,x1,x2,x3,y1,y2,y3,colo):
    ax.plot([x1,x2,x3,x3,x2,x2,x1,x1],[y1,y1,y1,y3,y3,y2,y2,y1],color=colo,linewidth=4,zorder=300)





fig=figure(figsize=(8,6))
axx=partmean('2014-1-01','2016-12-01','Mean 2014-2016')
axx.text(-20,-50,'Coastal current',fontsize=20)
axx.text(30,-50,'Slope current',fontsize=20)
axx.text(40,150,'EGC',fontsize=20)
axx.text(60,400,'IC',fontsize=20)
# savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/velpropmean.pdf',bbox_inches='tight')
# plotbox(axx,-12,15,20,200,'r')
# savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/velpropmean_wcc.pdf',bbox_inches='tight')
plotkinkybox(axx,15,35.5,99.5,20,200,1750,'red')
axx.text(80,1900,'DWBC',fontsize=20)
axx.text(-10,1750,'$\sigma_0 = 27.8 kg/m^3$',color='white')
savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/velpropmean_wegic.pdf',bbox_inches='tight')

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


def partmean_wden(date1,date2,textit,axx=0,cbarito=1,cwidth=0.025,nomoorlines=1):
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
    axx.text(-10,1250,textit,color='white',fontsize=20,zorder=200)

    return axx




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


univec['sal'][-2]= array([33.  , 34.  ,34.4, 34.8, 34.9,34.92, 34.94, 34.96, 34.98, 35.  ])
univec['tmp'][-2]=arange(-1,9,1)

def meanOSM(field,textit,ax1=0,ypos=1250,cwidth=0.025,cpos=0.95):
    if ax1==0:
        ax1=subplot(111)
    axvel=onecont(dat[univec[field][0]].mean(dim='date').T,field,univec[field][1],univec[field][2],univec[field][3],ax1,addcont=1,nomoorlines=0)
    cbaxes = fig.add_axes([cpos, 0.15,cwidth , 0.7])
    colorbar(axvel,extend='both',label=univec[field][-1],cax=cbaxes,ticks=univec[field][-2])
    ax1.text(-10,ypos,textit,color='white',fontsize=20,zorder=100)
    ax1.set_title('')
    plotinstpos(ax1,field)
    if field=='uacross':
        ax1.legend(loc=(0.05,0.08),fontsize=18).set_zorder(100)
    savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/'+field+'_mean_winst.pdf',bbox_inches='tight')

    return ax1

univec['uacross'][2]=cm.Blues

fig=figure(figsize=(8,6))
ax1=meanOSM('uacross','Mean velocity\n\n2014-2016')
ax1.text(-20,-50,'Coastal current',fontsize=20)
ax1.text(30,-50,'Slope current',fontsize=20)
savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/uacross_mean_winst.pdf',bbox_inches='tight')

fig=figure(figsize=(16,6))
ax1=subplot(121)
meanOSM('sal','Salinity',ax1,ypos=1750,cwidth=0.015,cpos=0.49)
ax2=subplot(122)
ax2=meanOSM('tmp','Pot. Temperature',ax2,ypos=1750,cwidth=0.015,cpos=0.915)
ax2.set_yticks([])
ax2.set_yticklabels('')
ax2.set_ylabel('')
savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/saltmp_mean_winst.pdf',bbox_inches='tight')
