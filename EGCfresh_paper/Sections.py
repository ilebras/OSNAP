#################################################################################
#################################################################################
#################################################################################
######################## Make sections ###################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

newgrid=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))

savename='_1810JHcal'

dat=newgrid.copy()
dat['across track velocity']=-1*dat['across track velocity']

# bathf=interpolate.interp1d(bathdist,bathbath)
# bathonmygrid=bathf(newgrid.distance)

# dat=newgrid.copy()
# for vv in newgrid:
#     if vv[0]!='d':
#         print(vv
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
                    axchoice.plot(distvec[rr],adcpdp[rr],'^',markersize=12,label='ADCP',color='lightgreen',zorder=25,mec='k')
                else:
                    axchoice.plot(distvec[rr],adcpdp[rr],'^',markersize=12,color='lightgreen',zorder=25,mec='k')
        mm=0
        for key in depths:
            for dd in range(len(depths[key])):
                if ('sal' in savename) | ('tmp' in savename) | ('pden' in savename):
                    if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
                        if (dd==0) & (key=='CF1'):
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=4,label='Point CTD')
                        else:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=4)

                elif 'uac' in savename:
                    if ('AQ' in inst[key][dd]):
                        if dd==0:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=6,label='Current meter')
                        else:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',markersize=6)
            mm+=1


def onecont(field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,solidcont=1,addcont=1,alone=1,addlab=1,skipna=1):
    if axx==0:
        axx=subplot(111)
    if solidcont==1:
        ax1=axx.contourf(dat.distance,dat.depth[::skipna],field,vrange,cmap=coloor,extend='both')
    else:
        ax1=0
    # if 'sal' in tit:
    #     axx.contour(dat.distance,dat.depth,dat['salinity'].mean(dim='date').T,levels=[34.9],colors='w',linewidths=12)
    if addcont==1:
        ax2=axx.contour(dat.distance,dat.depth[::skipna],field,levels=hlevs,colors='k')
    if addlab==1:
        if 'temp' in tit:
            fstr='%1.0f'
        if  ('vel' in tit):
            manual_locations = [(60, 1200),(50,300),(50, 100),]
            clabels=clabel(ax2, fmt='%1.1f', fontsize=10, manual=manual_locations)
        elif 'den' in tit:
            manual_locations = [(0, 50),(52,50), (52,500),(70,1000), (70,1500),(85,2000)]
            clabels=clabel(ax2, fmt='%1.2f', fontsize=10, manual=manual_locations)
        elif 'sal' in tit:
            manual_locations = [(0, 50),(25,50), (52,100),(70,300), (88,600)]
            clabels=clabel(ax2, fmt='%1.2f', fontsize=10, manual=manual_locations)
        elif ('uac' in tit):
            manual_locations = [(0,50),(60, 1200),(60,300),(60, 100),]
            clabels=clabel(ax2, fmt='%1.1f', fontsize=10, manual=manual_locations)
        else:
            clabels=clabel(ax2,fmt='%1.2f')
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

def saldenline(d1,d2,axit,tit):
    field='sal'
    ax1,ax2=onecont(dat[univec[field][0]].sel(date=slice(d1,d2)).mean(dim='date').T,'',univec[field][1],univec[field][2],univec[field][3],axx=axit,alone=0,nomoorlines=1,addcont=0,addlab=0)
    # axit.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(d1,d2)).mean(dim='date').T,[34.9],colors='limegreen',linewidth=2)
    axit.contour(dat.distance,dat.depth,dat['potential density'].sel(date=slice(d1,d2)).mean(dim='date').T,arange(25,28,0.1),colors='k')
    axit.axvline(40,color='purple',linewidth=4)
    axit.axvline(60,color='purple',linewidth=4)
    axit.text(-10,400,tit,color='white',fontsize=16,zorder=100)
    return ax1

def plotsaldenseas():
    f,axx=subplots(1,3,sharex=True,sharey=True,figsize=(11,2))
    salcont=saldenline('2014-10-01','2015-1-1',axx[0],'Fall')
    salcont=saldenline('2015-1-01','2015-4-1',axx[1],'Winter')
    salcont=saldenline('2015-4-01','2015-10-1',axx[2],'Spring/ \nSummer')
    axx[0].set_ylim(500,50)
    axx[0].set_xlim(-15,95)
    axx[0].set_yticks(arange(50,500,200))
    axx[0].set_ylabel('depth [m]')
    axx[1].set_xlabel('distance [km]')
    cbax=f.add_axes([0.925,0.15,0.015,0.7])

    colorbar(salcont,cax=cbax,label='Salinity',ticks=(34,34.9,34.96,35))
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/SalDenSecsSeas.pdf',bbox_inches='tight')


plotsaldenseas()


CTD=pickle.load(open(datadir+'OSNAP2016recovery/pickles/Shipboard/CTD_xarray.pickle','rb'))


def schemplot_CTDtest():
    f,ax1=subplots(2,1,figsize=(4,3),sharex=True,gridspec_kw={'height_ratios':[1,4]})
    var='sal'
    occ=-1
    for axx in ax1:
        ret=axx.contourf(hstack((CTD['distance [km]'][:3],CTD['distance [km]'][4:]))-10,CTD['pressure [db]'],hstack((CTD.salinity[:,:3,-1],CTD.salinity[:,4:,-1])),levels=univec[var][1],cmap=univec[var][2],extend='both')
        axx.contour(hstack((CTD['distance [km]'][:3],CTD['distance [km]'][4:]))-10,CTD['pressure [db]'],hstack((CTD.density[:,:3,-1],CTD.density[:,4:,-1])),levels=arange(25,28,0.1),colors='k',linewidth=2)
        axx.fill_between(bathdistnew,bathbathnew,2500*ones(len(bathbathnew)),color='k',zorder=22)
        axx.axvline(40,color='purple',linewidth=4)
        axx.axvline(60,color='purple',linewidth=4)
    ax1[0].set_ylim(50,20)
    ax1[1].set_ylim(500,50)
    ax1[1].set_ylabel('depth [m]')
    ax1[1].set_xlabel('distance [km]')

    xlim(5,120)
    cbax=f.add_axes([0.95,0.15,0.03,0.7])
    colorbar(ret,cax=cbax,label='Salinity',ticks=(34,34.9,34.96,35))
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/SalDenSecShipboard.pdf',bbox_inches='tight')

schemplot_CTDtest()


#####################################################################
# Make 4 panel plot of all means for paper with instrument positions
#################################################################################




def eachpan(field,axit,atit):
    ax1,ax2=onecont(dat[univec[field][0]].mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3],axx=axit,alone=0,nomoorlines=0)
    plotinstpos(axit,field)
    if ('sal' in field):
        cticks=[34, 34.6,34.9, 34.94,34.96, 35]
    elif 'uac' in field:
        cticks=arange(0,0.601,0.1)
    elif ('tmp' in field):
        cticks=range(0,10,2)
    else:
        cticks=univec[field][3]
    colorbar(ax1,ticks=cticks,ax=axit)
    axit.text(-10,2100,atit,color='white',zorder=1e3,fontsize=16)


colors = [(12,44,132),(158,202,225) ,(255,237,160),(217,95,14),(240,59,32)]#(237,248,177),,
sal_cmap = make_cmap(colors,position=[0,0.9,0.96,0.99,1],bit=True)#0.9,

# univec['uacross'][2]=cm.BuPu
# univec['uacross'][1]=arange(-0.6,0.001,0.025)
univec['sal']=['salinity',array([32, 34, 34.6, 34.9,34.94,34.96, 35,35.02]),sal_cmap,array([33,34, 34.6,34.9,34.94,34.96,34.98,35]),'']
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



def proposal_sections():
    f, (ax11, ax22) = plt.subplots(1, 2, sharex=True, sharey=True,figsize=(12,4))
    eachpan('uacross',ax11,'')
    fss=13
    ax11.legend(loc=(0.05,0.03),fontsize=fss).set_zorder(100)
    eachpan('sal',ax22,'')
    ax22.legend(loc=(0.05,0.03),fontsize=fss).set_zorder(100)
    f.text(0, 1, 'a) Mean mooring downstream velocity [m/s]',fontsize=16)
    f.text(0.5, 1, 'b) Mean mooring salinity (2014-2016)',fontsize=16)
    # f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
    # f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
    ax11.set_ylabel('depth [m]')
    ax11.set_xlabel('distance [km]')
    ax22.set_xlabel('distance [km]')
    [ax11.text(distvec[ii]-4,-100,'CF'+str(ii+1),fontsize=fss) for ii in range(7)]
    ax11.text(distvec[-1]-4,-100,'M1',fontsize=fss)
    [ax22.text(distvec[ii]-4,-100,'CF'+str(ii+1),fontsize=fss) for ii in range(7)]
    ax22.text(distvec[-1]-4,-100,'M1',fontsize=fss)
    ax11.plot(distvec[:-1],[210,210,210,410,1350,1850,1900],'ro',markersize=10,zorder=100)
    ax11.plot(distvec[:-1],[210,210,210,410,1350,1850,1900],'ko',markersize=4,zorder=101)
    ax22.plot(distvec[3],40,'ro',markersize=10)
    ax22.plot(distvec[3],40,'ko',markersize=4)
    plt.tight_layout()
    savefig('/home/isabela/Documents/proposals/OSNAP_mixing/OSNAPmixing-proposal/figures/Moorsecs.pdf',bbox_inches='tight')

proposal_sections()



def mkbox(x1,x2,y1,y2,axx):
    axx.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'k-',linewidth=5)

def partmean(date1,date2,textit,axx,field,textit2):
    d1tot='2014-10-01'
    d2tot='2015-10-1'
    totmean=dat[univec[field][0]].sel(date=slice(d1tot,d2tot)).mean(dim='date').T
    if 'uac' in field:
        meanie=mean(totmean).values
        axvel,na=onecont((dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T-totmean)[::10,:],field,arange(-0.2,0.25,0.05),
        cm.RdBu_r,arange(-0.2,0.25,0.05),axx,alone=0,addlab=0,skipna=10)
    elif 'sal' in field:
        axvel,na=onecont(dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T-totmean,field,arange(-0.25,0.3,0.05),
        cm.RdBu_r,arange(-0.25,0.3,0.05),axx,alone=0,addlab=0)
    elif 'tmp' in field:
        axvel,na=onecont(dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T-totmean,field,arange(-2,2.2,0.5),
        cm.RdBu_r,arange(-2,2.1,0.5),axx,alone=0,addlab=0)
    if 'sal' in field:
    #     # mkbox(45,62,5,500,axx)
        ax349=axx.contour(dat.distance,dat.depth.isel(depth=slice(0,300)),
                    dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').isel(depth=slice(0,300)).T,
                    levels=[34.9],colors='w',linewidths=6)
        ax349=axx.contour(dat.distance,dat.depth.isel(depth=slice(0,300)),
                        dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').isel(depth=slice(0,300)).T,
                        levels=[34.9],colors='purple',linewidths=3)
        # clabels=clabel(ax349, fmt='%1.1f', )
    # # if 'uac' in field:
    # #     ax349=axx.contour(dat.distance,dat.depth.isel(depth=slice(0,300)),dat['across track velocity'].sel(date=slice(date1,date2)).mean(dim='date').isel(depth=slice(0,300)).T,
    # #     levels=[-0.4],colors='k',linewidths=4)
    axx.text(-100,350,textit,color='k',fontsize=25,zorder=200)
    # axx.text(-10,425,textit2,color='w',fontsize=25,zorder=200)
    axx.set_ylim(500,0)
    axx.set_xlim(-15,95)
    axx.set_yticks([100,300,500])

    return axvel

def plotvel_sal_seas_wmean_wtmp():
    f, axx = plt.subplots(4,3, sharex=True, sharey=True,figsize=(9,9))
    field='uacross'
    axvmn1,axvmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],univec[field][3],univec[field][2],univec[field][3],axx=axx[0,0],alone=0)
    axx[0,0].text(-100,350,'Mean',color='k',fontsize=25,zorder=200)
    # axx[0,0].text(-10,425,'a)',color='w',fontsize=25,zorder=200)

    partmean('2014-10-01','2015-1-1','Fall',axx[1,0],'uacross','d)')
    partmean('2015-01-01','2015-04-1','Winter',axx[2,0],'uacross','g)')
    axvel=partmean('2015-04-01','2015-10-1','Spring/\nSummer',axx[3,0],'uacross','j)')

    field='sal'
    axsmn1,axsmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                            arange(32.5,35.1,0.5),cm.BuPu_r,sort(hstack((arange(32.5,35.1,0.5),34.85,34.95))),axx=axx[0,1],alone=0)
    # axx[0,1].text(-10,425,'b)',color='w',fontsize=25,zorder=200)
    fs=14
    partmean('2014-10-01','2015-1-1','',axx[1,1],'sal','e)')
    sal1=axx[1,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal1.set_backgroundcolor('white')
    # sal1.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    partmean('2015-01-01','2015-04-1','',axx[2,1],'sal','h)')
    sal2=axx[2,1].text(50,450,'34.9',color='purple',fontsize=fs)
    # sal2.set_backgroundcolor('white')
    # sal2.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    axsal=partmean('2015-04-01','2015-10-1','',axx[3,1],'sal','k)')
    sal3=axx[3,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal3.set_backgroundcolor('white')
    # sal3.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    field='tmp'
    axtmn1,axtmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                          range(-1,6,1),cm.BuPu,range(-1,8),axx=axx[0,2],alone=0)
    # axx[0,2].text(-10,425,'c)',color='w',fontsize=25,zorder=200)
    partmean('2014-10-01','2015-1-1','',axx[1,2],'tmp','c)')
    partmean('2015-01-01','2015-04-1','',axx[2,2],'tmp','i)')
    axtmp=partmean('2015-04-01','2015-10-1','',axx[3,2],'tmp','l)')

    fclab=18
    f.tight_layout()
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=20)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=20)
    f.text(0.11,1.06,'Velocity [m/s]',fontsize=fclab)
    f.text(0.03,-0.15,'Velocity anomaly [m/s]',fontsize=fclab)
    cwid=0.26
    cbaxes1 = f.add_axes([0.07, 1.02, cwid, 0.02])
    cbaxes1b = f.add_axes([0.07, -0.08,cwid, 0.02])
    f.colorbar(axvmn1,ticks=arange(0,0.7,0.2),orientation='horizontal',cax=cbaxes1)
    f.colorbar(axvel,ticks=arange(-0.2,0.2,0.1),orientation='horizontal',cax=cbaxes1b)
    f.text(0.47,1.06,'Salinity',fontsize=fclab)
    f.text(0.4,-0.15,'Salinity anomaly',fontsize=fclab)
    cbaxes2 = f.add_axes([0.39, 1.02, cwid, 0.02])
    cbaxes2b = f.add_axes([0.39, -0.08, cwid, 0.02])
    f.colorbar(axsmn1,ticks=arange(32,35.1,),orientation='horizontal',cax=cbaxes2)
    f.colorbar(axsal,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes2b)
    f.text(0.71,1.06,'Temperature [$^\circ$C]',fontsize=fclab)
    f.text(0.66,-0.15,'Temperature anomaly [$^\circ$C]',fontsize=fclab)
    cbaxes3 = f.add_axes([0.71, 1.02, cwid, 0.02])
    cbaxes3b = f.add_axes([0.71, -0.08, cwid, 0.02])
    f.colorbar(axtmn1,ticks=arange(0,8,2),orientation='horizontal',cax=cbaxes3)
    f.colorbar(axtmp,ticks=arange(-2,2.1,1),orientation='horizontal',cax=cbaxes3b)

    savefig('../figures/paperfigs/vel_sal_seasonal_wmean_wtmp.pdf',bbox_inches='tight')

plotvel_sal_seas_wmean_wtmp()


def plotvel_sal_seas_nospring():
    f, axx = plt.subplots(3,3, sharex=True, sharey=True,figsize=(9,7))
    field='uacross'
    axvmn1,axvmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],univec[field][3],univec[field][2],univec[field][3],axx=axx[0,0],alone=0)
    axx[0,0].text(-100,350,'Mean',color='k',fontsize=25,zorder=200)
    # axx[0,0].text(-10,425,'a)',color='w',fontsize=25,zorder=200)
    partmean('2014-10-01','2015-1-1','Fall',axx[1,0],'uacross','d)')
    axvel=partmean('2015-01-01','2015-04-1','Winter',axx[2,0],'uacross','g)')
    # axvel=partmean('2015-04-01','2015-10-1','Spring/\nSummer',axx[3,0],'uacross','j)')

    field='sal'
    axsmn1,axsmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                            arange(32.5,35.1,0.5),cm.BuPu_r,hstack((arange(32.5,35.1,0.5),34.85,34.95)),axx=axx[0,1],alone=0)
    # axx[0,1].text(-10,425,'b)',color='w',fontsize=25,zorder=200)
    fs=14
    partmean('2014-10-01','2015-1-1','',axx[1,1],'sal','e)')
    sal1=axx[1,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal1.set_backgroundcolor('white')
    # sal1.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    axsal=partmean('2015-01-01','2015-04-1','',axx[2,1],'sal','h)')
    sal2=axx[2,1].text(50,450,'34.9',color='purple',fontsize=fs)
    # sal2.set_backgroundcolor('white')
    # sal2.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    # axsal=partmean('2015-04-01','2015-10-1','',axx[3,1],'sal','k)')
    # sal3=axx[3,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # # sal3.set_backgroundcolor('white')
    # # sal3.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    field='tmp'
    axtmn1,axtmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                          range(-1,6,1),cm.BuPu,range(-1,8),axx=axx[0,2],alone=0)
    # axx[0,2].text(-10,425,'c)',color='w',fontsize=25,zorder=200)
    partmean('2014-10-01','2015-1-1','',axx[1,2],'tmp','c)')
    axtmp=partmean('2015-01-01','2015-04-1','',axx[2,2],'tmp','i)')
    # axtmp=partmean('2015-04-01','2015-10-1','',axx[3,2],'tmp','l)')

    fclab=18
    f.tight_layout()
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=20)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=20)
    f.text(0.11,1.06,'Velocity [m/s]',fontsize=fclab)
    f.text(0.03,-0.15,'Velocity anomaly [m/s]',fontsize=fclab)
    cwid=0.26
    cbaxes1 = f.add_axes([0.07, 1.02, cwid, 0.02])
    cbaxes1b = f.add_axes([0.07, -0.08,cwid, 0.02])
    f.colorbar(axvmn1,ticks=arange(0,0.7,0.2),orientation='horizontal',cax=cbaxes1)
    f.colorbar(axvel,ticks=arange(-0.2,0.2,0.1),orientation='horizontal',cax=cbaxes1b)
    f.text(0.47,1.06,'Salinity',fontsize=fclab)
    f.text(0.4,-0.15,'Salinity anomaly',fontsize=fclab)
    cbaxes2 = f.add_axes([0.39, 1.02, cwid, 0.02])
    cbaxes2b = f.add_axes([0.39, -0.08, cwid, 0.02])
    f.colorbar(axsmn1,ticks=arange(32,35.1,),orientation='horizontal',cax=cbaxes2)
    f.colorbar(axsal,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes2b)
    f.text(0.71,1.06,'Temperature [$^\circ$C]',fontsize=fclab)
    f.text(0.66,-0.15,'Temperature anomaly [$^\circ$C]',fontsize=fclab)
    cbaxes3 = f.add_axes([0.71, 1.02, cwid, 0.02])
    cbaxes3b = f.add_axes([0.71, -0.08, cwid, 0.02])
    f.colorbar(axtmn1,ticks=arange(0,8,2),orientation='horizontal',cax=cbaxes3)
    f.colorbar(axtmp,ticks=arange(-2,2.1,1),orientation='horizontal',cax=cbaxes3b)

    savefig('/home/isabela/Documents/conferences/1810_OSU/presentation/figures/vel_sal_seasonal_nospring.pdf',bbox_inches='tight')

plotvel_sal_seas_nospring()



XXXXXXXXXXXXXX

# univec['uacross'][1]=arange(-0.7,0.1,0.025)
# d1tot='2014-10-01'
# d2tot='2015-9-1'
# field='sal'
# totmean=dat[univec[field][0]].sel(date=slice(d1tot,d2tot)).mean(dim='date').T

# meanie=mean(totmean).values
# meanie
# salcont=array([32, 34, 34.6,34.9,34.92,34.94,34.96, 35,35.02])
### commented from here 180806
# salcont=array([32,33,34,34.3,34.6,34.9,34.92,34.94,34.96, 35,35.02])-34.9
# salcont
#
# univec['sal'][1]=sort(hstack((salcont,
#                               salcont[:-1]+diff(salcont)/8,
#                               salcont[:-1]+diff(salcont)/4,
#                               salcont[:-1]+diff(salcont)*3/8,
#                               salcont[:-1]+diff(salcont)/2,
#                               salcont[:-1]+diff(salcont)*5/8,
#                               salcont[:-1]+diff(salcont)*3/4,
#                               salcont[:-1]+diff(salcont)*7/8,
#                               )))
# univec['sal'][3]=array([33,34,34.6,34.8,34.9,34.92,34.94,34.96,34.98,35])-34.9
#
# colors = [(5,48,97),(67,147,195),(247,247,247),(214,96,77) ,(103,0,31)]#(237,248,177),,
# sal_cmap = make_cmap(colors,position=[0,0.9,0.96,0.99,1],bit=True)

# def plotvel_sal_seas():
#     f, axx = plt.subplots(3,2, sharex=True, sharey=True,figsize=(8,6.5))
#     partmean('2014-10-01','2015-1-1','Fall',axx[0,0],'uacross')
#     partmean('2015-01-01','2015-04-1','Winter',axx[1,0],'uacross')
#     axvel=partmean('2015-04-01','2015-10-1','Spring/\nSummer',axx[2,0],'uacross')
#
#     partmean('2014-10-01','2015-1-1','',axx[0,1],'sal')
#     partmean('2015-01-01','2015-04-1','',axx[1,1],'sal')
#     axsal=partmean('2015-04-01','2015-10-1','',axx[2,1],'sal')
#     f.tight_layout()
#     f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
#     f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
#     f.text(0.1,1.06,'Velocity anomaly [m/s]',fontsize=20)
#     cbaxes1 = f.add_axes([0.07, 1.02, 0.42, 0.02])
#     f.colorbar(axvel,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes1)
#     f.text(0.6,1.06,'Salinity anomaly',fontsize=20)
#     cbaxes2 = f.add_axes([0.54, 1.02, 0.42, 0.02])
#     f.colorbar(axsal,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes2)
#
#     savefig('../figures/paperfigs/vel_sal_seasonal.pdf',bbox_inches='tight')
#
# plotvel_sal_seas()

def plotvel_sal_seas_wmean():
    f, axx = plt.subplots(4,2, sharex=True, sharey=True,figsize=(8,10))
    field='uacross'
    axvmn1,axvmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3],axx=axx[0,0],alone=0)
    axx[0,0].text(-82,350,'Mean',color='k',fontsize=25,zorder=200)
    axx[0,0].text(-10,425,'a)',color='w',fontsize=25,zorder=200)

    partmean('2014-10-01','2015-1-1','Fall',axx[1,0],'uacross','c)')
    partmean('2015-01-01','2015-04-1','Winter',axx[2,0],'uacross','e)')
    axvel=partmean('2015-04-01','2015-10-1','Spring/\nSummer',axx[3,0],'uacross','g)')

    field='sal'
    axsmn1,axsmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                            arange(32.5,35.1,0.05),cm.BuPu_r,hstack((arange(32.5,35.1,0.5),34.85,34.95)),axx=axx[0,1],alone=0)
    axx[0,1].text(-10,425,'b)',color='w',fontsize=25,zorder=200)
    fs=14
    partmean('2014-10-01','2015-1-1','',axx[1,1],'sal','d)')
    sal1=axx[1,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal1.set_backgroundcolor('white')
    # sal1.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    partmean('2015-01-01','2015-04-1','',axx[2,1],'sal','f)')
    sal2=axx[2,1].text(50,450,'34.9',color='purple',fontsize=fs)
    # sal2.set_backgroundcolor('white')
    # sal2.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    axsal=partmean('2015-04-01','2015-10-1','',axx[3,1],'sal','h)')
    sal3=axx[3,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal3.set_backgroundcolor('white')
    # sal3.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    f.tight_layout()
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
    f.text(0.15,1.06,'Velocity [m/s]',fontsize=20)
    f.text(0.07,-0.15,'Velocity anomaly [m/s]',fontsize=20)
    cbaxes1 = f.add_axes([0.07, 1.02, 0.42, 0.02])
    cbaxes1b = f.add_axes([0.07, -0.08, 0.42, 0.02])
    f.colorbar(axvmn1,ticks=arange(-0.6,0.1,0.2),orientation='horizontal',cax=cbaxes1)
    f.colorbar(axvel,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes1b)
    f.text(0.7,1.06,'Salinity',fontsize=20)
    f.text(0.6,-0.15,'Salinity anomaly',fontsize=20)
    cbaxes2 = f.add_axes([0.54, 1.02, 0.42, 0.02])
    cbaxes2b = f.add_axes([0.54, -0.08, 0.42, 0.02])
    f.colorbar(axsmn1,ticks=arange(32.5,35.1,0.5),orientation='horizontal',cax=cbaxes2)
    f.colorbar(axsal,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes2b)

    savefig('../figures/paperfigs/vel_sal_seasonal_wmean.pdf',bbox_inches='tight')

plotvel_sal_seas_wmean()

univec['tmp'][1]
def plotvel_sal_seas_wmean_wtmp():
    f, axx = plt.subplots(4,3, sharex=True, sharey=True,figsize=(8,10))
    field='uacross'
    axvmn1,axvmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3],axx=axx[0,0],alone=0)
    axx[0,0].text(-82,350,'Mean',color='k',fontsize=25,zorder=200)
    axx[0,0].text(-10,425,'a)',color='w',fontsize=25,zorder=200)

    partmean('2014-10-01','2015-1-1','Fall',axx[1,0],'uacross','c)')
    partmean('2015-01-01','2015-04-1','Winter',axx[2,0],'uacross','e)')
    axvel=partmean('2015-04-01','2015-10-1','Spring/\nSummer',axx[3,0],'uacross','g)')

    field='sal'
    axsmn1,axsmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                            arange(32.5,35.1,0.05),cm.BuPu_r,hstack((arange(32.5,35.1,0.5),34.85,34.95)),axx=axx[0,1],alone=0)
    axx[0,1].text(-10,425,'b)',color='w',fontsize=25,zorder=200)
    fs=14
    partmean('2014-10-01','2015-1-1','',axx[1,1],'sal','d)')
    sal1=axx[1,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal1.set_backgroundcolor('white')
    # sal1.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    partmean('2015-01-01','2015-04-1','',axx[2,1],'sal','f)')
    sal2=axx[2,1].text(50,450,'34.9',color='purple',fontsize=fs)
    # sal2.set_backgroundcolor('white')
    # sal2.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    axsal=partmean('2015-04-01','2015-10-1','',axx[3,1],'sal','h)')
    sal3=axx[3,1].text(50,250,'34.9',color='purple',fontsize=fs)
    # sal3.set_backgroundcolor('white')
    # sal3.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    field='tmp'
    axsmn1,axsmn2=onecont(dat[univec[field][0]].sel(date=slice('2014-10-01','2015-10-01')).mean(dim='date').T,'Mean '+univec[field][0],
                          univec[field][1],cm.BuPu_r,univec[field][3],axx=axx[0,2],alone=0)

    partmean('2014-10-01','2015-1-1','',axx[1,2],'tmp','c)')
    partmean('2015-01-01','2015-04-1','',axx[2,2],'tmp','e)')
    axvel=partmean('2015-04-01','2015-10-1','',axx[3,2],'tmp','g)')


    f.tight_layout()
    f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
    f.text(0.15,1.06,'Velocity [m/s]',fontsize=20)
    f.text(0.07,-0.15,'Velocity anomaly [m/s]',fontsize=20)
    cbaxes1 = f.add_axes([0.07, 1.02, 0.42, 0.02])
    cbaxes1b = f.add_axes([0.07, -0.08, 0.42, 0.02])
    f.colorbar(axvmn1,ticks=arange(0,0.7,0.2),orientation='horizontal',cax=cbaxes1)
    f.colorbar(axvel,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes1b)
    f.text(0.7,1.06,'Salinity',fontsize=20)
    f.text(0.6,-0.15,'Salinity anomaly',fontsize=20)
    cbaxes2 = f.add_axes([0.54, 1.02, 0.42, 0.02])
    cbaxes2b = f.add_axes([0.54, -0.08, 0.42, 0.02])
    f.colorbar(axsmn1,ticks=arange(32.5,35.1,0.5),orientation='horizontal',cax=cbaxes2)
    f.colorbar(axsal,ticks=arange(-0.2,0.21,0.1),orientation='horizontal',cax=cbaxes2b)

    savefig('../figures/paperfigs/vel_sal_seasonal_wmean_wtmp.pdf',bbox_inches='tight')

plotvel_sal_seas_wmean_wtmp()

XXXXXXXXXXXXXX
## Residual overlay scripting
    # axden=axx.contour(dat.distance,dat.depth,dat['potential density'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[27.4,27.5,27.6],colors='g',linewidths=2)
    # ax34=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.2],colors='orange',linewidths=2)
    # ax349=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.9],colors='k',linewidths=4)
    # clab34=clabel(ax34,fmt='%1.1f',manual=[(25,100)],colors='k')
    # clab349=clabel(ax349,fmt='%1.1f',manual=[(50,250)])
    # if 'WINTER' not in textit:
    #     ax3495=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.95],colors='darkred',linewidths=2)
    #     clab3495=clabel(ax3495,fmt='%1.2f',manual=[(80,500)],colors='k')
    #     clablist=[clab34,clab349,clab3495]
    # else:
    #     clablist=[clab34,clab349]
    # for clabels in clablist:
    #     [txt.set_backgroundcolor('white') for txt in clabels]
    #     [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0,alpha=0.8)) for txt in clabels]

#
#
#
# def partmean_pden(date1,date2,textit,axx):
#     field='pden'
#     axvel,na=onecont(dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T,'vel',univec[field][1],
#     univec[field][2],univec[field][3],axx,addcont=0,alone=0)
#     ax34=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.2],colors='orange',linewidths=2)
#     ax349=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.9],colors='k',linewidths=4)
#
#     clab34=clabel(ax34,fmt='%1.1f',manual=[(25,100)],colors='k')
#     clab349=clabel(ax349,fmt='%1.1f',manual=[(50,250)])
#     if 'WINTER' not in textit:
#         ax3495=axx.contour(dat.distance,dat.depth,dat['salinity'].sel(date=slice(date1,date2)).mean(dim='date').T,levels=[34.95],colors='darkred',linewidths=2)
#         clab3495=clabel(ax3495,fmt='%1.2f',manual=[(80,500)],colors='k')
#         clablist=[clab34,clab349,clab3495]
#     else:
#         clablist=[clab34,clab349]
#     for clabels in clablist:
#         [txt.set_backgroundcolor('white') for txt in clabels]
#         [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0,alpha=0.8)) for txt in clabels]
#     axx.text(-10,2100,textit,color='white',fontsize=20,zorder=200)
#     axx.set_ylim(1000,0)
#
#     return axvel
#
# def plot4pden_seas_symm():
#     f, (ax11, ax22, ax33) = plt.subplots(3,1, sharex=True, sharey=True,figsize=(5,8))
#     partmean_pden('2014-10-01','2015-1-1','FALL',ax11)
#     partmean_pden('2015-01-01','2015-04-1','WINTER',ax22)
#     axvel=partmean_pden('2015-04-01','2015-9-1','SPRING/SUMMER',ax33)
#     # axvel=partmean('2015-7-01','2015-10-1','SUMMER',ax44)
#     f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
#     f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
#     plt.tight_layout()
#     cbaxes = f.add_axes([1.0, 0.3, 0.06, 0.45])
#     cbar=colorbar(axvel,ticks=univec['pden'][3],label='Potential density [kg/m$^3$]',cax=cbaxes)
#     savefig('../figures/paperfigs/pden_seasonal_symm.pdf',bbox_inches='tight')
#
#
# plot4pden_seas_symm()


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
univec['pden']

def partmeanQQ(date1,date2,tit,axx,field,field2):
    onecont(dat[univec[field][0]].sel(date=slice(date1,date2)).mean(dim='date').T,'vel',univec[field][1],
    univec[field][2],univec[field][3],axx,alone=0)
    onecont(dat[univec[field2][0]].sel(date=slice(date1,date2)).mean(dim='date').T,'vel',univec[field][1],
    univec[field][2],univec[field][3],axx,alone=0,solidcont=0)
    axx.text(-10,2100,tit,color='white',fontsize=20,zorder=200)

def seasplot(field,field2):
        f, (ax11, ax22, ax33) = plt.subplots(3,1, sharex=True, sharey=True,figsize=(5,10))
        partmeanQQ('2014-10-01','2015-1-1','FALL',ax11,field,field2)
        partmeanQQ('2015-01-01','2015-04-1','WINTER',ax22,field,field2)
        axvel=partmeanQQ('2015-04-01','2015-9-1','SPRING/SUMMER',ax33,field,field2)

seasplot('pden','sal')
# savefig('../figures/sections/pden_seas.pdf',bbox_inches='tight')
