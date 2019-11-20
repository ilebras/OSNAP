#################################################################################
#################################################################################
#################################################################################
######################## Make sections ###################################
#################################################################################
#################################################################################
#################################################################################

from firstfuncs_1618 import *

daily=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid.nc')
daily['across track velocity']=-1*daily['across track velocity']
daily=daily.where(daily['across track velocity']!=0)
daily=daily.sel(date=slice('2014-8-15','2018-9-1'))

dat=daily.copy()

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
    [axit.axvline(ss,color='purple',linewidth=3) for ss in distvec]
    # axit.axvline(60,color='purple',linewidth=4)
    axit.text(-10,400,tit,color='white',fontsize=16,zorder=100)
    return ax1


def plotsaldenseas():
    for ii in range(1,5):
        f,axx=subplots(1,3,sharex=True,sharey=True,figsize=(11,3))
        salcont=saldenline(str(2013+ii)+'-10-01',str(2014+ii)+'-1-1',axx[0],'Fall')
        salcont=saldenline(str(2014+ii)+'-1-01',str(2014+ii)+'-4-1',axx[1],'Winter')
        salcont=saldenline(str(2014+ii)+'-4-01',str(2014+ii)+'-10-1',axx[2],'Spring/ \nSummer')
        axx[0].set_ylim(1000,50)
        axx[0].set_xlim(-15,95)
        axx[0].set_yticks(arange(50,1000,300))
        axx[0].set_ylabel('depth [m]')
        axx[1].set_xlabel('distance [km]')
        cbax=f.add_axes([0.925,0.15,0.015,0.7])

        colorbar(salcont,cax=cbax,label='Salinity',ticks=(34,34.9,34.96,35))
        savefig(figdir+'ReAnalyze/egc_fresh/Saldensec_year'+str(ii)+'.png',bbox_inches='tight')


plotsaldenseas()


colors = [(12,44,132),(158,202,225) ,(255,237,160),(217,95,14),(240,59,32)]#(237,248,177),,
sal_cmap = make_cmap(colors,position=[0,0.9,0.96,0.99,1],bit=True)#0.9,
adcpdp=[160,160,160,340,90,90,90]



depths={}
for mm in range(1,8):
    depths['CF'+str(mm)]=ones(len(info['CF'+str(mm)+'_sensors'][0]))
    for ii in range(len(info['CF'+str(mm)+'_sensors'][0])):
        if info['CF'+str(mm)+'_sensors'][0][ii][1][0][0]=='B':
            depths['CF'+str(mm)][ii]=info['CF'+str(mm)+'_header']['depth_actual'][0][0][0][0]
        else:
            depths['CF'+str(mm)][ii]=array(info['CF'+str(mm)+'_sensors'][0][ii][1][0][0])



inst={}
for mm in range(1,8):
    inst['CF'+str(mm)]=[info['CF'+str(mm)+'_sensors'][0][ii][2][0] for ii in range(len(info['CF'+str(mm)+'_sensors'][0]))]

inst['CF5'][0]='lost'
depths['M1']=array([50,50,290,530,530,765,1000,1000,1250,1500,
                    1500,1730,1730,1950,1950])
inst['M1']=array(['AQ','MC','MC','AQ','MC','MC','AQ','MC','MC','AQ',
                  'MC','AQ','MC','MC','WH'])


univec['sal']=['salinity',array([32, 34, 34.6, 34.9,34.94,34.96, 35,35.02]),sal_cmap,array([33,34, 34.6,34.9,34.94,34.96,34.98,35]),'']




def eachpan(field,axit,atit,d1='2014-9-1',d2='2018-9-1'):
    ax1,ax2=onecont(dat[univec[field][0]].sel(date=slice(d1,d2)).mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3],axx=axit,alone=0,nomoorlines=0)
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
    savefig(figdir+'ReAnalyze/egc_fresh/MeanSections_1418.png',bbox_inches='tight')

plot4pan_means()

def plot4pan_means_yearly():
    for ii in range(1,5):
        f, ((ax11, ax22), (ax33, ax44)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(11,7))
        eachpan('uacross',ax11,'a) Velocity [m/s]',d1=str(2013+ii)+'-9-1',d2=str(2014+ii)+'-9-1')
        eachpan('pden',ax22,'b) Potential density [kg/m$^3$]',d1=str(2013+ii)+'-9-1',d2=str(2014+ii)+'-9-1')
        eachpan('sal',ax33,'c) Salinity',d1=str(2013+ii)+'-9-1',d2=str(2014+ii)+'-9-1')
        eachpan('tmp',ax44,'d) Potential temperature [$^\circ$C]',d1=str(2013+ii)+'-9-1',d2=str(2014+ii)+'-9-1')
        f.text(0.5, -0.03, 'distance [km]', ha='center',fontsize=16)
        f.text(-0.03, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=16)
        ax11.text(-20,-100,'Coastal current',fontsize=14)
        ax11.text(35,-100,'Slope current',fontsize=14)
        [ax22.text(distvec[ii]-4,-100,'CF'+str(ii+1),fontsize=13) for ii in range(7)]
        ax22.text(distvec[-1]-4,-100,'M1',fontsize=13)
        plt.tight_layout()
        savefig(figdir+'ReAnalyze/egc_fresh/MeanSections_1418_year'+str(ii)+'.png',bbox_inches='tight')

plot4pan_means_yearly()
