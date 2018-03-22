#################################################################################
#################################################################################
#################################################################################
######################## Compare Feili's product with mine #######################
#################################################################################
#################################################################################
#################################################################################
from aux_funcs import *

newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1803bathy.pickle','rb'))

[bathatmoor,newdpth]=pickle.load(open('../pickles/moordpths.pickle','rb'))

dat=newgrid.copy()
for vv in dat:
    if vv[0]!='d':
        print(vv)
        for dd,adist in enumerate(newgrid.distance):
            dat[vv][dd,:,:]=dat[vv][dd,:,:].where(dat[vv][dd,:,:].depth<=newdpth[dd])

#################################################################################
# Load Feili's data
#################################################################################
glob.glob(datadir+'Feili/*product_1803*.mat')

feili=io.loadmat(glob.glob(datadir+'Feili/*product_1803*.mat')[0])

distvec_fl=cumsum(vstack((0,sw.dist(feili['CF_gridded']['lat'][0][0],feili['CF_gridded']['lon'][0][0])[0])))+sw.dist([feili['CF_gridded']['lat'][0][0][0][0],CFlat[0]],[feili['CF_gridded']['lon'][0][0][0][0],CFlon[0]])[0][0]

feili['CF_gridded']['lat'][0][0][0][0]
sw.dist([feili['CF_gridded']['lat'][0][0][0][0],CFlat[0]],[feili['CF_gridded']['lon'][0][0][0][0],CFlon[0]])[0][0]

# figure()
# plot(CFlon,CFlat,'o-');
# plot(feili['CF_gridded']['lon'][0][0],feili['CF_gridded']['lat'][0][0],'o-');
#
# figure()
# plot(distvec,'o');
# plot(distvec_fl,'o');

def dateconv(matlab_datenum): # note, this is for integer version, acuracy is up to days
    python_datetime = datetime.datetime.fromordinal(int(matlab_datenum)) - datetime.timedelta(days=366)
    return python_datetime

fl_date=[dateconv(tt[0]) for tt in feili['CF_gridded']['time'][0][0]]

# dat_fl=xr.Dataset({'across track velocity': (['distance', 'depth', 'date'],  feili['CF_gridded']['velo'][0][0]),
#                 'salinity': (['distance', 'depth', 'date'],  feili['CF_gridded']['salt'][0][0]),
#                 'temperature': (['distance', 'depth', 'date'],  feili['CF_gridded']['ptmp'][0][0]),
#                 'potential density': (['distance', 'depth', 'date'],  feili['CF_gridded']['pden'][0][0]-1e3),},
#                 coords={'distance': distvec_fl,
#                         'depth':feili['CF_gridded']['depth'][0][0][0],
#                         'date': fl_date})

dat_fl=xr.Dataset({'across track velocity': (['distance', 'depth', 'date'],  feili['CF_gridded']['velo'][0][0]),
                'salinity': (['distance', 'depth', 'date'],  feili['CF_gridded']['psal'][0][0]),
                'temperature': (['distance', 'depth', 'date'],  feili['CF_gridded']['ptmp'][0][0]),
                'potential density': (['distance', 'depth', 'date'],  feili['CF_gridded']['pden'][0][0]-1e3),},
                coords={'distance': distvec_fl,
                        'depth':feili['CF_gridded']['depth'][0][0][0],
                        'date': fl_date})



dat_fl

#################################################################################
# Make mean sections of all fields
#################################################################################

def onecont(xchoose,field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,addcont=0):
    if axx==0:
        axx=subplot(111)
    ax1=axx.contourf(xchoose.distance,xchoose.depth,field,vrange,cmap=coloor,extend='both')
    # if 'sal' in tit:
    #     axx.contour(dat.distance,dat.depth,dat['salinity'].mean(dim='date').T,levels=[34.9],colors='w',linewidths=12)
    if addcont==1:
        ax2=contour(xchoose.distance,xchoose.depth,field,levels=hlevs,colors='k')
        clabel(ax2,fmt='%1.2f')
    axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance [km]',fontsize=14)
    ylabel('depth [m]',fontsize=14)
    xlim([-5,100])
    ylim([2200,0])
    title(tit,fontsize=22)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axx.set_yticks(range(500,2100,500))
    return ax1

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
                elif 'uac' in savename:
                    if ('AQ' in inst[key][dd]):
                        if dd==0:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='Current meter')
                        else:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6)
            mm+=1



def meanfeili(xchoose,field,tit,fig,ax1=0,ypos=1250,cwidth=0.02,cpos=0.95):
    if ax1==0:
        ax1=subplot(111)
    axvel=onecont(xchoose,xchoose[univec[field][0]].sel(date=slice('2014-8-22','2016-4-21')).mean(dim='date').T,
    field,univec[field][1],univec[field][2],univec[field][3],ax1,addcont=1,nomoorlines=0)
    cbaxes = fig.add_axes([cpos, 0.15,cwidth , 0.7])
    if field=='pden':
        cbartck=univec[field][1][::6]
    elif field=='uacross':
        cbartck=univec[field][-2]
    else:
        cbartck=univec[field][-2][::2]
    colorbar(axvel,extend='both',label=univec[field][-1],cax=cbaxes,ticks=cbartck)
    ax1.set_title(tit,fontsize=14)
    plotinstpos(ax1,field)
    if field=='uacross':
        ax1.legend(loc=(0.05,0.08),fontsize=14).set_zorder(100)
    return ax1

def plotcomp(field):
    fig=figure(figsize=(12,4))
    ax1=subplot(121)
    meanfeili(dat,field,'Isabela',fig,ax1)
    ax2=subplot(122)
    ax2=meanfeili(dat_fl,field,'Feili',fig,ax2)
    ax2.set_yticklabels('')
    ax2.set_ylabel('')
    suptitle(univec[field][0],fontsize=18)
    savefig('../figures/feili/meansection_comp_'+field+'.pdf',bbox_inches='tight')

plotcomp('sal')

plotcomp('tmp')


plotcomp('pden')

plotcomp('uacross')

#################################################################################
# Compare Transports
#################################################################################
daily=dat.copy()

mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3

ftrans=io.loadmat(glob.glob(datadir+'Feili/*Transport_1803*.mat')[0])
feili.keys()


fl_tottrans=ftrans['CF_transport']['total_d30'][0][0]
fl_uLSW=ftrans['CF_transport']['uLSW_d30'][0][0]
fl_dLSW=ftrans['CF_transport']['dLSW_d30'][0][0]


fl_tottrans_daily=ftrans['CF_transport']['total_d1'][0][0]
fl_uLSW_daily=ftrans['CF_transport']['uLSW_d1'][0][0]
fl_dLSW_daily=ftrans['CF_transport']['dLSW_d1'][0][0]

len(fl_tottrans)
len(fl_tottrans_daily)
584/30
def quickfilt(fvers):
    fout=zeros(len(fl_tottrans))
    for ii in range(len(fl_tottrans)):
        fout[ii]=mean(fvers[30*ii:30*ii+30])
    return fout

fl_tottrans_dmean=quickfilt(fl_tottrans_daily)
fl_uLSW_dmean=quickfilt(fl_uLSW_daily)
fl_dLSW_dmean=quickfilt(fl_dLSW_daily)

fl_date

def transcomp(tit,d1,d2,fvers,fvers2):
    figure(figsize=(12,3))
    Itrans=daily['xport'].sel(date=slice('2014-8-30','2016-4-20')).where((daily['potential density']<=d2)&(daily['potential density']>d1)).sum(dim='distance').sum(dim='depth').resample('30D',dim='date',base=1)
    Itrans.plot.line(marker='o',label='Isabela; '+str('{:5.2f}'.format(mean(Itrans.values)))+' Sv')
    plot(fl_date,fvers,'o-',label='Feili (from monthly prof); '+str('{:5.2f}'.format(mean(fvers)))+' Sv')
    plot(fl_date,fvers2,'o-',label='Feili (from daily prof); '+str('{:5.2f}'.format(mean(fvers2)))+' Sv')
    ylabel('Transport [Sv]')
    legend()
    title(tit)
    savefig('../figures/feili/transcomp_'+tit[:5]+'.pdf',bbox_inches='tight')


transcomp('Total transport',0,40,fl_tottrans,fl_tottrans_dmean)
transcomp('uLSW transport',27.68,27.74,fl_uLSW,fl_uLSW_dmean)
transcomp('dLSW transport',27.74,27.8,fl_dLSW,fl_dLSW_dmean)
