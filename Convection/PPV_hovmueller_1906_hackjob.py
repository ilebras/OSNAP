
###############################################################################################################################
###############################################################################################################################
############################################## PV timeseries with MLDs   ####################################################
###############################################################################################################################
###############################################################################################################################
from aux_funcs import *

from map_funcs import *

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))
ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb'))
oom=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_denmerged_xray.pickle','rb'))

months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')

d1=27.65
d2=27.73
d3=27.77

oom_ml=io.loadmat(datadir+'IrmingerSea/OOImf_MLD_PD_orig.mat')['MLD_PD'][0]
oom_time=io.loadmat(datadir+'IrmingerSea/OOImf_mtime_orig.mat')['time'][0]
oom_date=date_from_matlab(oom_time)

ooi_ML=xr.Dataset({'ML':('date',oom_ml)},coords={'date':oom_date})


ooi_ML_weekly=ooi_ML.resample(date='1W').max(dim='date')
ooi_ML_weekly=ooi_ML_weekly.fillna(75)

ML_weekly=dat.ML.resample(date='1W').max(dim='date')
ML_weekly=ML_weekly.fillna(75)

#### smooth out the PV and density fields

dat['PV_sm']=NaN*dat['PV']
dat['pden_sm']=NaN*dat['PV']

for tt,na in enumerate(dat.date.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['PV'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.013, output='ba')
            dat['PV_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['PV'][mm,nanind,tt].values)
            # dat['PV_sm'][mm,nanind,tt]=dat['PV'][mm,nanind,tt].values
        else:
            dat['PV_sm'][mm,:,tt]=NaN
        nanind=~isnan(dat['potential density'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.01, output='ba')
            dat['pden_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['potential density'][mm,nanind,tt].values)
        else:
            dat['pden_sm'][mm,:,tt]=NaN

for dd,na in enumerate(dat.depth.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['PV_sm'][mm,dd,:])
        if sum(nanind)>10:
            B, A = sig.butter(2,0.1, output='ba')
            dat['PV_sm'][mm,dd,nanind]=sig.filtfilt(B,A,dat['PV_sm'][mm,dd,nanind].values)
            # dat['PV_sm'][mm,dd,nanind]=dat['PV_sm'][mm,dd,nanind].values
        else:
            dat['PV_sm'][mm,dd,:]=NaN
        nanind=~isnan(dat['pden_sm'][mm,dd,:])
        if sum(nanind)>10:
            B, A = sig.butter(2,0.05, output='ba')
            dat['pden_sm'][mm,dd,nanind]=sig.filtfilt(B,A,dat['pden_sm'][mm,dd,nanind].values)
        else:
            dat['pden_sm'][mm,dd,:]=NaN

##smooth the OOI products too
ooi['PV_sm']=NaN*ooi['PV']
for tt,na in enumerate(ooi.date.values):
    nanind=~isnan(ooi['PV'][:,tt])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.05, output='ba')
        ooi['PV_sm'][nanind,tt]=sig.filtfilt(B,A,ooi['PV'][nanind,tt].values)
    else:
        ooi['PV_sm'][:,tt]=NaN

for dd,na in enumerate(ooi.prs_mid.values):
    nanind=~isnan(ooi['PV'][dd,:])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.05, output='ba')
        ooi['PV_sm'][dd,nanind]=sig.filtfilt(B,A,ooi['PV_sm'][dd,nanind].values)
    else:
        ooi['PV_sm'][dd,:]=NaN

oom['pden_sm']=NaN*oom['pden']
for pp,na in enumerate(oom.prs.values):
    nanind=~isnan(oom['pden'][:,pp])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.05, output='ba')
        oom['pden_sm'][nanind,pp]=sig.filtfilt(B,A,oom['pden'][nanind,pp].values)
    else:
        oom['pden_sm'][:,pp]=NaN

MLcol='#41ab5d'
PVcbar=cm.plasma_r


def plotPV(moornum,axx,tit):
    hh=axx.contourf(dat.date.values,dat.depth.values,log10(dat.PV_sm[moornum,:,:].values),51,vmin=-11.25,vmax=-10,cmap=PVcbar)
    axx.contour(dat.date.values,dat.depth.values,log10(dat.PV_sm[moornum,:,:].values),levels=[-11],colors='red',linewidths=3)
    axx.plot(ML_weekly.date,ML_weekly[:,moornum].T,color=MLcol,linewidth=4)
    if moornum==4:
        axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2],colors='k',zorder=100,linewidths=3)
    else:
        axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2,d3],colors='k',zorder=100,linewidths=3)
    axx.set_title(tit,fontsize=14)
    return hh


def plotPV_ooi(axx,tit):
    axx.contourf(ooi.date.values,ooi.prs_mid.values,log10(ooi.PV_sm.values),51,vmin=-11.25,vmax=-10,cmap=PVcbar)
    contour(ooi.date.values,ooi.prs_mid.values,log10(ooi.PV_sm.values),levels=[-11],colors='red',linewidths=3)
    axx.contour(oom.date.values,oom.depth.values,oom['pden_sm'].values.T,levels=[d1,d2,d3],colors='k',zorder=99,linewidths=3)
    axx.plot(ooi_ML_weekly.date,ooi_ML_weekly.ML,color=MLcol,linewidth=3)

    axx.set_title(tit,fontsize=14)

def pltveldensec(axx,draw_moors='yes'):
        vel=axx.contourf(-osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
        # axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[27.53],colors='r',linewidths=3)
        axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=dbnds,colors='k',linewidths=3)
        axx.fill_between(-osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)

        axx.set_ylabel('depth [m]')
        if draw_moors=='yes':
            axx.axvline(-CFlon[4],color='w',linewidth=4)
            axx.axvline(-CFlon[4],color=cf5col,linewidth=2)
            # axx.axvline(-CFlon[5],color='w',linewidth=4)
            # axx.axvline(-CFlon[5],color=cf6col,linewidth=2)
            axx.axvline(-CFlon[-1],color='w',linewidth=4)
            axx.axvline(-CFlon[-1],color=m1col,linewidth=2)
            axx.axvline(-ooi_lon['fla'],color='w',linewidth=4)
            axx.axvline(-ooi_lon['fla'],color=oomcol,linewidth=2)
        return vel

osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

def PV_presplot():
    f=figure(figsize=(9,13))
    outer = gridspec.GridSpec(2, 1, height_ratios = [1,3.5])
    gs0=gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = outer[0], width_ratios=[1.5,1],wspace=0.3)
    gs1=gridspec.GridSpecFromSubplotSpec(3,1, subplot_spec = outer[1],hspace=0.2)
    ax=plt.subplot(gs0[1])
    pltveldensec(ax)
    ax.set_ylim(3100,0)
    ax.set_xlim(43.5,39)
    ax.set_ylabel('depth [m]')

    tf=14
    ax.text(42.75,1200,'CF5',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    # ax.text(42.5,1800,'CF6',color=cf6col,fontsize=tf,zorder=101,backgroundcolor='w')
    ax.text(41.7,2300,'M1',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.5,2950,'OOI',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(38.9,350,'upper IIW',fontsize=tf-2,zorder=101)
    ax.text(38.9,950,'deep IIW',fontsize=tf-2,zorder=101)
    ax.set_xlabel('Longitude [$^\circ$W]')

    bx=plt.subplot(gs0[0])
    map=SimpleMap()
    x1,y1=map(-40,53.5)
    bx.text(x1,y1,'OSNAP East',backgroundcolor='white',fontsize=14)
    mlon=-43
    xlon=-39.5
    mlat=59.5
    xlat=60.5
    map.plot([mlon,xlon,xlon,mlon,mlon],[xlat,xlat,mlat,mlat,xlat],latlon=True,color='r',linewidth=3)

    # f,(ax1,ax2,ax3)=subplots(3,1,figsize=(9,8),sharex=True)
    ax1=plt.subplot(gs1[0])
    ax2=plt.subplot(gs1[1])
    ax3=plt.subplot(gs1[2])
    plotPV_ooi(ax3,'Irminger gyre interior, OOI')
    hh=plotPV(7,ax2,'Offshore of the boundary current, M1')
    # hh=plotPV(5,ax3,'Deep boundary current mooring, CF6')
    hh=plotPV(4,ax1,'Boundary current maximum, CF5')
    cbar_ax = f.add_axes([0.93, 0.2, 0.02, 0.4])
    cbar=colorbar(hh, cax=cbar_ax,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]',ticks=arange(-9.75,-11.5,-0.25))
    ax1.set_ylim(1300,0)
    ax2.set_ylim(1500,0)
    ax1.set_yticks([1000,500,0])
    ax2.set_yticks([1500,750,0])
    ax3.set_ylim(1500,100)
    ax3.set_yticks([1500,1000,500])
    ax2.set_ylabel('depth [m]')
    ax1.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
    ax2.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
    ax3.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
    ax1.set_xticklabels('')
    ax2.set_xticklabels('')
    ax3.xaxis.set_major_locator(years)
    ax3.xaxis.set_minor_locator(threemonth)
    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    ax1.text(datetime.datetime(2015,7,1),800,'upper IIW',fontsize=14,color='k')
    ax1.text(datetime.datetime(2015,4,1),1200,'deep IIW',fontsize=14,color='k')
    sday=20
    ax2.text(datetime.datetime(2014,9,sday),600,'upper IIW',fontsize=14,color='k')
    ax2.text(datetime.datetime(2014,9,sday),1000,'deep IIW',fontsize=14,color='k')
    ax3.text(datetime.datetime(2014,9,sday),600,'upper IIW',fontsize=14,color='k')
    ax3.text(datetime.datetime(2014,9,sday),1150,'deep IIW',fontsize=14,color='k')
    ax3.text(datetime.datetime(2016,4,1),1100,'Mixed\nLayer\nDepth',fontsize=16,color=MLcol)
    savefig(figdir+'MixedLayer/paperfigs/Fig1.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig1.pdf',bbox_inches='tight')


PV_presplot()
