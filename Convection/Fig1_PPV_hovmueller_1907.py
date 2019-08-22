###############################################################################################################################
###############################################################################################################################
############################################## PV timeseries with MLDs   ####################################################
###############################################################################################################################
###############################################################################################################################
from aux_funcs import *
from map_funcs import *

dat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_10m_wML.nc')
osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

# osnap_locs=xr.Dataset(coords={'lon':osnap.LONGITUDE.values,'lat':osnap.LATITUDE.values})
# osnap_locs.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/osnap_locs.nc','w',format='netCDF4')



#### smooth out the PV and density fields for plotting purposes
dat['PV_sm']=NaN*dat['PV']
dat['pden_sm']=NaN*dat['pden']

### smooth out in depth first
for tt,na in enumerate(dat.date.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['PV'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.02, output='ba')
            dat['PV_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['PV'][mm,nanind,tt].values)
        else:
            dat['PV_sm'][mm,:,tt]=NaN
        nanind=~isnan(dat['mden'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.02, output='ba')
            dat['pden_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['mden'][mm,nanind,tt].values)
        else:
            dat['pden_sm'][mm,:,tt]=NaN

## then smooth out in time
for dd,na in enumerate(dat.depth.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['pden_sm'][mm,dd,:])
        if sum(nanind)>10:
            B, A = sig.butter(2,0.05, output='ba')
            dat['pden_sm'][mm,dd,nanind]=sig.filtfilt(B,A,dat['pden_sm'][mm,dd,nanind].values)
        else:
            dat['pden_sm'][mm,dd,:]=NaN

for dd,na in enumerate(dat.mid_depth.values):
    for mm,na in enumerate(dat.distance.values):
            nanind=~isnan(dat['PV_sm'][mm,dd,:])
            if sum(nanind)>10:
                B, A = sig.butter(2,0.1, output='ba')
                dat['PV_sm'][mm,dd,nanind]=sig.filtfilt(B,A,dat['PV_sm'][mm,dd,nanind].values)
                # dat['PV_sm'][mm,dd,nanind]=dat['PV_sm'][mm,dd,nanind].values
            else:
                dat['PV_sm'][mm,dd,:]=NaN

###### Smooth the OOI PV from profiler mooring as well!
ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb'))

#### First, fill in NaNs vertically, then horizontally
ooi['PV_fill']=NaN*ooi['PV']
for tt,na in enumerate(ooi.date.values):
    nanind=~isnan(ooi['PV'][:,tt])
    if sum(nanind)>2:
        fooi=interp1d(where(nanind)[0],ooi['PV'][nanind,tt].values,bounds_error=False)
        ooi['PV_fill'][:,tt]=fooi(range(len(ooi.prs_mid)))


for dd,na in enumerate(ooi.prs_mid.values):
    nanind=~isnan(ooi['PV_fill'][dd,:])
    if na>150:
        fooi=interp1d(where(nanind)[0],ooi['PV_fill'][dd,nanind].values,bounds_error=False)
        ooi['PV_fill'][dd,:]=fooi(range(len(ooi.date)))



ooi['PV_sm']=NaN*ooi['PV']
for tt,na in enumerate(ooi.date.values):
    nanind=~isnan(ooi['PV_fill'][:,tt])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.01, output='ba')
        ooi['PV_sm'][nanind,tt]=sig.filtfilt(B,A,ooi['PV_fill'][nanind,tt].values)
    else:
        ooi['PV_sm'][:,tt]=NaN


for dd,na in enumerate(ooi.prs_mid.values):
    nanind=~isnan(ooi['PV_sm'][dd,:])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.05, output='ba')
        ooi['PV_sm'][dd,nanind]=sig.filtfilt(B,A,ooi['PV_sm'][dd,nanind].values)
    else:
        ooi['PV_sm'][dd,:]=NaN


MLcol='#41ab5d'
PVcbar=cm.plasma_r

ML_weekly=dat.ML.resample(date='3D').mean(dim='date')

dat.ML[2,:].sel(date=slice('2015-2-1','2015-5-1')).mean(dim='date')
dat.ML[2,:].sel(date=slice('2016-2-1','2016-5-1')).mean(dim='date')

ML_weekly[2,:].plot()
[axvline(dd) for dd in [datetime.datetime(2015,2,1),datetime.datetime(2015,5,1),datetime.datetime(2016,2,1),datetime.datetime(2016,5,1)] ]

def plotPV(moornum,axx,tit,moor='na'):
    if moor=='ooi':
        log10ooi=log10(ooi.PV_sm.values)
        log10ooi[log10ooi<-11.25]=-11.5
        log10ooi[isnan(log10ooi)]=-11.5
        hh=axx.contourf(ooi.date.values,ooi.prs_mid.values,log10ooi,51,vmin=-11.25,vmax=-10,cmap=PVcbar)
        contour(ooi.date.values,ooi.prs_mid.values,log10(ooi.PV_sm.values),levels=[-11],colors='red',linewidths=3)
    else:
        hh=axx.contourf(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),51,vmin=-11.25,vmax=-10,cmap=PVcbar)
        axx.contour(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),levels=[-11],colors='red',linewidths=3)
    axx.plot(ML_weekly.date,ML_weekly[moornum,:].T,color=MLcol,linewidth=3)
    # axx.plot(dat.date,dat.ML[moornum,:].T,color=MLcol,linewidth=2)
    if moornum==4:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2],colors='k',zorder=100,linewidths=3)
    else:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2,d3],colors='k',zorder=100,linewidths=3)
    axx.set_title(tit,fontsize=14)
    return hh,denlab

def pltveldensec(axx,draw_moors='yes'):
        vel=axx.contourf(-osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
        # axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[27.53],colors='r',linewidths=3)
        axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=dbnds,colors='k',linewidths=3)
        axx.fill_between(-osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
        axx.set_ylim(3100,0)
        axx.set_xlim(43.5,39)
        axx.set_ylabel('depth [m]')
        if draw_moors=='yes':
            for dd in dat.lon:
                axx.axvline(-dd,color='w',linewidth=4)
                axx.axvline(-dd,color='r',linewidth=2)
                for zz in [50,250,500,750,1000,1500]:
                    axx.plot(-dd,zz,'k.')
                if dd>-41.5:
                    for zz in [1250,1750,2000,2250,2500]:
                        axx.plot(-dd,zz,'k.')
        return vel


def PV_presplot():
    f=figure(figsize=(9,12))
    outer = gridspec.GridSpec(2, 1, height_ratios = [3,1],hspace=0.3)
    gs0=gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = outer[1], width_ratios=[1.5,1],wspace=0.35)
    gs1=gridspec.GridSpecFromSubplotSpec(3,1, subplot_spec = outer[0],hspace=0.2)
    ax=plt.subplot(gs0[1])
    vel=pltveldensec(ax)
    caxvel=f.add_axes([0.93, 0.135, 0.01, 0.15])
    colorbar(vel,ax=ax,ticks=arange(-0.4,0.5,0.2),label='velocity [m/s]',cax=caxvel)


    tf=14
    ax.text(42.8,1200,'CF5',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(42.5,2000,'CF6',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(41.7,2300,'M1',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.5,2950,'OOI',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.8,350,'uIIW',fontsize=tf-2,zorder=101)
    ax.text(40.8,950,'dIIW',fontsize=tf-2,zorder=101)
    ax.set_xlabel('Longitude [$^\circ$W]')
    ax.set_title('e) Boundary current region')

    bx=plt.subplot(gs0[0])
    map=SimpleMap()
    x1,y1=map(-45,54.5)
    bx.set_title('d) OSNAP East')
    # bx.text(x1,y1,'OSNAP East',backgroundcolor='white',fontsize=14)
    mlon=-43
    xlon=-39.5
    mlat=59.5
    xlat=60.5
    map.plot([mlon,xlon,xlon,mlon,mlon],[xlat,xlat,mlat,mlat,xlat],latlon=True,color='r',linewidth=3)

    ax1=plt.subplot(gs1[0])
    ax2=plt.subplot(gs1[1])
    ax3=plt.subplot(gs1[2])
    hh,denlab=plotPV(3,ax3,'c) OOI: Irminger gyre interior',moor='ooi')
    labspot=datetime.datetime(2014,11,1).toordinal()
    # ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,225),(labspot,800),(labspot,1250)])
    hh,denlab=plotPV(2,ax2,'b) M1: Offshore of the boundary current')
    # labspot=datetime.datetime(2015,7,1).toordinal()
    ax2.clabel(denlab,fmt='%1.2f',manual=[(labspot,300),(labspot,600),(labspot,1100)])
    hh,denlab=plotPV(0,ax1,'a) CF5: Boundary current maximum')
    labspot=datetime.datetime(2015,7,1).toordinal()
    ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,600),(labspot,1100)])
    cbar_ax = f.add_axes([0.93, 0.425, 0.02, 0.4])
    cbar=colorbar(hh, cax=cbar_ax,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]',ticks=arange(-9.75,-11.5,-0.25))
    vmx=11.22
    cbar.ax.plot([0, 1], [(-11+vmx)/(-10+vmx)]*2, 'r',linewidth=3)

    dpmin=150
    ax1.set_ylim(1300,dpmin)
    ax2.set_ylim(1500,dpmin)
    ax1.set_yticks([1000,500])
    ax2.set_yticks([1500,1000,500])
    ax3.set_ylim(1500,dpmin)
    ax3.set_yticks([1500,1000,500])
    # ax4.set_ylim(1500,dpmin)
    # ax4.set_yticks([1500,1000,500])

    ax2.set_ylabel('depth [m]')
    for axx in [ax1,ax2,ax3]:
        axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,5)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)

    for axx in [ax1,ax2]:
        axx.set_xticklabels('')
    fs=12

    ax3.xaxis.set_minor_formatter(monthFMT)
    ax3.xaxis.set_major_formatter(yearFMT)
    labcol='k'
    ax1.text(datetime.datetime(2015,7,1),800,'upper IIW',fontsize=14,color=labcol,weight='bold')
    ax1.text(datetime.datetime(2015,4,1),1200,'deep IIW',fontsize=14,color=labcol,weight='bold')
    sday=20
    ax2.text(datetime.datetime(2014,9,sday),600,'upper IIW',fontsize=14,color=labcol,weight='bold')
    ax2.text(datetime.datetime(2014,9,sday),1000,'deep IIW',fontsize=14,color=labcol,weight='bold')
    ax3.text(datetime.datetime(2014,9,sday),600,'upper IIW',fontsize=14,color=labcol,weight='bold')
    ax3.text(datetime.datetime(2014,9,sday),1000,'deep IIW',fontsize=14,color=labcol,weight='bold')
    ax3.text(datetime.datetime(2016,4,1),1100,'Mixed\nLayer\nDepth',fontsize=16,color=MLcol,weight='bold')

    savefig(figdir+'MixedLayer/paperfigs/Fig1.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig1.pdf',bbox_inches='tight')


PV_presplot()
