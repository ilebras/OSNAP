###############################################################################################################################
###############################################################################################################################
############################################## PV timeseries with MLDs   ####################################################
###############################################################################################################################
###############################################################################################################################
from aux_funcs import *
from map_funcs import *

dat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_5m_wML.nc')

dat

# osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))
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
            else:
                dat['PV_sm'][mm,dd,:]=NaN

###### Smooth the OOI PV from profiler mooring as well!
# ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb')) pickle issues!!
# Load the netcdf now (saved in Load_Irminger_Femke.py)
ooi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/OOI_HYPM_xray_fromFemke.nc')



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



ooi['PV_sm1']=NaN*ooi['PV']
ooi['PV_sm']=NaN*ooi['PV']
for tt,na in enumerate(ooi.date.values):
    nanind=~isnan(ooi['PV_fill'][:,tt])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.01, output='ba')
        ooi['PV_sm1'][nanind,tt]=sig.filtfilt(B,A,ooi['PV_fill'][nanind,tt].values)
    else:
        ooi['PV_sm1'][:,tt]=NaN


for dd,na in enumerate(ooi.prs_mid.values):
    nanind=~isnan(ooi['PV_sm'][dd,:])
    if sum(nanind)>10:
        B, A = sig.butter(2,0.05, output='ba')
        ooi['PV_sm'][dd,nanind]=sig.filtfilt(B,A,ooi['PV_sm'][dd,nanind].values)
    else:
        ooi['PV_sm'][dd,:]=NaN


MLcol='#41ab5d' #'#10b2dd'#
PVcbar=cm.plasma_r


## Try smoothing ML, too
dat['ML_sm']=NaN*dat['ML']
for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['ML'][mm,:])
        if sum(nanind)>10:
            B, A = sig.butter(2,0.2, output='ba')
            dat['ML_sm'][mm,nanind]=sig.filtfilt(B,A,dat['ML'][mm,nanind].values)
        else:
            dat['ML_sm'][mm,:]=NaN

dat['ML_sm'][2,:].plot()

ML_weekly=dat.ML.resample(date='3D').mean(dim='date')
dencol='k'

en4=xr.open_dataset(datadir+'aux_data/EN4_surfden.nc')

en4_1416=en4.pden.sel(date=slice('2015-1-1','2016-12-1'))
en1416_janmar=en4_1416[(en4_1416['date.month']>=1)&(en4_1416['date.month']<=3)].mean(dim='date').T
en1416_febmar=en4_1416[(en4_1416['date.month']>=2)&(en4_1416['date.month']<=3)].mean(dim='date').T

en1416_janmar.plot()
en1416_febmar.plot()

osnap=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.nc')

def mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed):

    lon_0= (lon_end + lon_start)/2.0
    lat_0= - (abs(lat_end)+abs(lat_start))/2.0


    map = Basemap(lat_0=lat_0,lon_0=lon_0,llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='h',projection='stere')

    x, y = map(lon,lat)
    lon4,lat4=np.meshgrid(en4.lon.values,en4.lat.values)
    x4,y4=map(lon4,lat4)
    CS1 = map.contourf(x4,y4,en1416_janmar.values,51,cmap=cm.YlGn,vmin=27.3,vmax=27.8,extend='both')

    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000],colors='grey')

    # map.drawcoastlines()
    map.fillcontinents('darkgrey')

    # eastind=osnap.LONGITUDE>-43
    # map.plot(osnap.LONGITUDE.values[eastind],osnap.LATITUDE.values[eastind],color='k',latlon=True,linewidth=4)

    return map,x,y,CS0,CS1

def SimpleMap():

    lat_start=55
    lat_end  =67
    lon_start=-53
    lon_end  =-18

    """Get the etopo1 data"""
    etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]

    res = findSubsetIndices(lat_start-5,lat_end+5,lon_start-40,lon_end+10,lats,lons)

    lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)

    map,x,y,CS0,CS1=mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed)
    map.drawmeridians(range(lon_start+3,lon_end+10,10),labels=[0,0,0,1],linewidth=0.0001,fontsize=12)
    map.drawparallels(arange(lat_start,lat_end+2,5),labels=[1,0,0,0],linewidth=0.0001,fontsize=12)

    # [map.plot(ooi_lon[ll],ooi_lat[ll],'o',latlon=True) for ll in ooi_lat]
    map.plot(dat.lon.values,dat.lat.values,'ko',latlon=True)
    return map


def pltveldensec(axx,draw_moors='yes'):
        vel=axx.contourf(-osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.sel(TIME=slice('2014-8-1','2015-8-1')).mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4)
        # axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=dbnds,colors='k',linewidths=3)
        axx.contour(-dat.lon,dat.depth,dat.pden.sel(date=slice('2014-8-1','2015-8-1')).mean(dim='date').T,levels=dbnds,colors='k',linewidths=3)
        axx.fill_between(-osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
        axx.set_ylim(3000,0)
        axx.set_yticks([0,1000,2000,3000])
        axx.set_xlim(43.5,39.5)
        axx.set_ylabel('depth [m]')
        if draw_moors=='yes':
            for kk in CFlon:
                axx.axvline(-kk,color='k',linewidth=1)
            for dd in dat.lon:
                axx.axvline(-dd,color='w',linewidth=4)
                axx.axvline(-dd,color='r',linewidth=2)
                for zz in [50,250,500,750,1000,1500]:
                    axx.plot(-dd,zz,'k.')
                if dd>-41.5:
                    for zz in [1250,1750,2000,2250,2500]:
                        axx.plot(-dd,zz,'k.')
        return vel

def plotPV(moornum,axx,tit,moor='na',vminooi=-11.25):
    if moor=='ooi':
        log10ooi=log10(ooi.PV_sm.values)
        log10ooi[log10ooi<-11.25]=-11.5
        log10ooi[isnan(log10ooi)]=-11.5
        hh=axx.contourf(ooi.date.values[2:-2],ooi.prs_mid.values,log10ooi[:,2:-2],51,vmin=vminooi,vmax=-10,cmap=PVcbar)
        contour(ooi.date.values[2:-6],ooi.prs_mid.values,log10(ooi.PV_sm.values)[:,2:-6],levels=[-11],colors='red',linewidths=3)
    else:
        hh=axx.contourf(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),51,vmin=-11.25,vmax=-10,cmap=PVcbar)
        axx.contour(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),levels=[-11],colors='red',linewidths=3)
    # axx.plot(ML_weekly.date,ML_weekly[moornum,:].T,color=MLcol,linewidth=3)
    # axx.plot(dat.date,dat.ML_sm[moornum,:].T,color='grey',linewidth=3)
    ML_real=dat.ML_sm[moornum,:].T.copy()
    ML_real[ML_real<=190]=NaN
    if moornum==4:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2],colors=dencol,zorder=100,linewidths=3)
    elif moor=='ooi':
        denlab=axx.contour(dat.date.values[6:-22],dat.depth.values,dat['pden_sm'][moornum,:,6:-22].values,levels=[d1,d2,d3],colors=dencol,zorder=100,linewidths=3)
        axx.plot(dat.date,ML_real,color=MLcol,linewidth=3)
    else:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2,d3],colors=dencol,zorder=100,linewidths=3)
        axx.plot(dat.date,ML_real,color=MLcol,linewidth=3)
    axx.set_title(tit,fontsize=14)
    return hh,denlab

def PV_presplot():
    f=figure(figsize=(9,12))
    outer = gridspec.GridSpec(2, 1, height_ratios = [2.75,1],hspace=0.25)
    gs0=gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = outer[1], width_ratios=[1,1.1],wspace=0.8)
    gs1=gridspec.GridSpecFromSubplotSpec(3,1, subplot_spec = outer[0],hspace=0.2)
    ax=plt.subplot(gs0[0])
    vel=pltveldensec(ax)
    caxvel=f.add_axes([0.41, 0.135, 0.01, 0.15])
    colorbar(vel,ax=ax,ticks=arange(-0.4,0.5,0.2),label='velocity [m/s]',cax=caxvel)

    tf=14
    ax.text(43.2,500,'CF1-4',color='w',fontsize=tf-3,zorder=101)
    ax.text(42.8,1200,'CF5',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(42.5,2000,'CF6',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(42,2250,'CF7',color='w',fontsize=tf-3,zorder=101)
    ax.text(41.5,2400,'M1',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.5,2950,'OOI',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.8,450,'uISIW',fontsize=tf-2,zorder=101)
    ax.text(40.8,1000,'dISIW',fontsize=tf-2,zorder=101)
    ax.set_xlabel('Longitude [$^\circ$W]')
    fs=14
    ax.set_title('d) Annual mean cross-section',fontsize=fs)

    bx=plt.subplot(gs0[1])
    map=SimpleMap()
    m = plt.cm.ScalarMappable(cmap=cm.YlGn)
    m.set_array(en1416_janmar.values)
    m.set_clim(27.3,27.8)
    caxden=f.add_axes([0.93, 0.135, 0.01, 0.15])
    cbar=colorbar(m,cax=caxden,label='surface $\sigma_{\Theta}$ [kg m$^{-3}$]',boundaries=np.linspace(27.3,27.8,51),ticks=arange(27.3,28.1,0.1))

    x1,y1=map(-45,54.5)
    bx.set_title('e) Regional context',fontsize=fs)
    bx.set_xlabel('Longitude',labelpad=20)
    bx.set_ylabel('Latitude',labelpad=35)

    ax1=plt.subplot(gs1[0])
    ax2=plt.subplot(gs1[1])
    ax3=plt.subplot(gs1[2])

    hh,denlab=plotPV(3,ax3,'c) OOI: Irminger gyre interior',moor='ooi')
    labspot=datetime.datetime(2014,11,1).toordinal()
    hh,denlab=plotPV(2,ax2,'b) M1: Edge of the boundary current')
    ax2.clabel(denlab,fmt='%1.2f',manual=[(labspot,300),(labspot,600),(labspot,1100)])
    hh,denlab=plotPV(0,ax1,'a) CF5: Boundary current core')
    labspot=datetime.datetime(2015,7,1).toordinal()
    # hh,denlab=plotPV(1,ax2,'b) CF6: Deep boundary current core',moor='ooi')
    ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,600),(labspot,1100)])
    make_hh=plt.cm.ScalarMappable(cmap=PVcbar)
    make_hh.set_array(log10(dat.PV_sm[0,:,:].values))
    make_hh.set_clim(-10,-11.25)
    cbar_axhh = f.add_axes([0.93, 0.45, 0.02, 0.35])
    cbarhh=plt.colorbar(make_hh, cax=cbar_axhh,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]',ticks=arange(-9.75,-11.5,-0.25))
    cbarhh.ax.plot([0, 1], [-11]*2, 'r',linewidth=3)

    dpmin=150
    dpmax=1400
    for axx in [ax1,ax2,ax3]:
        axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,5)])
        axx.xaxis.set_major_locator(threemonth)
        axx.set_ylim(dpmax,dpmin)
        axx.set_yticks([dpmax,700,dpmin])
    for axx in [ax1,ax2]:
        axx.set_xticklabels('')
    fs=12
    ax2.set_ylabel('depth [m]',fontsize=fs)
    ax3.set_xticklabels(['October','January\n2015','April','July','October','January\n2016','April','July',])
    labcol=dencol
    ecoli='w'
    # ax1.text(datetime.datetime(2015,7,1),800,'upper ISIW',fontsize=14,color=ecoli,weight='bold')
    ax1.text(datetime.datetime(2015,7,1),800,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    ax1.text(datetime.datetime(2015,4,1),1200,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    sday=15
    ax2.text(datetime.datetime(2014,9,sday),600,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    ax2.text(datetime.datetime(2014,9,sday),1000,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    ax3.text(datetime.datetime(2014,9,sday),700,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    ax3.text(datetime.datetime(2014,9,sday),975,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    ax3.text(datetime.datetime(2016,4,1),1100,'Mixed\nLayer\nDepth',fontsize=16,color=MLcol,weight='bold')

    savefig(figdir+'MixedLayer/paperfigs/Fig1.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig1.pdf',bbox_inches='tight')


PV_presplot()

# def PPV_only():
#     f=figure(figsize=(10,9))
#     gs1=gridspec.GridSpec(3,1, hspace=0.2)
#
#     ax1=plt.subplot(gs1[0])
#     ax2=plt.subplot(gs1[1])
#     ax3=plt.subplot(gs1[2])
#     hh,denlab=plotPV(3,ax3,'OOI: Irminger gyre interior',moor='ooi')
#     labspot=datetime.datetime(2014,11,1).toordinal()
#     # ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,225),(labspot,800),(labspot,1250)])
#     hh,denlab=plotPV(2,ax2,'M1: Edge of the boundary current')
#     # labspot=datetime.datetime(2015,7,1).toordinal()
#     ax2.clabel(denlab,fmt='%1.2f',manual=[(labspot,300),(labspot,600),(labspot,1100)])
#     hh,denlab=plotPV(0,ax1,'CF5: Boundary current core')
#     labspot=datetime.datetime(2015,7,1).toordinal()
#     ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,600),(labspot,1100)])
#     make_hh=plt.cm.ScalarMappable(cmap=PVcbar)
#     make_hh.set_array(log10(dat.PV_sm[0,:,:].values))
#     make_hh.set_clim(-10,-11.25)
#     cbar_axhh = f.add_axes([0.93, 0.25, 0.02, 0.5])
#     cbarhh=colorbar(make_hh, cax=cbar_axhh,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]',ticks=arange(-9.75,-11.5,-0.25))
#     cbarhh.ax.plot([0, 1], [-11]*2, 'r',linewidth=3)
#
#     dpmin=150
#     ax1.set_ylim(1300,dpmin)
#     ax2.set_ylim(1500,dpmin)
#     ax1.set_yticks([1000,500])
#     ax2.set_yticks([1500,1000,500])
#     ax3.set_ylim(1500,dpmin)
#     ax3.set_yticks([1500,1000,500])
#
#
#     ax2.set_ylabel('depth [m]')
#     for axx in [ax1,ax2,ax3]:
#         axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,5)])
#         axx.xaxis.set_major_locator(years)
#         axx.xaxis.set_minor_locator(threemonth)
#
#     for axx in [ax1,ax2]:
#         axx.set_xticklabels('')
#     fs=12
#
#     ax3.xaxis.set_minor_formatter(monthFMT)
#     ax3.xaxis.set_major_formatter(yearFMT)
#     labcol=dencol
#     ecoli='w'
#     ax1.text(datetime.datetime(2015,7,1),800,'upper ISIW',fontsize=14,color=labcol,weight='bold')
#     ax1.text(datetime.datetime(2015,4,1),1200,'deep ISIW',fontsize=14,color=labcol,weight='bold')
#     sday=15
#     ax2.text(datetime.datetime(2014,9,sday),600,'upper ISIW',fontsize=14,color=labcol,weight='bold')
#     ax2.text(datetime.datetime(2014,9,sday),1000,'deep ISIW',fontsize=14,color=labcol,weight='bold')
#     ax3.text(datetime.datetime(2014,9,sday),700,'upper ISIW',fontsize=14,color=labcol,weight='bold')
#     ax3.text(datetime.datetime(2014,9,sday),975,'deep ISIW',fontsize=14,color=labcol,weight='bold')
#     ax3.text(datetime.datetime(2016,4,1),1100,'Mixed\nLayer\nDepth',fontsize=16,color=MLcol,weight='bold')
#
#     savefig(figdir+'MixedLayer/paperfigs/PPVonly_Fig1.png',bbox_inches='tight',dpi=300)
#     savefig(figdir+'MixedLayer/paperfigs/PPVonly_Fig1.pdf',bbox_inches='tight')
#
#
# PPV_only()
#
#
def CF5_only():
    f=figure(figsize=(10,3))
    gs1=gridspec.GridSpec(1,1, hspace=0.2)

    ax1=plt.subplot(gs1[0])
    # ax2=plt.subplot(gs1[1])
    # ax3=plt.subplot(gs1[2])
    # hh,denlab=plotPV(3,ax3,'OOI: Irminger gyre interior',moor='ooi')
    # labspot=datetime.datetime(2014,11,1).toordinal()
    # # ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,225),(labspot,800),(labspot,1250)])
    # hh,denlab=plotPV(2,ax2,'M1: Edge of the boundary current')
    # # labspot=datetime.datetime(2015,7,1).toordinal()
    # ax2.clabel(denlab,fmt='%1.2f',manual=[(labspot,300),(labspot,600),(labspot,1100)])
    hh,denlab=plotPV(0,ax1,'CF5: Boundary current core')
    labspot=datetime.datetime(2015,7,1).toordinal()
    ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,600),(labspot,1100)])
    make_hh=plt.cm.ScalarMappable(cmap=PVcbar)
    make_hh.set_array(log10(dat.PV_sm[0,:,:].values))
    make_hh.set_clim(-10,-11.25)
    cbar_axhh = f.add_axes([0.93, 0.25, 0.02, 0.5])
    cbarhh=colorbar(make_hh, cax=cbar_axhh,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]',ticks=arange(-9.75,-11.5,-0.25))
    cbarhh.ax.plot([0, 1], [-11]*2, 'r',linewidth=3)

    dpmin=150
    ax1.set_ylim(1300,dpmin)
    # ax2.set_ylim(1500,dpmin)
    ax1.set_yticks([1000,500])
    # ax2.set_yticks([1500,1000,500])
    # ax3.set_ylim(1500,dpmin)
    # ax3.set_yticks([1500,1000,500])
    #

    ax1.set_ylabel('depth [m]')
    # for axx in [ax1,ax2,ax3]:
    axx=ax1
    axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,5)])
    axx.xaxis.set_major_locator(years)
    axx.xaxis.set_minor_locator(threemonth)
    #
    # for axx in [ax1,ax2]:
    #     axx.set_xticklabels('')
    # fs=12
    #
    ax1.xaxis.set_minor_formatter(monthFMT)
    ax1.xaxis.set_major_formatter(yearFMT)
    ax1.set_xticklabels(['October','January\n2015','April','July','October','January\n2016','April','July',])
    labcol=dencol
    # ecoli='w'
    ax1.text(datetime.datetime(2015,7,1),800,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    ax1.text(datetime.datetime(2015,4,1),1200,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    # sday=15
    # ax2.text(datetime.datetime(2014,9,sday),600,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    # ax2.text(datetime.datetime(2014,9,sday),1000,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    # ax3.text(datetime.datetime(2014,9,sday),700,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    # ax3.text(datetime.datetime(2014,9,sday),975,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    # ax3.text(datetime.datetime(2016,4,1),1100,'Mixed\nLayer\nDepth',fontsize=16,color=MLcol,weight='bold')

    savefig(figdir+'MixedLayer/paperfigs/CF5only_Fig1.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/CF5only_Fig1.pdf',bbox_inches='tight')

CF5_only()

PVcbar=cm.YlGnBu

def plotPV(moornum,axx,tit,moor='na',vminooi=-11.25):
    if moor=='ooi':
        log10ooi=log10(ooi.PV_sm1.values)
        # log10ooi[log10ooi<-11.25]=-11.5
        # log10ooi[isnan(log10ooi)]=-11.5
        log10ooi[log10ooi<-12]=-12
        log10ooi[isnan(log10ooi)]=-12
        hh=axx.contourf(ooi.date.values[2:-2],ooi.prs_mid.values,log10ooi[:,2:-2],51,vmin=vminooi,vmax=-10,cmap=PVcbar)
        contour(ooi.date.values[2:-6],ooi.prs_mid.values,log10(ooi.PV_sm.values)[:,2:-6],levels=[-11.25],colors='red',linewidths=3)
    else:
        hh=axx.contourf(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),51,vmin=-11.25,vmax=-10,cmap=PVcbar)
        axx.contour(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),levels=[-11],colors='red',linewidths=3)
    # axx.plot(ML_weekly.date,ML_weekly[moornum,:].T,color=MLcol,linewidth=3)
    # axx.plot(dat.date,dat.ML_sm[moornum,:].T,color='grey',linewidth=3)
    ML_real=dat.ML_sm[moornum,:].T.copy()
    ML_real[ML_real<=190]=NaN
    if moornum==4:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2],colors=dencol,zorder=100,linewidths=3)
    elif moor=='ooi':
        # denlab=axx.contour(dat.date.values[6:-22],dat.depth.values,dat['pden_sm'][moornum,:,6:-22].values,levels=[d1,d2,d3],colors=dencol,zorder=100,linewidths=3)
        axx.plot(dat.date,ML_real,color='r',linewidth=3)
    else:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2,d3],colors=dencol,zorder=100,linewidths=3)
        axx.plot(dat.date,ML_real,color=MLcol,linewidth=3)
    axx.set_title(tit,fontsize=14)
    return hh#,denlab

def OOI_only():
    f=figure(figsize=(10,2))
    gs1=gridspec.GridSpec(1,1, hspace=0.2)

    ax1=plt.subplot(gs1[0])
    # ax2=plt.subplot(gs1[1])
    # ax3=plt.subplot(gs1[2])
    vmo=-12
    hh=plotPV(3,ax1,'',moor='ooi',vminooi=vmo)
    labspot=datetime.datetime(2014,11,1).toordinal()
    # # ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,225),(labspot,800),(labspot,1250)])
    # hh,denlab=plotPV(2,ax2,'M1: Edge of the boundary current')
    # # labspot=datetime.datetime(2015,7,1).toordinal()
    # ax2.clabel(denlab,fmt='%1.2f',manual=[(labspot,300),(labspot,600),(labspot,1100)])
    # hh,denlab=plotPV(0,ax1,'CF5: Boundary current core')
    # labspot=datetime.datetime(2015,7,1).toordinal()
    # ax1.clabel(denlab,fmt='%1.2f',manual=[(labspot,600),(labspot,1100)])
    make_hh=plt.cm.ScalarMappable(cmap=PVcbar)
    make_hh.set_array(log10(dat.PV_sm[0,:,:].values))
    make_hh.set_clim(-10,vmo)
    cbar_axhh = f.add_axes([0.93, 0.25, 0.02, 0.5])
    cbarhh=colorbar(make_hh, cax=cbar_axhh,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]',ticks=arange(-9.5,vmo-1,-0.5))
    cbarhh.ax.plot([0, 1], [-11]*2, 'r',linewidth=3)

    dpmin=150
    ax1.set_ylim(1300,dpmin)
    # ax2.set_ylim(1500,dpmin)
    ax1.set_yticks([1000,500])
    # ax2.set_yticks([1500,1000,500])
    # ax3.set_ylim(1500,dpmin)
    # ax3.set_yticks([1500,1000,500])
    #

    ax1.set_ylabel('depth [m]')
    # for axx in [ax1,ax2,ax3]:
    axx=ax1
    axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,5)])
    # axx.xaxis.set_major_locator(years)
    axx.xaxis.set_major_locator(threemonth)
    #
    # for axx in [ax1,ax2]:
    #     axx.set_xticklabels('')
    # fs=12
    #
    # ax1.xaxis.set_major_formatter(monthFMT)
    # ax1.xaxis.set_major_formatter(yearFMT)
    ax1.set_xticklabels(['October','January\n2015','April','July','October','January\n2016','April','July',])
    labcol=dencol
    # ecoli='w'
    # ax1.text(datetime.datetime(2015,7,1),800,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    # ax1.text(datetime.datetime(2015,4,1),1200,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    sday=15
    # ax2.text(datetime.datetime(2014,9,sday),600,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    # ax2.text(datetime.datetime(2014,9,sday),1000,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    # ax1.text(datetime.datetime(2014,9,sday),700,'upper ISIW',fontsize=14,color=labcol,weight='bold')
    # ax1.text(datetime.datetime(2014,9,sday),975,'deep ISIW',fontsize=14,color=labcol,weight='bold')
    ax1.text(datetime.datetime(2016,4,1),1000,'Mixed\nLayer\nDepth',fontsize=16,color='r',weight='bold')

    savefig(figdir+'MixedLayer/paperfigs/OOIonly_Fig1.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/OOIonly_Fig1.pdf',bbox_inches='tight')

OOI_only()

def velsec_only():
    f,ax=subplots(1,1,figsize=(4,3))
    vel=pltveldensec(ax)
    caxvel=f.add_axes([0.93, 0.25, 0.025, 0.5])
    colorbar(vel,ax=ax,ticks=arange(-0.4,0.5,0.2),label='velocity [m/s]',cax=caxvel,extend='both')

    tf=14
    ax.text(42.6,1200,'CF5',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(42.3,2000,'CF6',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(41.5,2300,'M1',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.3,2950,'OOI',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(41.05,450,'upper ISIW',fontsize=tf-2,zorder=101)
    ax.text(41.05,1000,'deep ISIW',fontsize=tf-2,zorder=101)
    ax.set_xlabel('Longitude [$^\circ$W]')

    savefig(figdir+'MixedLayer/paperfigs/velseconly_Fig1.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/velseconly_Fig1.pdf',bbox_inches='tight')

velsec_only()


def pltvelsec(axx,draw_moors='yes'):
        vel=axx.contourf(-osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.sel(TIME=slice('2014-8-1','2015-8-1')).mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4)
        # axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=dbnds,colors='k',linewidths=3)
        # axx.contour(-dat.lon,dat.depth,dat.pden.sel(date=slice('2014-8-1','2015-8-1')).mean(dim='date').T,levels=dbnds,colors='k',linewidths=3)
        axx.fill_between(-osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
        axx.set_ylim(3000,0)
        axx.set_yticks([0,1000,2000,3000])
        axx.set_xlim(43.5,39.5)
        axx.set_ylabel('depth [m]')
        if draw_moors=='yes':
            for kk in CFlon:
                axx.axvline(-kk,color='k',linewidth=1)
            for dd in dat.lon:
                axx.axvline(-dd,color='w',linewidth=4)
                axx.axvline(-dd,color='r',linewidth=2)
                for zz in [50,250,500,750,1000,1500]:
                    axx.plot(-dd,zz,'k.')
                if dd>-41.5:
                    for zz in [1250,1750,2000,2250,2500]:
                        axx.plot(-dd,zz,'k.')
        return vel

def velsec_only_noden():
    f,ax=subplots(1,1,figsize=(4,3))
    vel=pltvelsec(ax)
    caxvel=f.add_axes([0.93, 0.25, 0.025, 0.5])
    colorbar(vel,ax=ax,ticks=arange(-0.4,0.5,0.2),label='velocity [m/s]',cax=caxvel,extend='both')

    tf=14
    ax.text(42.6,1200,'CF5',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(42.3,2000,'CF6',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(41.5,2300,'M1',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    ax.text(40.3,2950,'OOI',color='w',fontsize=tf,zorder=101)
    ax.text(44.5,-100,'Greenland',fontsize=tf,fontweight='bold')
    # ax.text(40,-100,'Central',fontsize=tf,fontweight='bold')
    # ax.text(39,150,'Irminger\nSea',fontsize=tf,fontweight='bold')
    # ax.text(41.05,450,'upper ISIW',fontsize=tf-2,zorder=101)
    # ax.text(41.05,1000,'deep ISIW',fontsize=tf-2,zorder=101)
    ax.set_xlabel('Longitude [$^\circ$W]')

    savefig('/home/isabela/Documents/conferences_seminars/2002_OSM/presentation/figures/velseconly_Fig1_noden.pdf',bbox_inches='tight',dpi=300)

velsec_only_noden()

#
# def Denmaponly():
#     f=figure(figsize=(4,3))
#
#     map=SimpleMap()
#     m = plt.cm.ScalarMappable(cmap=cm.YlGn)
#     m.set_array(en1416_janmar.values)
#     m.set_clim(27.3,27.8)
#     caxden=f.add_axes([0.93, 0.2, 0.02, 0.6])
#     cbar=colorbar(m,cax=caxden,label='surface $\sigma_{\Theta}$ [kg m$^{-3}$]',boundaries=np.linspace(27.3,27.8,51),ticks=arange(27.3,28.1,0.1))
#
#     x1,y1=map(-45,54.5)
#
#     savefig(figdir+'MixedLayer/paperfigs/Denmaponly_Fig1.png',bbox_inches='tight',dpi=300)
#     savefig(figdir+'MixedLayer/paperfigs/Denmaponly_Fig1.pdf',bbox_inches='tight')
#
#
# Denmaponly()
