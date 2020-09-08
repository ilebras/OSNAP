from aux_funcs_2020 import *

dat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_rotnewfieldnames_2m.nc')

# extend velocity data downwards so I can get relative vorticity in the bottom triangles
dat['across track velocity (extend)']=dat['across track velocity'].copy()
# dat['potential density (extend)']=dat['potential density'].copy()
for ii,dd in enumerate(dat.distance.values):
    for jj,tt in enumerate(dat.date.values):
        velvec=dat['across track velocity'].sel(distance=dd).sel(date=tt)
        # pdenvec=dat['potential density'].sel(distance=dd).sel(date=tt)
        if sum(~isnan(velvec))!=0:
            velvec[isnan(velvec)]=velvec[~isnan(velvec)][-1]
            dat['across track velocity (extend)'][ii,:,jj]=velvec
        # if sum(~isnan(pdenvec))!=0:
        #     pdenvec[isnan(pdenvec)]=pdenvec[~isnan(pdenvec)][-1]
        #     dat['potential density (extend)'][ii,:,jj]=pdenvec

#smooth velocity data
dat['across track velocity (smooth)']=(['distance','depth','date'],dat['across track velocity (extend)'].groupby_bins('depth',range(-50,2200,100),labels=range(0,2200,100)).mean(dim='depth').interp(depth_bins=dat.depth.values))
dat['across track velocity (smooth plot)']=(['distance','depth','date'],dat['across track velocity'].groupby_bins('depth',range(-50,2200,100),labels=range(0,2200,100)).mean(dim='depth').interp(depth_bins=dat.depth.values))


#smooth pden data
dat['potential density (smooth)']=(['distance','depth','date'],dat['potential density'].groupby_bins('depth',range(-50,2200,100),labels=range(0,2200,100)).mean(dim='depth').interp(depth_bins=dat.depth.values))

## q/f = (1 + zeta/f) b_z - (u_z^2 + v_z^2)

# calculate each of these terms at mid-depth points.
dat=dat.assign_coords(mid_depth=dat.depth.values[:-1]+diff(dat.depth.values)/2)

#b_z=g/rho_0(-drho/dz)
g=9.8
rho0=1027
dat['bz_I']=(('distance','mid_depth','date'),(dat['potential density'].diff(dim='depth')/dat.depth.diff(dim='depth')).values*g/rho0)
# dat['uz2_I']=(('distance','mid_depth','date'),((dat['along track velocity'].diff(dim='depth')/dat.depth.diff(dim='depth')).values)**2)
dat['vz2_I']=(('distance','mid_depth','date'),((dat['across track velocity (smooth)'].diff(dim='depth')/dat.depth.diff(dim='depth')).values)**2)

dat['bz']=(('distance','depth','date'),dat['bz_I'].interp(mid_depth=dat.depth.values))
dat['bz']=dat['bz'].where(dat['bz']>0)

# dat['uz2']=(('distance','depth','date'),dat['uz2_I'].interp(mid_depth=dat.depth.values))
dat['vz2']=(('distance','depth','date'),dat['vz2_I'].interp(mid_depth=dat.depth.values))

dat=dat.assign_coords(mid_distance=dat.distance.values[:-1]+diff(dat.distance.values)/2)


dat['zeta_f_I']=(('mid_distance','depth','date'),(dat['across track velocity (smooth)'].diff(dim='distance')/dat.distance.diff(dim='distance')/1e3/sw.f(60)).values)
dat['zeta_f']=(('distance','depth','date'),(dat['zeta_f_I'].interp(mid_distance=dat.distance.values)))
dat['zeta_f'][-1,:,:]=dat['zeta_f_I'][-1,:,:]
dat['zeta']=sw.f(60)*dat['zeta_f']

dat['bx_I']=(('mid_distance','depth','date'),(dat['potential density'].diff(dim='distance')/dat.distance.diff(dim='distance')*g/rho0/1e3).values)
dat['bx']=(('distance','depth','date'),(dat['bx_I'].interp(mid_distance=dat.distance.values)))

dat['bz_norm']=dat['bz']*(1+dat['zeta_f'])
dat['bz_norm_f']=dat['bz_norm']*sw.f(60)
dat['bz_anom']=dat['bz_norm']-dat['bz']
dat['bz_anom_f']=dat['bz_norm_f']-dat['bz']*sw.f(60)

dat['PV/f']=(dat['bz_norm']-dat['vz2'])
dat['PV']=sw.f(60)*(dat['bz_norm']-dat['vz2'])
dat['log(PV)']=log10(dat['PV'])
dat['log(fbz)']=log10(sw.f(60)*dat['bz'])

dat['zeta_f'].mean(dim='date').plot()
dat.PV.mean(dim='date').plot()

era5=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/hourly_downloaded_200618/ERA5_hourly_SI_ilebras.nc')

era5=era5.rename({'time':'date'})

tau_daily=era5['tau'].mean(dim='longitude').resample(date='1D').mean(dim='date')

dat['FBwind']=tau_daily*dat['bx'].sel(depth=0)/rho0/sw.f(60)

def plot_PV_windcomp_3panel(moornum):
    f,axx=subplots(3,1,figsize=(12,9),sharex=True)
    # (dat['FBwind'][:,moornum]*4900/2.5e-6).plot(ax=axx[0],label='Wind',color='grey',linewidth=2)
    FBwind_sm=(dat['FBwind'][:,moornum]).resample(date='2D').mean(dim='date')
    FBwind_sm.plot(ax=axx[0],label='Wind',color='black')
    (era5['hf'].mean(dim='longitude')/4900*2.5e-6).resample(date='1D').mean(dim='date').plot(ax=axx[0],label='Buoyancy',color='C2') #Scaled as in Thomas and Lee (2005)
    axx[0].set_title('CF'+str(moornum+1))
    axx[0].set_xlabel('')
    axx[0].set_ylabel('PV extraction [m$^{2}$ s$^{-3}$]')
    axx[0].axhline(0,color='k')
    thresh=-1e-6
    axx[0].axhline(thresh,color='grey',linestyle='--')
    PVdates=dat.date.resample(date='2D').mean(dim='date')[FBwind_sm<thresh]
    for ii,dd in enumerate([450,650]):
        jj=ii+1
        dind=int(dd/2)
        axx[jj].plot(dat.date,dat['PV'][moornum,dind,:],label='PV$=(f+\zeta) b_z-f\ v_z^2$',color='red',linewidth=2,zorder=2)
        axx[jj].plot(dat.date,sw.f(60)*(dat.bz[moornum,dind,:]),label='$f\ b_z$',color='black',zorder=2)
        axx[jj].plot(dat.date,(dat.bz_anom_f[moornum,dind,:]),label='$\zeta \ b_z$',color='C1',zorder=2)
        axx[jj].plot(dat.date,-sw.f(60)*dat.vz2[moornum,dind,:],label='$-f\ v_z^2$',color='blue',zorder=2)
        axx[jj].axhline(0,color='grey')
        axx[jj].set_ylabel('[s$^{-3}$] at '+str(dd)+'m' )

    axx[0].legend(loc='lower left',ncol=2,framealpha=1)
    axx[1].legend(loc='upper left',ncol=4,framealpha=1)
    for axi in axx:
        for pp in PVdates.values:
            axi.axvline(pp,color='grey',zorder=1)
        # [axi.axvline(dd,color='orange',zorder=1) for dd in [datetime.datetime(2017,1,18),datetime.datetime(2017,1,26)]]
        axi.axvspan(datetime.datetime(2017,1,1),datetime.datetime(2017,2,1),color='yellow',zorder=1)
    axx[1].set_ylim(-2e-10,1.5e-9)
    axx[2].set_ylim(-2e-10,8e-10)
    xlim([datetime.datetime(2014,8,15),datetime.datetime(2018,8,15)])
    savefig(figdir+'SI/FBwind_PVcomp_CF'+str(moornum+1)+'_3panel.png',bbox_inches='tight')
    savefig(figdir+'SI/FBwind_PVcomp_CF'+str(moornum+1)+'_3panel.pdf',bbox_inches='tight')
    if moornum==4:
        axx[0].set_title('PV extraction and PV between instruments at CF5',fontweight='bold',fontsize=14)
        savefig(figdir+'SI/paper_prep/FBwind_PVcomp_CF'+str(moornum+1)+'_3panel.png',bbox_inches='tight')
        savefig(figdir+'SI/paper_prep/FBwind_PVcomp_CF'+str(moornum+1)+'_3panel.pdf',bbox_inches='tight')
    return PVdates

pvdates=plot_PV_windcomp_3panel(4)

# for mm in range(5,7):
#     plot_PV_windcomp_3panel(mm)

pvcheck=array([pvdates.date[ii].values for ii in [0,4,8,12]],dtype='datetime64[D]')

def make_uni(varit,cmapo,vmini,vmaxi,levs):
    uni[varit]={}
    uni[varit]['cmap']=cmapo
    uni[varit]['vmin']=vmini
    uni[varit]['vmax']=vmaxi
    uni[varit]['levels']=levs

    return uni

colors = [(12,44,132),(78,179,211),(237,248,177) ,(247,104,161),(129,15,124)]
sal_cmap = make_cmap(colors,position=[0,0.9,0.92,0.99,1],bit=True)
sal_levs=[34,34.2,34.4,34.6,34.8,34.85, 34.9,34.91, 34.92,34.93, 34.94,34.95, 34.96,34.97, 34.98,34.99,35]
uni=make_uni('sal',sal_cmap,34,35,sal_levs)
uni=make_uni('uacross',cm.Purples_r,-0.6,0.6,arange(-0.6,0.05,0.05))


def plot_vel_logPVsec():
    d1='2017-1-15'
    d2='2017-2-15'
    distlim=slice(25,120)
    ymax=750
    depthlim=slice(0,ymax)
    colvec=['red','orange']
    f,axi=subplots(1,1,sharex=True,sharey=True,figsize=(7,4))
    for ii,var in enumerate(['log(PV)','log(fbz)']):
        col=axi.contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                         dat[var].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                         linewidths=3,colors=colvec[ii],levels=[-9.9])
    col=axi.contourf(dat.distance.sel(distance=distlim),dat['depth'],
                dat['across track velocity (smooth plot)'].sel(distance=distlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                cmap=cm.Purples_r,levels=arange(-0.6,0.05,0.05),linewidths=3,extend='min')
    axi.contour(dat.distance.sel(distance=distlim),dat['depth'],
                dat['potential density'].sel(distance=distlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                colors='k',levels=[27.5,27.6,27.65,27.73],linewidths=3)
    # axi.set_title(d1+' --- '+d2)
    for ii,dd in enumerate(dat.distance[3:8]):
        if ii==0:
            for zz in [50,100,250,400]:
                axi.plot(dd,zz,'o',color='w',markersize=13)
                axi.plot(dd,zz,'o',color='k',markersize=10)
        else:
            for zz in [50,100,250,500,750]:
                axi.plot(dd,zz,'o',color='w',markersize=13)
                axi.plot(dd,zz,'o',color='k',markersize=10)
        axi.axvline(dd,color='k')
        if ii==4:
            axi.text(dd-1.5,-20,'M1',color='k',fontsize=12)
        else:
            axi.text(dd-1.5,-20,'CF'+str(ii+4),color='k',fontsize=12)
    colorbar(col,label='velocity [m/s]')
    axi.set_ylabel('depth [m]')
    axi.set_xlabel('distance [km]')
    axi.set_facecolor('grey')
    ylim(ymax,0)
    # xlim(35,80)
    savefig(figdir+'SI/paper_prep/Vel_PV_sec.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/Vel_PV_sec.pdf',bbox_inches='tight')

plot_vel_logPVsec()


uni=make_uni('log',cm.inferno,-11,-8,arange(-11,-8,0.25))


def plot_vel_logPVsec_invertcol():
    d1='2017-1-15'
    d2='2017-2-15'
    distlim=slice(25,120)
    ymax=750
    depthlim=slice(0,ymax)
    f,axi=subplots(1,1,sharex=True,sharey=True,figsize=(7,4))
    var='log(PV)'
    col=axi.contourf(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                         dat[var].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                         cmap=uni['log']['cmap'],levels=uni['log']['levels'])
    axi.contour(dat.distance.sel(distance=distlim),dat['depth'],
                dat['across track velocity (smooth plot)'].sel(distance=distlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                colors='k',levels=arange(-0.6,0.05,0.05),linewidths=3)
    for ii,dd in enumerate(dat.distance[3:8]):
        if ii==0:
            for zz in [50,100,250,400]:
                axi.plot(dd,zz,'o',color='w',markersize=13)
                axi.plot(dd,zz,'o',color='k',markersize=10)
        else:
            for zz in [50,100,250,500,750]:
                axi.plot(dd,zz,'o',color='w',markersize=13)
                axi.plot(dd,zz,'o',color='k',markersize=10)
        axi.axvline(dd,color='k')
        if ii==4:
            axi.text(dd-1.5,-20,'M1',color='k',fontsize=12)
        else:
            axi.text(dd-1.5,-20,'CF'+str(ii+4),color='k',fontsize=12)
    colorbar(col,label='log(PV) [s$^{-3}$]')
    axi.set_ylabel('depth [m]')
    axi.set_xlabel('distance [km]')
    axi.set_facecolor('grey')
    ylim(ymax,0)
    # xlim(35,80)
    savefig(figdir+'SI/paper_prep/Vel_PV_sec_invert.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/Vel_PV_sec_invert.pdf',bbox_inches='tight')

plot_vel_logPVsec_invertcol()



cf5=xr.open_dataset(datadir+'OSNAP2018recovery/Hourly_netcdf/CF5_mcat_2018recovery_hourlymerged.nc').sortby('DEPTH').sel(TIME=slice('2017-1-18','2017-1-27'))
# cf6=xr.open_dataset(datadir+'OSNAP2018recovery/Hourly_netcdf/CF6_mcat_2018recovery_hourlymerged.nc').sortby('DEPTH').sel(TIME=slice('2017-1-18','2017-1-26'))

def cf_hourly_meat(cfvar,axi):
        ymax=345
        for dd in cfvar.DEPTH.sel(DEPTH=slice(0,ymax)):
            col=axi.scatter(cfvar.TIME.values,-gsw.z_from_p(cfvar.PRES.sel(DEPTH=dd).values,60),c=cfvar.PSAL.sel(DEPTH=dd).values,vmin=uni['sal']['vmin'],vmax=uni['sal']['vmax'],s=60,zorder=3,cmap=uni['sal']['cmap'])
        axi.set_xticks([datetime.datetime(2017,1,18),datetime.datetime(2017,1,20),datetime.datetime(2017,1,22),datetime.datetime(2017,1,24),datetime.datetime(2017,1,26)])
        axi.axvspan(datetime.datetime(2017,1,20),datetime.datetime(2017,1,21),color='gray',zorder=1,alpha=0.3)
        axi.axvspan(datetime.datetime(2017,1,24),datetime.datetime(2017,1,25),color='gray',zorder=1,alpha=0.3)
        axi.set_ylim(ymax,40)
        axi.set_ylabel('depth [m]')
        axi.set_xlim(datetime.datetime(2017,1,18),datetime.datetime(2017,1,27))
        return col

def cf_hourly(cfvar,mnum):
    f,axi=subplots(1,1,figsize=(9,5))
    col=cf_hourly_meat(cfvar,axi)
    colorbar(col,label='salinity')
    title('CF'+mnum)
    savefig(figdir+'SI/paper_prep/CF'+mnum+'_hourlysalinity_scatter.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/CF'+mnum+'_hourlysalinity_scatter.pdf',bbox_inches='tight')

cf_hourly(cf5,'5')
# cf_hourly(cf6,'6')


uni['uacross']['vmin']=-0.75
uni['uacross']['vmax']=0
uni['uacross']['levels']=arange(-0.75,0.01,0.05)

def plot_sect_meat(axi,var,var2,d1,d2):
    distlim=slice(30,80)
    ymax=750
    depthlim=slice(0,ymax)
    col=axi.contourf(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                         dat[var].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                         vmin=uni[var2]['vmin'],vmax=uni[var2]['vmax'],cmap=uni[var2]['cmap'],levels=uni[var2]['levels'],extend='both')
    if var2=='sal':
        axi.contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                    dat['potential density'].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                    colors='k',levels=arange(26,29,0.1),linewidths=3)
    depthlim=slice(400,ymax)
    # colvec=['red','orange']
    # for ii,var in enumerate(['log(PV)','log(fbz)']):
    #     axi.contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
    #                      dat[var].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
    #                      linewidths=3,colors=colvec[ii],levels=[-9.9])
    for ii,dd in enumerate(dat.distance[3:7]):
        if ii==0:
            for zz in [50,100,250,400]:
                axi.plot(dd,zz,'o',color='w',markersize=13)
                axi.plot(dd,zz,'o',color='k',markersize=10)
        else:
            for zz in [50,100,250,500,750]:
                axi.plot(dd,zz,'o',color='w',markersize=13)
                axi.plot(dd,zz,'o',color='k',markersize=10)
        axi.axvline(dd,color='k')
        if var2=='sal':
            axi.text(dd-1.5,-40,'CF'+str(ii+4),color='k',fontsize=12)
    axi.set_ylim(750,0)
    axi.set_facecolor('grey')
    return col


def zoom_panels():
    fig=figure(figsize=(10,10))
    outer=gridspec.GridSpec(2,1, height_ratios=[1,1.5],hspace=0.5)
    gs1=gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
    ax1=fig.add_subplot(gs1[0])
    ax1.set_title('CF5 instrument record',fontsize=14,fontweight='bold')
    cf_hourly_meat(cf5,ax1)
    gs2=gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec = outer[1],wspace=0.1)
    axx={}
    for ii in range(2):
        for jj in range(2):
            axx[ii,jj]=fig.add_subplot(gs2[ii,jj])
    axx[0,0].text(48,-175,'Before 2017 event',fontsize=14,fontweight='bold')
    axx[0,1].text(48,-175,'After 2017 event',fontsize=14,fontweight='bold')
    plot_sect_meat(axx[0,0],'salinity','sal','2017-1-20','2017-1-20')
    plot_sect_meat(axx[1,0],'across track velocity (smooth plot)','uacross','2017-1-20','2017-1-20')
    salcol=plot_sect_meat(axx[0,1],'salinity','sal','2017-1-24','2017-1-24')
    velcol=plot_sect_meat(axx[1,1],'across track velocity (smooth plot)','uacross','2017-1-24','2017-1-24')
    sax=fig.add_axes([0.95,0.4,0.02,0.3])
    vax=fig.add_axes([0.95,0.125,0.02,0.15])
    colorbar(salcol,cax=sax,label='salinity')
    colorbar(velcol,cax=vax,label='velocity [m/s]',ticks=arange(-0.8,0.1,0.2))
    [axi.set_ylabel('depth [m]') for axi in [axx[0,0],axx[1,0]]]
    [axi.set_xlabel('distance [km]') for axi in [axx[1,0],axx[1,1]]]
    [axi.set_xticklabels('') for axi in [axx[0,0],axx[0,1]]]
    [axi.set_yticklabels('') for axi in [axx[0,1],axx[1,1]]]
    savefig(figdir+'SI/paper_prep/Zoom_SI_event.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/Zoom_SI_event.pdf',bbox_inches='tight')


zoom_panels()

#################################################################################
#########################  Transport analysis ###################################
#################################################################################

daily=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid_2m.nc')
daily=daily.where(daily['across track velocity']!=0)
daily=daily.sel(date=slice('2014-8-15','2018-9-1')) #just to avoid some zeros at the beginning which mess up filtering.

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))
daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3

onesxr=daily.salinity/daily.salinity

sep=9
d1=27.65
d2=27.73
d3=27.77
xport={}
xport['egic']=daily['xport'][sep:,:,:].sum('depth').sum('distance')

xport['upper limb']=daily['xport'].where((daily['potential density']<27.56)).sum('distance').sum('depth').where(xport['egic']<-5).interpolate_na(dim='date')
xport['lower limb']=daily['xport'].where((daily['potential density']>=27.56)).sum('distance').sum('depth').where(xport['egic']<-5).interpolate_na(dim='date')
xport['egic']=xport['egic'].where(xport['egic']<-5).interpolate_na(dim='date')


# xport['upper']=daily['xport over 27.8'].where((daily['potential density']<d2)&(daily['potential density']>=d1)).sum('distance').sum('depth')
# xport['deep']=daily['xport over 27.8'].where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')
#
# area={}
# area['upper']=(onesxr.where((daily['potential density']<d2)&(daily['potential density']>=d1))[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
# area['deep']=(onesxr.where((daily['potential density']<d3)&(daily['potential density']>=d2))[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')

depthdiff=int(unique(diff(dat.depth)))
distdiff=diff(dat.distance)
distdiff_i=hstack((distdiff[0]/2,(distdiff[:-1]+distdiff[1:])/2,distdiff[-1]/2))

logPVlim=-9.9

low_PV_thickness=(onesxr.where((dat['log(PV)']<logPVlim)&(dat.depth>=500)&(dat.depth<=750))*depthdiff).sum(dim='depth')
low_PV_area=low_PV_thickness.T*distdiff_i

dat

CF5_lowPV=(low_PV_area*dat['across track velocity (smooth)'].sel(depth=slice(500,750)).mean(dim='depth')/1e3).isel(distance=4)
CF5_lowPV

uppercol='#ae017e'
deepcol='#1c9099'
pvcol='#fc4e2a'

def trans_fig():
    f,axx=subplots(1,1,figsize=(12,3),sharex=True)
    axx.plot(daily.date,xport['egic'],color='k')
    axi=axx
    # for axi in axx:
    axi.set_xlabel('')
    axi.set_ylabel('')
    axi.set_xlim(datetime.datetime(2014,8,15),datetime.datetime(2018,8,15))
    for dd in pvdates.values:
        axi.axvline(dd,color='grey')
        axi.axvspan(datetime.datetime(2017,1,1),datetime.datetime(2017,2,1),color='yellow',zorder=1)
    axx.set_ylabel('Volume transport [Sv]')
    savefig(figdir+'SI/paper_prep/EGIC_trans.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/EGIC_trans.pdf',bbox_inches='tight')

trans_fig()

def trans_fig_overturning():
    f,axx=subplots(2,1,figsize=(12,6),sharex=True)
    axx[0].plot(daily.date,xport['egic'],color='k',label='')
    axx[1].plot(CF5_lowPV.date,CF5_lowPV,color='r',label='upper limb')
    # axx[1].plot(daily.date,xport['lower limb'],color='b',label='lower limb')
    for axi in axx:
        axi.set_xlabel('')
        axi.set_ylabel('')
        axi.set_xlim(datetime.datetime(2014,8,15),datetime.datetime(2018,8,15))
        for dd in pvdates.values:
            axi.axvline(dd,color='grey')
            axi.axvspan(datetime.datetime(2017,1,1),datetime.datetime(2017,2,1),color='yellow',zorder=1)
        axi.set_ylabel('Volume transport [Sv]')
    # savefig(figdir+'SI/paper_prep/EGIC_trans.png',bbox_inches='tight')
    # savefig(figdir+'SI/paper_prep/EGIC_trans.pdf',bbox_inches='tight')

trans_fig_overturning()

def trans_fig_wind():
    f,axx=subplots(2,1,figsize=(12,6),sharex=True)
    # tau_daily.plot(ax=axx[0])
    tau_daily.resample(date='2D').mean(dim='date').plot(ax=axx[0])
    axx[1].plot(daily.date,xport['egic'],color='k')
    # [axx[1].plot(low_PV_area.date,low_PV_area[:,ii]) for ii in range(4,5)]
    # xport_pv['deep'].plot(ax=axx[2],color=deepcol)
    for axi in axx:
        axi.set_xlabel('')
        axi.set_ylabel('')
        axi.set_xlim(datetime.datetime(2014,8,15),datetime.datetime(2018,8,15))
        for dd in pvdates.values:
            axi.axvline(dd,color='grey')
            axi.axvspan(datetime.datetime(2017,1,1),datetime.datetime(2017,2,1),color='yellow',zorder=1)
    axx[0].set_ylabel('Along-current\nwind stress [N m$^{-2}$]')
    axx[1].set_ylabel('Volume transport [Sv]')
    axx[0].set_title('')
    axx[0].axhline(0,color='grey')
    # savefig(figdir+'SI/paper_prep/EGIC_trans.png',bbox_inches='tight')
    # savefig(figdir+'SI/paper_prep/EGIC_trans.pdf',bbox_inches='tight')

trans_fig_wind()



def den_plus(d1,d2,titword):
    distlim=slice(30,125)
    ymax=750
    depthlim=slice(0,ymax)
    fig=figure()
    grd=fig.add_gridspec(2,1, height_ratios=[1,5],hspace=0.1)
    axx={}
    for ii in range(2):
        axx[ii]=fig.add_subplot(grd[ii])
        axx[ii].set_xlim(35,80)
    dat['bx'].sel(distance=distlim).sel(depth=0).sel(date=slice(d1,d2)).mean(dim='date').plot(color='k',marker='o',ax=axx[0],linewidth=3)
    axx[0].set_ylim(0,6e-7)
    axx[0].set_xticklabels('')
    for dd in dat.distance[4:7]:
        for zz in [100,250,500,750]:
            axx[1].plot(dd,zz,'o',color='w',markersize=13)
            axx[1].plot(dd,zz,'o',color='k',markersize=10)
        axx[1].axvline(dd,color='k')
    # var2='log'
    col=axx[1].contourf(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=slice(0,ymax)),
                         dat['salinity'].sel(distance=distlim).sel(depth=slice(0,ymax)).sel(date=slice(d1,d2)).mean(dim='date').T,
                         # vmin=uni[var2]['vmin'],vmax=uni[var2]['vmax'],cmap=uni[var2]['cmap'])
                         levels=arange(34,35,0.05),vmin=34,vmax=35,extend='both')
                         # cmap=cm.RdYlBu_r)
    axx[1].contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                dat['across track velocity (smooth plot)'].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                colors='magenta',levels=[-0.5,-0.4,-0.3],linewidths=3,zorder=3)
    fax=fig.add_axes([0.95,0.15,0.02,0.5])
    colorbar(col,cax=fax,label='salinity')
    axx[1].contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
            dat['potential density'].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
            colors='k',levels=[27.2,27.3,27.4,27.5,27.6,27.65,27.73],linewidths=2)
    axx[1].set_ylim(ymax,0)
    axx[1].set_xlabel('distance [km]')
    axx[1].set_ylabel('depth [m]')
    axx[0].set_title(d1)
    axx[1].set_facecolor('k')
    savefig(figdir+'SI/paper_prep/Zoomsec_'+titword+'.png',bbox_inches='tight')
    savefig(figdir+'SI/paper_prep/Zoomsec_'+titword+'.pdf',bbox_inches='tight')

den_plus('2017-1-20','2017-1-20','before')
den_plus('2017-1-24','2017-1-24','after')

def plot_sect(var,var2,d1,d2):
    distlim=slice(30,100)
    depthlim=slice(0,750)
    f,axi=subplots(1,1,sharex=True,sharey=True,figsize=(7,4))
    if var2=='sal':
        col=axi.contourf(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                         dat[var].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                         levels=sal_levs,cmap=uni[var2]['cmap'])
    else:
        col=axi.contourf(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                         dat[var].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                         vmin=uni[var2]['vmin'],vmax=uni[var2]['vmax'],cmap=uni[var2]['cmap'],levels=linspace(uni[var2]['vmin'],uni[var2]['vmax'],10))
    # axi.contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
    #             dat['across track velocity (smooth)'].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
    #             colors='k',levels=arange(-1,1,0.1),linewidths=3)
    axi.contour(dat.distance.sel(distance=distlim),dat['depth'].sel(depth=depthlim),
                dat['potential density'].sel(distance=distlim).sel(depth=depthlim).sel(date=slice(d1,d2)).mean(dim='date').T,
                colors='k',levels=[27.5,27.6,27.65,27.73],linewidths=3)
    axi.set_title(d1+' --- '+d2)
    for dd in dat.distance[4:7]:
        for zz in [100,250,500,750]:
            axi.plot(dd,zz,'o',color='w',markersize=13)
            axi.plot(dd,zz,'o',color='k',markersize=10)
        axi.axvline(dd,color='k')
    colorbar(col,label=var)
    axi.set_ylabel('depth [m]')
    axi.set_xlabel('distance [km]')
    ylim(750,0)
    xlim(35,80)
    savefig(figdir+'SI/Section_'+var+d1+'--'+d2+'.png',bbox_inches='tight')



uni=make_uni('PV',cm.RdBu_r,-1e-5,1e-5)

d1='2017-1-20'
d2='2017-1-22'
var2s=['PV','log','log']
for ii,var in enumerate(['PV','log(PV)','log(fbz)']):
    plot_sect(var,var2s[ii],d1,d2)

# def plot_PV_windcomp(moornum):
#     f,axx=subplots(2,1,figsize=(20,8),sharex=True)
#     # (dat['FBwind'][:,moornum]*4900/2.5e-6).plot(ax=axx[0],label='Wind',color='grey',linewidth=2)
#     FBwind_sm=(dat['FBwind'][:,moornum]).resample(date='2D').mean(dim='date')
#     FBwind_sm.plot(ax=axx[0],label='Wind',color='black')
#     (era5['hf'].mean(dim='longitude')/4900*2.5e-6).resample(date='1D').mean(dim='date').plot(ax=axx[0],label='Buoyancy',color='C2') #Scaled as in Thomas and Lee (2005)
#     axx[0].set_title('CF'+str(moornum+1))
#     axx[0].set_xlabel('')
#     axx[0].set_ylabel('PV extraction [m$^{2}$ s$^{-3}$]')
#     axx[0].axhline(0,color='k')
#     thresh=-1e-6
#     axx[0].axhline(thresh,color='k',linestyle='--')
#     PVdates=dat.date.resample(date='2D').mean(dim='date')[FBwind_sm<thresh]
#     for pp in PVdates.values:
#         axx[1].axvline(pp,color='grey')
#     dd=650
#     dind=int(dd/2)
#     axx[1].plot(dat.date,(dat.bz[moornum,dind,:]-dat.vz2[moornum,dind,:]),label='$b_z(1+\zeta/f)-v_z^2$',color='red',linewidth=2)
#     axx[1].plot(dat.date,(dat.bz[moornum,dind,:]),label='$b_z$',color='black')
#     axx[1].plot(dat.date,(dat.bz_anom[moornum,dind,:]),label='$b_z \ \zeta /f$',color='C1')
#     axx[1].plot(dat.date,-dat.vz2[moornum,dind,:],label='$-v_z^2$',color='blue')
#     axx[1].axhline(0,color='grey')
#     axx[1].set_ylabel('PV at 650m [s$^{-2}$]')
#     for axi in axx:
#         axi.legend(loc=(1.01,0.4))
#     axx[1].set_ylim(-0.1e-5,0.7e-5)
#     xlim([datetime.datetime(2014,8,15),datetime.datetime(2018,8,15)])
#     savefig(figdir+'SI/FBwind_PVcomp_CF'+str(moornum+1)+'.png',bbox_inches='tight')
#     return PVdates

# def plot_PV_windcomp_4panel(moornum):
#     f,axx=subplots(4,1,figsize=(20,16),sharex=True)
#     # (dat['FBwind'][:,moornum]*4900/2.5e-6).plot(ax=axx[0],label='Wind',color='grey',linewidth=2)
#     FBwind_sm=(dat['FBwind'][:,moornum]).resample(date='2D').mean(dim='date')
#     FBwind_sm.plot(ax=axx[0],label='Wind',color='black')
#     (era5['hf'].mean(dim='longitude')/4900*2.5e-6).resample(date='1D').mean(dim='date').plot(ax=axx[0],label='Buoyancy',color='C2') #Scaled as in Thomas and Lee (2005)
#     axx[0].set_title('CF'+str(moornum+1))
#     axx[0].set_xlabel('')
#     axx[0].set_ylabel('PV extraction [m$^{2}$ s$^{-3}$]')
#     axx[0].axhline(0,color='k')
#     thresh=-1e-6
#     axx[0].axhline(thresh,color='k',linestyle='--')
#     PVdates=dat.date.resample(date='2D').mean(dim='date')[FBwind_sm<thresh]
#     for ii,dd in enumerate([200,450,700]):
#         jj=ii+1
#         dind=int(dd/2)
#         axx[jj].plot(dat.date,(dat.bz[moornum,dind,:]-dat.vz2[moornum,dind,:]),label='$b_z(1+\zeta/f)-v_z^2$',color='red',linewidth=2)
#         axx[jj].plot(dat.date,(dat.bz[moornum,dind,:]),label='$b_z$',color='black')
#         axx[jj].plot(dat.date,(dat.bz_anom[moornum,dind,:]),label='$b_z \ \zeta /f$',color='C1')
#         axx[jj].plot(dat.date,-dat.vz2[moornum,dind,:],label='$-v_z^2$',color='blue')
#         axx[jj].axhline(0,color='grey')
#         axx[jj].set_ylabel('PV at '+str(dd)+'m [s$^{-2}$]')
#         for pp in PVdates.values:
#             axx[jj].axvline(pp,color='grey')
#     for axi in [axx[0],axx[1]]:
#         axi.legend(loc=(1.01,0.4))
#     axx[1].set_ylim(-2e-5,4e-5)
#     axx[2].set_ylim(-0.2e-5,1.25e-5)
#     axx[3].set_ylim(-0.1e-5,0.6e-5)
#     xlim([datetime.datetime(2014,8,15),datetime.datetime(2018,8,15)])
#     savefig(figdir+'SI/FBwind_PVcomp_CF'+str(moornum+1)+'_3panel.png',bbox_inches='tight')
#     return PVdates
#
# def plot_FBwind():
#     f,axx=subplots(4,1,figsize=(20,12),sharex=True)
#     dat['bx'].sel(depth=200).isel(distance=4).plot(ax=axx[0],label='CF5')
#     dat['bx'].sel(depth=200).isel(distance=5).plot(ax=axx[0],label='CF6')
#     dat['bx'].sel(depth=200).isel(distance=6).plot(ax=axx[0],label='CF7')
#     axx[0].legend()
#     axx[0].set_ylabel('bx at 200m')
#     dat['bx'].sel(depth=0).isel(distance=4).plot(ax=axx[1])
#     dat['bx'].sel(depth=0).isel(distance=5).plot(ax=axx[1])
#     dat['bx'].sel(depth=0).isel(distance=6).plot(ax=axx[1])
#     axx[1].set_ylabel('bx at 0m')
#     axx[2].axhline(0,color='grey')
#     tau_daily.plot(ax=axx[2],color='C2')
#     dat['FBwind'][:,4].plot(ax=axx[3])
#     dat['FBwind'][:,5].plot(ax=axx[3])
#     dat['FBwind'][:,6].plot(ax=axx[3])
#     for axi in axx:
#         axi.set_xlabel('')
#         axi.set_title('')
#         axi.set_xlim(dt.datetime(2014,8,15),dt.datetime(2018,8,15))
#     savefig(figdir+'SI/FBwind_breakdown.png',bbox_inches='tight')
#
# plot_FBwind()
# PVmap=cm.YlGnBu
# #### plot time series of PV (components at moorings CF4-CF6)
# def plot_pv_comps(moornum):
#     f,axx=subplots(4,1,sharex=True,sharey=True,figsize=(10,10))
#     dlim=150
#     vmaxi=0.5
#     (dat.bz[moornum,dlim:,:]*1e5).plot(ax=axx[0],vmin=0,vmax=vmaxi,cbar_kwargs={"label":'$b_z$'},cmap=PVmap)
#     (dat.bz_anom[moornum,dlim:,:]*1e5).plot(ax=axx[1],vmin=0,vmax=vmaxi/2,cbar_kwargs={"label":'$b_z \ \zeta /f$'},cmap=PVmap)
#     ((dat.vz2[moornum,dlim:,:])*1e5).plot(ax=axx[2],vmin=0,vmax=vmaxi/2,cbar_kwargs={"label":'$v_z^2$'},cmap=PVmap)
#     ((dat.bz_norm[moornum,dlim:,:]-dat.vz2[moornum,dlim:,:])*1e5).plot(ax=axx[3],vmin=0,vmax=vmaxi,cbar_kwargs={"label":'$b_z (1 + \zeta/f)-v_z^2$'},cmap=PVmap)
#     ylim(1500,0)
#     for xx in axx:
#         xx.set_xlabel('')
#         xx.set_ylabel('depth [m]')
#         xx.set_title('')
#         xx.contour(dat.date.resample(date='1W').mean(dim='date').values,dat.depth.values,dat['potential density'][moornum,:,:].resample(date='1W').mean(dim='date'),levels=[27.65,27.73],colors='k')
#     savefig(figdir+'SI/PV_hovmueller_CF'+str(moornum+1)+'.png',bbox_inches='tight')
#
#     # #plot time series of PV at a few different depths
#     # for dd in [1150]:#[300,624,750,1000]:
#     #     dind=int(dd/2)
#     #     figure(figsize=(20,3))
#     #     plot(dat.date,(dat.bz[moornum,dind,:]-dat.vz2[moornum,dind,:])*1e5,label='$b_z(1+\zeta/f)-v_z^2$',color='red',linewidth=2)
#     #     plot(dat.date,(dat.bz_anom[moornum,dind,:])*1e5,label='$b_z x \zeta /f$',color='orange')
#     #     plot(dat.date,-dat.vz2[moornum,dind,:]*1e5,label='$-v_z^2$',color='blue')
#     #     plot(dat.date,(dat.bz[moornum,dind,:])*1e5,label='$b_z$',color='black')
#     #     legend()
#     #     axhline(0,color='grey')
#     #     ylim(-0.2,0.9)
#     #     xlim([datetime.datetime(2014,8,15),datetime.datetime(2018,8,15)])
#     #     savefig(figdir+'SI/PV_tseries_CF'+str(moornum+1)+'_'+str(dind*2)+'m.png',bbox_inches='tight')
#
# for mm in range(4,7):
#     plot_pv_comps(mm)
# def plot_zeta_f():
#     (dat.zeta_f.mean(dim='date').T).plot()
#     contour(dat.distance.values,dat.depth.values,dat['across track velocity'].mean(dim='date').values.T,colors='k',levels=arange(-0.4,0.8,0.05))
#     ylim(2000,0)
#     ylabel('distance [km]')
#     text(45,-20,'CF5')
#     text(60,-20,'CF6')
#     savefig(figdir+'SI/zeta_f_mean.png',bbox_inches='tight')
#
# plot_zeta_f()
#
# def plot_prof(moor,d1,d2):
#     plot(dat['PV'][moor,:,:].sel(date=slice(d1,d2)).mean(dim='date'),dat.depth)
#     plot(dat['bz'][moor,:,:].sel(date=slice(d1,d2)).mean(dim='date'),dat.depth)
#     plot(dat['bz_norm'][moor,:,:].sel(date=slice(d1,d2)).mean(dim='date'),dat.depth)
#     ylim(750,0)
#     PVlim=5e-5
#     xlim(-PVlim,PVlim)
#
# plot_prof(4,'2016-11-20','2016-11-30')
# plot_prof(4,'2017-01-15','2017-01-15')
# plot_prof(4,'2017-01-20','2017-01-25')
# /plot_prof(6,'2017-01-20','2017-01-25')
#
#
# pvcheck_all=array(pvdates,dtype='datetime64[D]')
#
# def plot_slice(var,var2):
#     dat[var].sel(depth=50).sel(distance=slice(35,80)).plot(vmin=uni[var2]['vmin'],vmax=uni[var2]['vmax'],cmap=uni[var2]['cmap'])
#     savefig(figdir+'SI/DSlice_'+var2+'.png',bbox_inches='tight')
#
# plot_slice('bx','bx')
# plot_slice('salinity','sal')
# plot_slice('temperature','tmp')
#
# def plot_lineden(var):
#     figure(figsize=(20,3))
#     plot(dat.date,dat[var].sel(depth=50).sel(distance=slice(40,80)).T)
#     [axvline(dd,color='grey',alpha=0.3,linewidth=5) for dd in pvcheck_all]
#     [axvline(dd,color='k',alpha=0.3,linewidth=5) for dd in pvcheck]
#     if var=='salinity':
#         axhline(34,color='C0')
#     ylabel(var)
#     xlim(dt.datetime(2014,8,15),dt.datetime(2018,8,15))
#     savefig(figdir+'SI/LineSlice_'+var+'.png',bbox_inches='tight')
#
# plot_lineden('salinity')
#
# for var in ['bx','salinity','temperature','potential density']:
#     plot_lineden(var)
# # Could also find out where the (surface) density at CF6 would stabilize in CF5 column
#
# def plot_TS(datech,axi):
#     axi.plot(dat.sel(date=datech).sel(distance=slice(40,80)).sel(depth=slice(50,2000)).salinity.T,dat.sel(date=datech).sel(distance=slice(40,80)).sel(depth=slice(50,2000)).temperature.T,'.')
#     axi.plot(dat.sel(date=datech).sel(distance=slice(40,80)).sel(depth=200).salinity.T,dat.sel(date=datech).sel(distance=slice(40,80)).sel(depth=200).temperature.T,'k*')
#     axi.contour(salvec,tmpvec,pdenmat2,colors='grey',levels=arange(26.5,28.5,0.1))
#     axi.set_xlim(33.6,35.2)
#     axi.set_ylim(2,6.5)
#
# #make a composite of really fresh times right before pvchecks, or just pick each of them out!
#
# [(2*ii-6) for ii in range(6)]
#
# for pp in pvcheck:
#     f,axx=subplots(1,7,sharex=True,sharey=True,figsize=(24,3))
#     for ii,axi in enumerate(axx):
#         eldato=pp+timedelta64(ii)
#         plot_TS(eldato,axi)
#         axi.set_title(eldato)
#     axx[3].set_xlabel('salinity')
#     axx[0].set_ylabel('pot. temperature')
#     savefig(figdir+'SI/TS_series_'+str(pp)+'.png')
#     for ii,axi in enumerate(axx):
#         axi.contour(salvec,tmpvec,pdenmat2,colors='grey',levels=arange(26.5,28.5,0.01))
#     axx[0].set_xlim()
