from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'

################################################################################################################################
#################
# # # #### load osnap data and cut out the eastern portion
lonbnd=-44
dat=xr.open_dataset(datadir+'OSNAP_Gridded_Fields_201408_201805_2020.nc')
osnap_obs=dat.sel(LONGITUDE=slice(lonbnd,0))
osnap_obs['LATITUDE']=dat['LATITUDE'].values[dat.LONGITUDE>lonbnd]

#add absolute salinity and conservative temperature
def add_SA_CT(osnap_obs):
    PRES_out=gsw.p_from_z(-osnap_obs['DEPTH'],osnap_obs['LATITUDE'].mean())
    SA_out=NaN*osnap_obs['PSAL'].copy()
    for ii,pp in enumerate(PRES_out):
        for jj,ll in enumerate(osnap_obs.LONGITUDE.values):
            SA_out[:,ii,jj]=gsw.SA_from_SP(osnap_obs['PSAL'][:,ii,jj],pp,ll,osnap_obs.LATITUDE[jj])
    if 'THETA' in osnap_obs:
        CT_out=gsw.CT_from_pt(SA_out,osnap_obs['THETA'])
        osnap_obs['PDEN']=osnap_obs['SIGMA_THETA']-1e3
    else:
        CT_out=gsw.CT_from_pt(SA_out,osnap_obs['PTMP'])
    osnap_obs['SA']=(('TIME','DEPTH','LONGITUDE'),SA_out)
    osnap_obs['CT']=(('TIME','DEPTH','LONGITUDE'),CT_out)
    return osnap_obs

osnap_obs=add_SA_CT(osnap_obs)

fs_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-et-al-2020/Tsubouchi2020_fs_xray_2008.nc')#.sel(TIME=slice('2005-09-01','2011-3-31'))
bso_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-et-al-2020/Tsubouchi2020_bso_xray_2008.nc')#.sel(TIME=slice('2005-09-01','2011-3-31'))
fs_obs=fs_obs.transpose('TIME','DEPTH','LONGITUDE','LON_LONG')
bso_obs=bso_obs.transpose('TIME','DEPTH','LONGITUDE','LON_LONG')

####Load NorESM sections
osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_18yrs_2004.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_18yrs_2004.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_18yrs_2004.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_18yrs_2004.nc')


for xrai in [osnap,fs,bso,ns]:
    xrai=add_SA_CT(xrai)



################################################################################################################################
################################################################################################################################
###########################################      GET overturning    #######################################################
################################################################################################################################
################################################################################################################################
denstep=0.01
denvec=arange(25.9,28.3,denstep)
def getpsi(east):
    east['TRANS']=east['VELO']*east['AREA']/1e6
    psimat=NaN*ones((len(denvec),len(east.TIME)))
    for ii,dd in enumerate(denvec):
            psimat[ii,:]=east['TRANS'].where(east['PDEN']>dd-denstep/2).where(east['PDEN']<=dd+denstep/2).sum(dim='DEPTH').sum(dim='LONGITUDE').values
    east=east.assign_coords(DENLEV=(denvec))
    east['TRANSDEN']=((['DENLEV','TIME']),psimat)
    east['TRANSDEN_adj']=(east.TRANSDEN-(east.TRANSDEN.mean(dim='TIME').sum(dim='DENLEV')/len(east.DENLEV)))
    east['PSI']=east['TRANSDEN'].cumsum(dim='DENLEV')
    east['PSI_adj']=east['TRANSDEN_adj'].cumsum(dim='DENLEV')
    east['PSI_inv']=east['TRANSDEN'][::-1,:].cumsum(dim='DENLEV')
    east['SIGMAX']=(('TIME'),east.DENLEV[east.PSI.argmax(dim='DENLEV').values])
    east['SIGMAX_adj']=(('TIME'),east.DENLEV[east.PSI_adj.argmax(dim='DENLEV').values])
    east['SIGMAX_inv']=(('TIME'),east.DENLEV[east.PSI_inv.argmin(dim='DENLEV').values])
    east['SIGMEAN']=east.DENLEV[east.PSI_adj.mean(dim='TIME').argmax(dim='DENLEV').values]
    east['MOC']=east.PSI.max(dim='DENLEV')
    east['MOC_adj']=east.PSI_adj.max(dim='DENLEV')
    east['MOC_inv']=-east.PSI_inv.min(dim='DENLEV')
    east['MOCMEAN']=east.PSI[east.PSI_adj.mean(dim='TIME').argmax(dim='DENLEV').values,:]
    return east

osnap_obs=getpsi(osnap_obs)
bso_obs=getpsi(bso_obs)
fs_obs=getpsi(fs_obs)

#model
osnap=getpsi(osnap)
bso=getpsi(bso)
fs=getpsi(fs)

def plot_comp_psi(which,whichobs,name,psitype='PSI'):
    figure()
    plot(whichobs[psitype].mean(dim='TIME'),whichobs.DENLEV,linewidth=3,label=name+' extended',color='limegreen')
    plot(which[psitype].mean(dim='TIME'),which.DENLEV,linewidth=3,label='NorESM (2010-2018)',color='purple')
    if name=='OSNAP':
        plot(which[psitype].sel(TIME=slice('2014-8-1','2018-5-11')).mean(dim='TIME'),which.DENLEV,linewidth=3,label='NorESM during OSNAP',color='purple',linestyle='--')
    xlabel('Streamfunction [Sv]')
    ylabel('pot. density anomaly [kg m$^{-3}$]')
    legend(fontsize=16,loc=(1.05,0.2))
    ylim(28.3,26)
    savefig(figdir+name+'_PsiComp.png',bbox_inches='tight')
    savefig(figdir+name+'_PsiComp.pdf',bbox_inches='tight')

def quick_MOC_comp():
    figure()
    osnap_obs.SIGMAX.plot()
    osnap_obs.SIGMAX_adj.plot()
    osnap_obs.SIGMAX_inv.plot()

    figure()
    osnap_obs.MOC.plot()
    osnap_obs.MOC_adj.plot()
    osnap_obs.MOC_inv.plot()
    osnap_obs.MOCMEAN.plot()

quick_MOC_comp()


plot_comp_psi(osnap,osnap_obs,'OSNAP')
plot_comp_psi(osnap,osnap_obs,'OSNAP',psitype='PSI_adj')
plot_comp_psi(osnap,osnap_obs,'OSNAP',psitype='PSI_inv')
plot_comp_psi(fs,fs_obs,'Fram Strait')
plot_comp_psi(bso,bso_obs,'Barents Sea Opening')
#
# def MOC_comp():
#     osnap.MOC.plot()
#     osnap.MOCMEAN.plot()
#     osnap_obs.MOC.plot()
#     osnap_obs.MOCMEAN.plot()
#     xlim(datetime.datetime(2014,1,1),datetime.datetime(2017,1,1))
#     figure()
#     osnap.SIGMAX.plot()
#     axhline(osnap.SIGMEAN)
#     osnap_obs.SIGMAX.plot()
#     axhline(osnap_obs.SIGMEAN,color='C1')
#     xlim(datetime.datetime(2014,1,1),datetime.datetime(2017,1,1))
#     figure()
#     osnap.TRANS.sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
#     osnap_obs.TRANS.sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
#
# MOC_comp()
osnap=getpsi(osnap)
bso=getpsi(bso)
fs=getpsi(fs)

################################################################################################################################
################################################################################################################################
########################################      PLOT MEAN SECTION PANELS    #######################################################
################################################################################################################################
################################################################################################################################
sigmax=osnap.SIGMEAN.values
sigmax_obs=27.56#osnap_obs.SIGMEAN.values GO WITH FEILI's value

sigmax
sigmax_obs

oslim=-38
oslim_obs=-40

fslim=2
fslim_obs=-2

vmax=0.3
univec['VELO']=['across track velocity',arange(-vmax,vmax+0.05,0.05),cm.RdBu_r,arange(-vmax,vmax+0.1,0.1),'[m/s]']
univec['PTMP']=['pot. temperature', arange(-1,11,0.5),cm.RdYlBu_r,range(0,11,2),'[$^\\circ$C]']
univec['PSAL']=['salinity',arange(34,35.6,0.1),cm.PiYG_r, array([34., 34.5 ,35,35.5,]),'']

prop_key=['VELO','PSAL','PTMP']
label_key=['velocity [m s$^{-1}$]','salinity','pot. temperature [$^\circ$C]']

titvec=['OSNAP East','Fram Strait','Barents Sea Opening']
obsvec=[osnap_obs,fs_obs,bso_obs]
norvec=[osnap,fs,bso]
yspl=[100,400,50]
ymax=[3200,3000,500]

# VELOCITY, SALINITY, TEMPERATURE AND DENSITY MEANS OVER PERIOD OF OVERLAP
def plot_each(dat,xxvar,axx,var,ymax,sigi):
    filled=axx.contourf(xxvar,dat.DEPTH,dat[var].groupby('TIME.month').mean('TIME').mean(dim='month'),univec[var][1],cmap=univec[var][2],extend='both')
    axx.contour(xxvar,dat.DEPTH,dat[var].groupby('TIME.month').mean('TIME').mean(dim='month'),levels=univec[var][1][::2],colors='k')
    axx.contour(xxvar,dat.DEPTH,dat['PDEN'].groupby('TIME.month').mean('TIME').mean(dim='month'),levels=[sigi],colors='k',linewidths=4,zorder=100) # add isopycnal of maximum overturning
    axx.set_facecolor('k')
    axx.set_ylim(ymax,0)
    return filled

def cont_partial(dat,xxvar,axx,ymin,ymax,sigi,var='VELO'):
    filled=axx.contourf(xxvar,dat.DEPTH,dat[var].groupby('TIME.month').mean('TIME').mean(dim='month'),arange(-0.3,0.3,0.025),cmap=univec[var][2],extend='both')
    axx.contour(xxvar,dat.DEPTH,dat['PDEN'].mean(dim='TIME'),levels=[sigi],colors='k',linewidths=4) # add isopycnal of maximum overturning
    axx.set_facecolor('k')
    axx.set_ylim(ymax,ymin)
    return filled#,conts

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}

fsdplim=250

def sections_obs_only():
    fig = plt.figure(figsize=(12, 9), constrained_layout=False)
    grd_outer = fig.add_gridspec(2,2, hspace=0.4)
    ax1=grd_outer[0,:]
    ax2=grd_outer[1,0]
    ax3=grd_outer[1,1]
    fsz=16
    for i,axi in enumerate([ax1,ax2,ax3]):
        tit=titvec[i]
        obs=obsvec[i]
        grd = axi.subgridspec(2,1, wspace=0.1, hspace=0.0,height_ratios=[1,5])
        axx={}
        axx[0,0]=fig.add_subplot(grd[0])
        axx[1,0]=fig.add_subplot(grd[1])
        if 'Barents' in tit:
            obs_x=obs.LATITUDE
        else:
            obs_x=obs.LONGITUDE
        velcont=cont_partial(obs,obs_x,axx[1,0],yspl[i],ymax[i],sigmax_obs)
        cont_partial(obs,obs_x,axx[0,0],0,yspl[i],sigmax_obs)
        if 'OSNAP' in tit:
            axx[0,0].plot([oslim_obs,oslim_obs],[0,50],color='k',linewidth=3)
            axx[0,0].text(-15,80,'AWS',color='k',fontsize=fsz,weight='bold')
            axx[0,0].text(-43,-10,'PWS',color='k',fontsize=fsz,weight='bold')
            axx[1,0].text(-40,1500,'DWS',color='k',fontsize=fsz,weight='bold')
            axx[1,0].set_xlabel('Longitude [$^\circ$E]',fontsize=fsz)
            axx[1,0].set_ylabel('depth [m]',fontsize=fsz)
            axx[1,0].text(-19,2750,tit,color='white',fontsize=fsz+5)
        elif 'Fram' in tit:
            axx[1,0].text(-18,2600,tit,color='white',fontsize=fsz+5)
            axx[1,0].text(0,700,'AWN',color='k',fontsize=fsz,weight='bold')
            axx[1,0].text(-8,-200,'PWN',color='k',fontsize=fsz,weight='bold')
            axx[1,0].set_xlabel('Longitude [$^\circ$E]',fontsize=fsz)
            axx[1,0].set_ylabel('depth [m]',fontsize=fsz)
            axx[0,0].contour(obs_x,fs_obs.DEPTH,fs_obs['PTMP'].where(~isnan(fs_obs['VELO'])).mean(dim='TIME'),levels=[0],colors='r',linewidths=4)
            # axx[0,0].axvline(-5)
            # axx[0,0].plot([fslim_obs,fslim_obs],[0,fsdplim],color='purple',linewidth=3)
            # axx[0,0].plot([-20,fslim_obs],[fsdplim,fsdplim],color='purple',linewidth=3)
            # axx[0,0].plot([-20,-20],[0,fsdplim],color='purple',linewidth=3)

        else:
            axx[1,0].text(74,450,'Barents Sea\nOpening',color='white',fontsize=fsz+5)
            axx[1,0].text(72,150,'AWN',color='k',fontsize=fsz,weight='bold')
            axx[1,0].set_xlabel('Latitude [$^\circ$N]',fontsize=fsz)
    xc=0.95
    cwi=0.02
    cle=0.4
    caxit_vel=fig.add_axes([xc,0.3,cwi,cle])
    colorbar(velcont,label='velocity [m/s]',cax=caxit_vel,ticks=arange(-0.3,0.3,0.1))
    savefig(figdir_paper+'Sections_obs_only.png',bbox_inches='tight')
    savefig(figdir_paper+'Sections_obs_only.pdf',bbox_inches='tight')

sections_obs_only()



# def sections_model_only():
#     fig = plt.figure(figsize=(12,9), constrained_layout=False)
#     grd_outer = fig.add_gridspec(2,2, hspace=0.4)
#     ax1=grd_outer[0,:]
#     ax2=grd_outer[1,0]
#     ax3=grd_outer[1,1]
#     fsz=16
#     for i,axi in enumerate([ax1,ax2,ax3]):
#         tit=titvec[i]
#         obs=norvec[i]
#         grd = axi.subgridspec(2,1, wspace=0.1, hspace=0.0,height_ratios=[1,5])
#         axx={}
#         axx[0,0]=fig.add_subplot(grd[0])
#         axx[1,0]=fig.add_subplot(grd[1])
#         if 'Barents' in tit:
#             obs_x=obs.LATITUDE
#         else:
#             obs_x=obs.LONGITUDE
#         velcont=cont_partial(obs,obs_x,axx[1,0],yspl[i],ymax[i],sigmax)
#         cont_partial(obs,obs_x,axx[0,0],0,yspl[i],sigmax)
#         if 'OSNAP' in tit:
#             axx[0,0].plot([oslim,oslim],[0,75],color='k',linewidth=3)
#             axx[0,0].text(-15,80,'AWS',color='k',fontsize=fsz,weight='bold')
#             axx[0,0].text(-43,-10,'PWS',color='k',fontsize=fsz,weight='bold')
#             axx[1,0].text(-40,1500,'DWS',color='k',fontsize=fsz,weight='bold')
#             axx[1,0].set_xlabel('Longitude [$^\circ$E]',fontsize=fsz)
#             axx[1,0].set_ylabel('depth [m]',fontsize=fsz)
#             axx[1,0].text(-19,2750,tit,color='white',fontsize=fsz+5)
#         elif 'Fram' in tit:
#             # axx[0,0].plot([fslim,fslim],[0,fsdplim],color='purple',linewidth=3)
#             # axx[0,0].plot([-20,fslim],[fsdplim,fsdplim],color='purple',linewidth=3)
#             # axx[0,0].plot([-20,-20],[0,fsdplim],color='purple',linewidth=3)
#             axx[1,0].text(-16,2600,tit,color='white',fontsize=fsz+5)
#             axx[1,0].text(0,700,'AWN',color='k',fontsize=fsz,weight='bold')
#             axx[1,0].text(-8,-150,'PWN',color='k',fontsize=fsz,weight='bold')
#             axx[1,0].set_xlabel('Longitude [$^\circ$E]',fontsize=fsz)
#             axx[1,0].set_ylabel('depth [m]',fontsize=fsz)
#             axx[0,0].contour(obs_x,fs.DEPTH,fs['PTMP'].mean(dim='TIME'),levels=[0],colors='r',linewidths=4)
#         else:
#             axx[1,0].text(73.5,450,'Barents Sea\nOpening',color='white',fontsize=fsz+5)
#             axx[1,0].text(72,150,'AWN',color='k',fontsize=fsz,weight='bold')
#             axx[1,0].set_xlabel('Latitude [$^\circ$N]',fontsize=fsz)
#     xc=0.95
#     cwi=0.02
#     cle=0.4
#     caxit_vel=fig.add_axes([xc,0.3,cwi,cle])
#     # caxit_sal=fig.add_axes([xc,0.125,cwi,cle])
#     colorbar(velcont,label='velocity [m/s]',cax=caxit_vel,ticks=arange(-0.3,0.3,0.1))
#     # colorbar(salcont,label='Salinity',cax=caxit_sal)
#     savefig(figdir_paper+'Sections_model_only.png',bbox_inches='tight')
#     savefig(figdir_paper+'Sections_model_only.pdf',bbox_inches='tight')
#
# sections_model_only()


############################################################################################
################  NORESM  WATER MASS PARTITIONING   #############################################

AWS={}
PWS={}
DWS={}
AWN={}
PWN={}

pws_depth=100

xray=osnap
for var in ['TRANS','SA','CT','PDEN']:
    if var=='TRANS':
        DWS[var]=xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS[var]=xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS[var]=xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')

    else:
        DWS[var]=(xray[var]*xray['TRANS']).where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')

xray=fs
for var in ['TRANS','SA','CT','PDEN']:
    if var=='TRANS':
        PWN[var]=-fs['TRANS'].where(fs.CT<0).where(fs.DEPTH<500).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN[var]=-fs['TRANS'].where(fs.CT>=0).sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].where(fs.PTMP<0).where(fs.DEPTH>=500).sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
        # PWN[var]=-fs['TRANS'].where(fs.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        # AWN[var]=-fs['TRANS'].where(fs.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
    else:
        PWN[var]=-(fs[var]*fs['TRANS']).where(fs.CT<0).where(fs.DEPTH<500).sum(dim='DEPTH').sum(dim='LONGITUDE')/PWN['TRANS']
        AWN[var]=-((fs[var]*fs['TRANS']).where(fs.CT>=0).sum(dim='DEPTH').sum(dim='LONGITUDE')+(fs[var]*fs['TRANS']).where(fs.PTMP<0).where(fs.DEPTH>=500).sum(dim='DEPTH').sum(dim='LONGITUDE')+(bso[var]*bso['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE'))/AWN['TRANS']
        # PWN[var]=-(fs[var]*fs['TRANS']).where(fs.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/PWN['TRANS']
        # AWN[var]=-((fs[var]*fs['TRANS']).where(fs.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')+(bso[var]*bso['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE'))/AWN['TRANS']

# for var in ['PSAL','PTMP','PDEN','TRANS']:
#     if var=='TRANS':
#         PWN[var]=-fs['TRANS'].where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')
#         AWN[var]=-fs['TRANS'].where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
#     else:
#         PWN[var]=(fs[var]*fs['TRANS']).where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')/fs['TRANS'].where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')
#         AWN[var]=((fs[var]*fs['TRANS']).where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')+(bso[var]*bso['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE'))/(fs['TRANS'].where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')+bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE'))

def plot_WMN():
    scatter(PWN['SA'].groupby('TIME.year').mean('TIME').values,PWN['CT'].groupby('TIME.year').mean('TIME').values,color='purple')
    scatter(AWN['SA'].groupby('TIME.year').mean('TIME').values,AWN['CT'].groupby('TIME.year').mean('TIME').values,color='orange')
    # xlim(34,35.5)
    # ylim(-2,12)

plot_WMN()


AWN['CT'].groupby('TIME.year').mean('TIME').plot()

############################################################################################
################  OBS WATER MASS PARTITIONING   #############################################

def get_obs_WMS(trans_var):
    AWS_obs={}
    PWS_obs={}
    DWS_obs={}
    AWN_obs={}
    AWN_obs_fs={}
    AWN_obs_bso={}
    PWN_obs={}

    xray=osnap_obs
    sigmax_obs=osnap_obs.SIGMEAN.values
    for var in ['SA','CT','PDEN','TRANS']:
        if var=='TRANS':
            DWS_obs[var]=xray[trans_var].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            PWS_obs[var]=xray[trans_var].where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWS_obs[var]=xray[trans_var].where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        else:
            DWS_obs[var]=(xray[var]*xray[trans_var]).where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray[trans_var].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            PWS_obs[var]=(xray[var]*xray[trans_var]).where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray[trans_var].where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWS_obs[var]=(xray[var]*xray[trans_var]).where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray[trans_var].where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')

    xray=fs_obs
    for var in ['TRANS','SA','CT','PDEN',]:
        if var=='TRANS':
            PWN_obs[var]=-fs_obs[trans_var].where(xray.CT<0).where(xray.DEPTH<500).sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWN_obs_fs[var]=-fs_obs[trans_var].where(xray.CT>=0).sum(dim='DEPTH').sum(dim='LONGITUDE')-fs_obs[trans_var].where(xray.CT<0).where(xray.DEPTH>=500).sum(dim='DEPTH').sum(dim='LONGITUDE')
            # PWN_obs[var]=-fs_obs[trans_var].where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            # AWN_obs_fs[var]=-fs_obs[trans_var].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWN_obs_bso[var]=-bso_obs[trans_var].sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWN_obs[var]=AWN_obs_fs[var]+AWN_obs_bso[var]
        else:
            PWN_obs[var]=-(fs_obs[var]*fs_obs[trans_var]).where(xray.CT<0).where(xray.DEPTH<500).sum(dim='DEPTH').sum(dim='LONGITUDE')/PWN_obs['TRANS']
            AWN_obs_fs[var]=-((fs_obs[var]*fs_obs[trans_var]).where(xray.CT>=0).sum(dim='DEPTH').sum(dim='LONGITUDE')+(fs_obs[var]*fs_obs[trans_var]).where(xray.CT<0).where(xray.DEPTH>500).sum(dim='DEPTH').sum(dim='LONGITUDE'))/AWN_obs_fs['TRANS']
            # PWN_obs[var]=-(fs_obs[var]*fs_obs[trans_var]).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/PWN_obs['TRANS']
            # AWN_obs_fs[var]=(fs_obs[var]*fs_obs[trans_var]).where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/fs_obs[trans_var].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWN_obs_bso[var]=(bso_obs[var]*bso_obs[trans_var]).sum(dim='DEPTH').sum(dim='LONGITUDE')/bso_obs[trans_var].sum(dim='DEPTH').sum(dim='LONGITUDE')
            AWN_obs[var]=(AWN_obs_fs[var]*AWN_obs_fs['TRANS']+AWN_obs_bso[var]*AWN_obs_bso['TRANS'])/AWN_obs['TRANS']

    return AWS_obs,PWS_obs,DWS_obs,PWN_obs,AWN_obs,AWN_obs_fs,AWN_obs_bso

AWS_obs,PWS_obs,DWS_obs,PWN_obs,AWN_obs,AWN_obs_fs,AWN_obs_bso=get_obs_WMS('TRANS')

axvline(-5)
axvline(-6)
axhline(-0.6)
fs_obs['TRANS'].where(fs_obs.CT<0).where(fs_obs.DEPTH<500).sum(dim='DEPTH').mean(dim='TIME').cumsum(dim='LONGITUDE').plot()

def plot_WMN():
    scatter(PWN_obs['SA'].groupby('TIME.year').mean('TIME').values,PWN_obs['CT'].groupby('TIME.year').mean('TIME').values,color='purple')
    scatter(AWN_obs['SA'].groupby('TIME.year').mean('TIME').values,AWN_obs['CT'].groupby('TIME.year').mean('TIME').values,color='orange')


plot_WMN()

####################################################################################
#########  OBS WATER MASS PARTITIONING with MASS BALANCE   #########################
####################################################################################
# apply mass balance at osnap
trans_tot=osnap_obs.TRANS.sum(dim='LONGITUDE').sum(dim='DEPTH')
vel_corr=trans_tot/osnap_obs.AREA.sum(dim='LONGITUDE').sum(dim='DEPTH')
osnap_obs['TRANS_adj']=osnap_obs['TRANS']-vel_corr*osnap_obs.AREA

# apply mass balance at fs + bso
trans_tot_north=fs_obs.TRANS.sum(dim='LONGITUDE').sum(dim='DEPTH')+bso_obs.TRANS.sum(dim='LONGITUDE').sum(dim='DEPTH')
vel_corr_north=trans_tot_north/(fs_obs.AREA.sum(dim='LONGITUDE').sum(dim='DEPTH')+bso_obs.AREA.sum(dim='LONGITUDE').sum(dim='DEPTH'))
(fs_obs.AREA.sum(dim='LONGITUDE').sum(dim='DEPTH')+bso_obs.AREA.sum(dim='LONGITUDE').sum(dim='DEPTH'))

fs_obs['TRANS_adj']=fs_obs['TRANS']-vel_corr_north*fs_obs.AREA
bso_obs['TRANS_adj']=bso_obs['TRANS']-vel_corr_north*bso_obs.AREA

fs_obs['TRANS'].sum(dim='LONGITUDE').sum(dim='DEPTH').plot()
fs_obs['TRANS_adj'].sum(dim='LONGITUDE').sum(dim='DEPTH').plot()
bso_obs['TRANS'].sum(dim='LONGITUDE').sum(dim='DEPTH').plot()
bso_obs['TRANS_adj'].sum(dim='LONGITUDE').sum(dim='DEPTH').plot()

AWS_obs_mb,PWS_obs_mb,DWS_obs_mb,PWN_obs_mb,AWN_obs_mb,AWN_obs_fs_mb,AWN_obs_bso_mb=get_obs_WMS('TRANS_adj')


# ################################################################################################################################
# ################################################################################################################################
# ###########################################      SAVE WM props    #######################################################
# ################################################################################################################################
# ################################################################################################################################

WM_obs={}
for ii,xrw in enumerate([AWS_obs,PWS_obs,DWS_obs,AWN_obs,PWN_obs]):
    WM_obs[ii]=xr.concat([xrw['PDEN'],xrw['SA'],xrw['CT'],xrw['TRANS']],pd.Index(['PDEN','SA','CT','TRANS'],name='var'))

WMall_obs=xr.concat([WM_obs[ww] for ww in WM_obs],pd.Index(['AWS','PWS','DWS','AWN','PWN'],name='WM'))

WMall_obs=WMall_obs.to_dataset('var')
WMall_obs.to_netcdf(datadir+'FW_WM/OSNAP2014-18_Tsub2020_WM_2008.nc','w')


WM={}
for ii,xrw in enumerate([AWS,PWS,DWS,AWN,PWN]):
    WM[ii]=xr.concat([xrw['PDEN'],xrw['SA'],xrw['CT'],xrw['TRANS']],pd.Index(['PDEN','SA','CT','TRANS'],name='var'))

WMall=xr.concat([WM[ww] for ww in WM],pd.Index(['AWS','PWS','DWS','AWN','PWN'],name='WM'))
WMall=WMall.to_dataset('var')

WMall.to_netcdf(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc','w')

WM_obs_mb={}
for ii,xrw in enumerate([AWS_obs_mb,PWS_obs_mb,DWS_obs_mb,AWN_obs_mb,PWN_obs_mb]):
    WM_obs_mb[ii]=xr.concat([xrw['PDEN'],xrw['SA'],xrw['CT'],xrw['TRANS']],pd.Index(['PDEN','SA','CT','TRANS'],name='var'))
WMall_obs_mb=xr.concat([WM_obs_mb[ww] for ww in WM_obs_mb],pd.Index(['AWS','PWS','DWS','AWN','PWN'],name='WM'))
WMall_obs_mb=WMall_obs_mb.to_dataset('var')
WMall_obs_mb.to_netcdf(datadir+'FW_WM/OSNAP2014-18_WM_mb_2008.nc','w')

### GRAVEYARD
##plot every time step of FS and BSO vel, including tmp and den boundaries

# def plot_all_sects(dat,xxvar,ymin,ymax,sigi,var='VELO',fs=1):
#     for ii in range(len(dat.TIME)):
#         f,axx=subplots(1,1)
#         # filled=axx.contourf(xxvar,dat.DEPTH,dat[var].isel(TIME=ii),arange(-0.3,0.3,0.025),cmap=univec[var][2],extend='both')
#         filled=axx.contourf(xxvar,dat.DEPTH,dat[var].isel(TIME=ii),arange(-2,10,0.5),extend='both')
#         # axx.contour(xxvar,dat.DEPTH,dat['PDEN'].isel(TIME=ii),levels=[sigi],colors='k',linewidths=4) # add isopycnal of maximum overturning
#         # axx.contour(xxvar,dat.DEPTH,dat['CT'].isel(TIME=ii),levels=[0],colors='r',linewidths=4) # add 0degree isotherm CT
#         # axx.contour(xxvar,dat.DEPTH,dat['PTMP'].isel(TIME=ii),levels=[0],colors='orange',linblueewidths=4) # add 0degree isotherm PTMP
#         axx.set_facecolor('k')
#         axx.set_ylim(ymax,ymin)
#         axx.set_title(AWN_obs_fs['CT'].isel(TIME=ii).values)
#         colorbar(filled)
#         if fs==1:
#             savefig('/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/FS_sects/FS_CT_'+str(ii)+'.png')
#         elif fs==0:
#             savefig('/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/BSO_sects/BSO_CT_'+str(ii)+'.png')
#
#
# plot_all_sects(fs_obs,fs_obs.LONGITUDE,0,3000,sigmax_obs,var='CT')
#
# plot_all_sects(bso_obs,bso_obs.LATITUDE,0,500,sigmax_obs,fs=0,var='CT')
#
# plot_all_sects(fs_obs,fs_obs.LONGITUDE,0,3000,sigmax_obs)
#
# plot_all_sects(bso_obs,bso_obs.LATITUDE,0,500,sigmax_obs,fs=0)
# def calc_cumsums(sec):
#     sec['xp_xcum']=(sec.TRANS).mean(dim='TIME').sum(dim='DEPTH').cumsum(dim='LONGITUDE')
#     sec['xp_zcum']=(sec.TRANS).mean(dim='TIME').sum(dim='LONGITUDE').cumsum(dim='DEPTH')
#     sec['Txp_xcum']=(sec.TRANS*sec.CT).mean(dim='TIME').sum(dim='DEPTH').cumsum(dim='LONGITUDE')
#     sec['Txp_zcum']=(sec.TRANS*sec.CT).mean(dim='TIME').sum(dim='LONGITUDE').cumsum(dim='DEPTH')
#     return sec
#
# fs_obs=calc_cumsums(fs_obs)
# bso_obs=calc_cumsums(bso_obs)
#
# def plot_transcum(sec):
#     f,axx=subplots(3,2,figsize=(10,8))
#     sec['xp_xcum'].plot(ax=axx[0,0])
#     sec['Txp_xcum'].plot(ax=axx[1,0])
#     (sec['Txp_xcum']/sec['xp_xcum']).plot(ax=axx[2,0])
#     sec['xp_zcum'].plot(ax=axx[0,1])
#     sec['Txp_zcum'].plot(ax=axx[1,1])
#     (sec['Txp_zcum']/sec['xp_zcum']).plot(ax=axx[2,1])
#     suptitle((sec['Txp_xcum']/sec['xp_xcum'])[-1].values)
#     for axvec in axx:
#         for axi in axvec:
#             axi.axhline(0,color='k')
#
# plot_transcum(bso_obs)
# plot_transcum(fs_obs)
