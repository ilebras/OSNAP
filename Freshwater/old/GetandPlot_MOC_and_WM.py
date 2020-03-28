from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'
import cmocean

################################################################################################################################
###########################################    Load NorESM and obs    #######################################################
osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')

(bso['VELO']*bso['AREA']).sum(dim='DEPTH').sum(dim='LONGITUDE').mean(dim='TIME')/1e6

# # #### load osnap data and cut out the eastern portion
lonbnd=-44
dat=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.nc')
osnap_obs=dat.sel(LONGITUDE=slice(lonbnd,0))
osnap_obs['LATITUDE']=dat['LATITUDE'].values[dat.LONGITUDE>lonbnd]

fs_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_fs_xray_2001.nc')
bso_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_bso_xray_2001.nc')
fs_obs=fs_obs.transpose('TIME','DEPTH','LONGITUDE','LON_LONG')
bso_obs=bso_obs.transpose('TIME','DEPTH','LONGITUDE','LON_LONG')


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
    east['PSI']=east['TRANSDEN'].cumsum(dim='DENLEV')
    east['SIGMAX']=(('TIME'),east.DENLEV[east.PSI.argmax(dim='DENLEV').values])
    east['SIGMEAN']=east.DENLEV[east.PSI.mean(dim='TIME').argmax(dim='DENLEV').values]
    east['MOC']=east.PSI.max(dim='DENLEV')
    east['MOCMEAN']=east.PSI[east.PSI.mean(dim='TIME').argmax(dim='DENLEV').values,:]

    return east

osnap=getpsi(osnap)
osnap_obs=getpsi(osnap_obs)
bso=getpsi(bso)
bso_obs=getpsi(bso_obs)
fs=getpsi(fs)
fs_obs=getpsi(fs_obs)

def plot_comp_psi(which,whichobs,name):
    figure()
    # plot(osnap_obs.PSI,osnap_obs.DENLEV,'C0',alpha=0.2,label='')
    # plot(osnap.PSI,osnap.DENLEV,'C1',alpha=0.2,label='')
    plot(whichobs.PSI.mean(dim='TIME'),whichobs.DENLEV,linewidth=3,label=name+' (2014-2016)',color='limegreen')
    plot(which.PSI.mean(dim='TIME'),which.DENLEV,linewidth=3,label='NorESM (2010-2018)',color='purple')
    if name=='OSNAP':
        plot(which.PSI.sel(TIME=slice('2014-8-1','2016-9-1')).mean(dim='TIME'),which.DENLEV,linewidth=3,label='NorESM during OSNAP',color='purple',linestyle='--')
    xlabel('Streamfunction [Sv]')
    ylabel('pot. density anomaly [kg m$^{-3}$]')
    legend(fontsize=16)
    ylim(28.3,26)
    savefig(figdir+name+'_PsiComp.png',bbox_inches='tight')
    savefig(figdir+name+'_PsiComp.pdf',bbox_inches='tight')

plot_comp_psi(osnap,osnap_obs,'OSNAP')
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




################################################################################################################################
################################################################################################################################
########################################      PLOT MEAN SECTION PANELS    #######################################################
################################################################################################################################
################################################################################################################################
sigmax=osnap.SIGMEAN.values
sigmax_obs=osnap_obs.SIGMEAN.values


oslim=-39
oslim_obs=-41

fslim=3
fslim_obs=5

vmax=0.3
univec['VELO']=['across track velocity',arange(-vmax,vmax+0.05,0.05),cm.RdBu_r,arange(-vmax,vmax+0.1,0.1),'[m/s]']
univec['PTMP']=['pot. temperature', arange(-1,11,0.5),cm.RdYlBu_r,range(0,11,2),'[$^\\circ$C]']
univec['PSAL']=['salinity',arange(34,35.6,0.1),cm.PiYG_r, array([34., 34.5 ,35,35.5,]),'']

prop_key=['VELO','PSAL','PTMP']
label_key=['velocity [m s$^{-1}$]','salinity','pot. temperature [$^\circ$C]']

# VELOCITY, SALINITY, TEMPERATURE AND DENSITY MEANS OVER PERIOD OF OVERLAP
def plot_each(dat,xxvar,axx,var,ymax,sigi):
    filled=axx.contourf(xxvar,dat.DEPTH,dat[var].groupby('TIME.month').mean('TIME').mean(dim='month'),univec[var][1],cmap=univec[var][2],extend='both')
    axx.contour(xxvar,dat.DEPTH,dat[var].groupby('TIME.month').mean('TIME').mean(dim='month'),levels=univec[var][1][::2],colors='k')
    axx.contour(xxvar,dat.DEPTH,dat['PDEN'].groupby('TIME.month').mean('TIME').mean(dim='month'),levels=[sigi],colors='k',linewidths=4) # add isopycnal of maximum overturning
    axx.set_facecolor('k')
    axx.set_ylim(ymax,0)
    return filled


# def comp_obs_model(nor,obs,ymax,tit,xlab,stit):
#     f,axx=subplots(3,2,figsize=(12,10),sharex=True,sharey=True)
#     cout={}
#     for ii in range(3):
#         if 'Barents' in tit:
#             nor_x=nor.LATITUDE
#             obs_x=obs.LATITUDE
#         else:
#             nor_x=nor.LONGITUDE
#             obs_x=obs.LONGITUDE
#         plot_each(obs,obs_x,axx[ii,0],prop_key[ii],ymax,sigmax_obs)
#         cout[ii]=plot_each(nor,nor_x,axx[ii,1],prop_key[ii],ymax,sigmax)
#     if 'Fram' in tit:
#         for ii in range(3):
#             axx[ii,0].axvline(fslim_obs,color='k',linewidth=3)
#             axx[ii,1].axvline(fslim,color='k',linewidth=3)
#     fs=14
#     axx[0,0].set_title('Observations',fontsize=fs)
#     axx[0,1].set_title('NorESM model',fontsize=fs)
#     f.text(0.05, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fs)
#     f.text(0.5, 0.05, xlab, ha='center',fontsize=fs)
#     suptitle(tit,fontsize=fs+2)
#     caxit={}
#     xc=0.95
#     cwi=0.01
#     cle=0.2
#     caxit[0]=f.add_axes([xc,0.67,cwi,cle])
#     caxit[1]=f.add_axes([xc,0.4,cwi,cle])
#     caxit[2]=f.add_axes([xc,0.13,cwi,cle])
#     for ii in range(3):
#         colorbar(cout[ii],label=label_key[ii],cax=caxit[ii])
#     savefig(figdir+'CompSect_Obs_NorESM_'+stit+'.png',bbox_inches='tight')
#
# comp_obs_model(fs,fs_obs,3000,'Fram Strait','Longitude [$^\circ$W]','fs')
#
# comp_obs_model(bso,bso_obs,500,'Barents Sea Opening','Latitude [$^\circ$N]','bso')
#
# comp_obs_model(osnap,osnap_obs,3200,'OSNAP East','Longitude [$^\circ$W]','osnap')
#
# def cont_partial(dat,xxvar,axx,ymin,ymax,sigi):
#     var='VELO'
#     filled=axx.contourf(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),arange(-0.2,0.2,0.025),cmap=univec[var][2],extend='both')
#     axx.contour(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),colors='k',levels=[0])#univec[var][1][::2])
#     var='PSAL'
#     conts=axx.contour(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),hstack((arange(33,34.5,0.5),arange(34.5,35.8,0.1))),cmap=cmocean.cm.haline,extend='both',linewidths=3)
#
#     axx.contour(xxvar,dat.DEPTH,dat['PDEN'].mean(dim='TIME'),levels=[sigi],colors='k',linewidths=4) # add isopycnal of maximum overturning
#     axx.set_facecolor('k')
#     axx.set_ylim(ymax,ymin)
#     return filled,conts
#
#
# def comp_obs_model_new(nor,obs,ymax,tit,xlab,stit,yspl):
#     fig = plt.figure(figsize=(14, 5), constrained_layout=False)
#     grd = fig.add_gridspec(2,2, wspace=0.1, hspace=0.0,height_ratios=[1,5])
#     fsz=14
#     axx={}
#     axx[0,0]=fig.add_subplot(grd[0,0])
#     axx[0,1]=fig.add_subplot(grd[0,1])
#     axx[1,0]=fig.add_subplot(grd[1,0],sharex=axx[0,0])
#     axx[1,1]=fig.add_subplot(grd[1,1],sharex=axx[0,1])
#     if 'Barents' in tit:
#         nor_x=nor.LATITUDE
#         obs_x=obs.LATITUDE
#     else:
#         nor_x=nor.LONGITUDE
#         obs_x=obs.LONGITUDE
#     velcont,salcont=cont_partial(obs,obs_x,axx[1,0],yspl,ymax,sigmax_obs)
#     cont_partial(obs,obs_x,axx[0,0],0,yspl,sigmax_obs)
#     cont_partial(nor,nor_x,axx[1,1],yspl,ymax,sigmax_obs)
#     cont_partial(nor,nor_x,axx[0,1],0,yspl,sigmax_obs)
#     if 'OSNAP' in tit:
#         axx[0,0].plot([oslim_obs,oslim_obs],[0,50],color='k',linewidth=3)
#         axx[0,1].plot([oslim,oslim],[0,30],color='k',linewidth=3)
#         axx[1,0].text(-19,2750,tit,color='white',fontsize=fsz)
#         axx[0,0].text(-15,80,'AWS',color='k',backgroundcolor='w',fontsize=fsz-2)
#         axx[0,1].text(-15,80,'AWS',color='k',backgroundcolor='w',fontsize=fsz-2)
#         axx[0,0].text(-45,-10,'PWS',color='k',fontsize=fsz-2)
#         axx[0,1].text(-43,-10,'PWS',color='k',fontsize=fsz-2)
#         axx[1,0].text(-40,1500,'DWS',color='k',backgroundcolor='w',fontsize=fsz-2)
#         axx[1,1].text(-40,1500,'DWS',color='k',backgroundcolor='w',fontsize=fsz-2)
#     elif 'Fram' in tit:
#         for ii in range(2):
#             axx[ii,0].axvline(fslim_obs,color='k',linewidth=3)
#             axx[ii,1].axvline(fslim,color='k',linewidth=3)
#             axx[1,0].text(8,1500,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
#             axx[1,1].text(9,1500,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
#             axx[1,0].text(-3,1500,'PWN',color='k',backgroundcolor='w',fontsize=fsz-2)
#             axx[1,1].text(-2,1500,'PWN',color='k',backgroundcolor='w',fontsize=fsz-2)
#     else:
#             axx[1,0].text(72,150,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
#             axx[1,1].text(72,150,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
#     axx[0,1].set_yticklabels('')
#     axx[1,1].set_yticklabels('')
#     axx[0,0].set_title('Observations',fontsize=fsz)
#     axx[0,1].set_title('NorESM model',fontsize=fsz)
#     fig.text(0.05, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fsz)
#     fig.text(0.5, 0, xlab, ha='center',fontsize=fsz)
#     suptitle(tit,fontsize=fsz+2)
#     cwi=0.025
#     cle=0.7
#     caxit_vel=fig.add_axes([0.95,0.15,cwi,cle])
#     caxit_sal=fig.add_axes([1.08,0.15,cwi,cle])
#     colorbar(velcont,label='velocity [m/s]',cax=caxit_vel)
#     colorbar(salcont,label='Salinity',cax=caxit_sal)
#     savefig(figdir+'CompSect_Obs_NorESM_synth_'+stit+'.png',bbox_inches='tight')
#     savefig(figdir+'CompSect_Obs_NorESM_synth_'+stit+'.pdf',bbox_inches='tight')
#
#
#
# comp_obs_model_new(osnap,osnap_obs,3200,'OSNAP East','Longitude [$^\circ$W]','osnap',100)
#
# comp_obs_model_new(fs,fs_obs,3000,'Fram Strait','Longitude [$^\circ$W]','fs',300)
#
# comp_obs_model_new(bso,bso_obs,500,'Barents Sea Opening','Latitude [$^\circ$N]','bso',50)

titvec=['OSNAP East','Fram Strait','Barents Sea Opening']
obsvec=[osnap_obs,fs_obs,bso_obs]
norvec=[osnap,fs,bso]
yspl=[100,350,50]
ymax=[3200,3000,500]

def comp_obs_model_all():
    fig = plt.figure(figsize=(12, 15), constrained_layout=False)
    grd_outer = fig.add_gridspec(3,1, hspace=0.4)
    fsz=16
    for i in range(3):
        tit=titvec[i]
        obs=obsvec[i]
        nor=norvec[i]
        grd = grd_outer[i].subgridspec(2,2, wspace=0.1, hspace=0.0,height_ratios=[1,5])
        axx={}
        axx[0,0]=fig.add_subplot(grd[0,0])
        axx[0,1]=fig.add_subplot(grd[0,1])
        axx[1,0]=fig.add_subplot(grd[1,0],sharex=axx[0,0])
        axx[1,1]=fig.add_subplot(grd[1,1],sharex=axx[0,1])
        if i==0:
                axx[0,0].set_title('Observations\n',fontsize=fsz)
                axx[0,1].set_title('NorESM model\n',fontsize=fsz)
        if 'Barents' in tit:
            nor_x=nor.LATITUDE
            obs_x=obs.LATITUDE
        else:
            nor_x=nor.LONGITUDE
            obs_x=obs.LONGITUDE
        velcont,salcont=cont_partial(obs,obs_x,axx[1,0],yspl[i],ymax[i],sigmax_obs)
        cont_partial(obs,obs_x,axx[0,0],0,yspl[i],sigmax_obs)
        cont_partial(nor,nor_x,axx[1,1],yspl[i],ymax[i],sigmax_obs)
        cont_partial(nor,nor_x,axx[0,1],0,yspl[i],sigmax_obs)
        if 'OSNAP' in tit:
            axx[0,0].plot([oslim_obs,oslim_obs],[0,50],color='k',linewidth=3)
            axx[0,1].plot([oslim,oslim],[0,30],color='k',linewidth=3)
            axx[1,0].text(-19,2750,tit,color='white',fontsize=fsz)
            axx[0,0].text(-15,80,'AWS',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[0,1].text(-15,80,'AWS',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[0,0].text(-45,-10,'PWS',color='k',fontsize=fsz-2)
            axx[0,1].text(-43,-10,'PWS',color='k',fontsize=fsz-2)
            axx[1,0].text(-40,1500,'DWS',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[1,1].text(-40,1500,'DWS',color='k',backgroundcolor='w',fontsize=fsz-2)
        elif 'Fram' in tit:
            for iii in range(2):
                axx[iii,0].axvline(fslim_obs,color='k',linewidth=3)
                axx[iii,1].axvline(fslim,color='k',linewidth=3)
            axx[1,0].text(-18,2600,tit,color='white',fontsize=fsz)
            axx[1,0].text(8,1500,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[1,1].text(9,1500,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[1,0].text(-3,1500,'PWN',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[1,1].text(-2,1500,'PWN',color='k',backgroundcolor='w',fontsize=fsz-2)
        else:
            axx[1,0].text(74.5,450,'Barents Sea\nOpening',color='white',fontsize=fsz)
            axx[1,0].text(72,150,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
            axx[1,1].text(72,150,'AWN',color='k',backgroundcolor='w',fontsize=fsz-2)
        axx[0,1].set_yticklabels('')
        axx[1,1].set_yticklabels('')
        axx[1,0].set_ylabel('depth [m]',fontsize=14)
    fig.text(0.5, 0.64, 'Longitude [$^\circ$W]', ha='center',fontsize=fsz)
    fig.text(0.5, 0.352, 'Longitude [$^\circ$W]', ha='center',fontsize=fsz)
    fig.text(0.5, 0.065, 'Latitude [$^\circ$N]', ha='center',fontsize=fsz)
    xc=0.95
    cwi=0.025
    cle=0.325
    caxit_vel=fig.add_axes([xc,0.5,cwi,cle])
    caxit_sal=fig.add_axes([xc,0.15,cwi,cle])
    colorbar(velcont,label='velocity [m/s]',cax=caxit_vel)
    colorbar(salcont,label='Salinity',cax=caxit_sal)
    savefig(figdir_paper+'CompSect_Obs_NorESM_all.png',bbox_inches='tight')
    savefig(figdir_paper+'CompSect_Obs_NorESM_all.pdf',bbox_inches='tight')


comp_obs_model_all()
############################################################################################
################  NORESM  WATER MASS PARTITIONING   #############################################

AWS={}
PWS={}
DWS={}
AWN={}
PWN={}

pws_depth=100

xray=osnap
for var in ['TRANS','PSAL','PTMP','PDEN']:
    if var=='TRANS':
        DWS[var]=xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS[var]=xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS[var]=xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')

    else:
        DWS[var]=(xray[var]*xray['TRANS']).where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')

#split up FS into PW and AW using definition in Tsubouchi et al. 2018, which is boundary between "EGC" + "Middle", 2W
# Here I'm using 4E, as it better maximizes transport...? Check on this...
#call all BSO water AW.
for var in ['PSAL','PTMP','PDEN','TRANS']:
    if var=='TRANS':
        PWN[var]=-fs['TRANS'].where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN[var]=-fs['TRANS'].where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
    else:
        PWN[var]=(fs[var]*fs['TRANS']).where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')/fs['TRANS'].where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN[var]=((fs[var]*fs['TRANS']).where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')+(bso[var]*bso['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE'))/(fs['TRANS'].where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')+bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE'))




def plot_NorESM_WMprops_tvar():
    f,axx=subplots(4,1,figsize=(12,20),sharex=True)
    for ii,var in enumerate(['TRANS','PSAL','PTMP','PDEN']):
            PWS[var].plot(label='PWS',ax=axx[ii])
            AWS[var].plot(label='AWS',ax=axx[ii])
            DWS[var].plot(label='DWS',ax=axx[ii])
            AWN[var].plot(label='AWN',ax=axx[ii])
            PWN[var].plot(label='PWN',ax=axx[ii])
            axx[ii].set_ylabel(var)
            axx[ii].set_xlabel('')
    axx[0].legend()
    savefig(figdir+'NorESM_WMtvar_all.png',bbox_inches='tight')


plot_NorESM_WMprops_tvar()

############################################################################################
################  OBS WATER MASS PARTITIONING   #############################################

AWS_obs={}
PWS_obs={}
DWS_obs={}
AWN_obs={}
AWN_obs_fs={}
AWN_obs_bso={}
PWN_obs={}

xray=osnap_obs
sigmax_obs=osnap_obs.SIGMEAN.values
for var in ['PSAL','PTMP','PDEN','TRANS']:
    if var=='TRANS':
        DWS_obs[var]=xray['TRANS'].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS_obs[var]=xray['TRANS'].where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS_obs[var]=xray['TRANS'].where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
    else:
        DWS_obs[var]=(xray[var]*xray['TRANS']).where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS_obs[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS_obs[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')


#split up FS into PW and AW using definition in Tsubouchi et al. 2018, which is boundary between "EGC" + "Middle", 2W
for var in ['TRANS','PSAL','PTMP','PDEN',]:
    if var=='TRANS':
        PWN_obs[var]=-fs_obs['TRANS'].where(fs_obs.LONGITUDE<fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs[var]=-fs_obs['TRANS'].where(fs_obs.LONGITUDE>=fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')-bso_obs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs_fs[var]=-fs_obs['TRANS'].where(fs_obs.LONGITUDE>=fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
    else:
        PWN_obs[var]=(fs_obs[var]*fs_obs['TRANS']).where(~isnan(fs_obs[var])).where(fs_obs.LONGITUDE<fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/fs_obs['TRANS'].where(~isnan(fs_obs[var])).where(fs_obs.LONGITUDE<fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs_fs[var]=(fs_obs[var]*fs_obs['TRANS']).where(~isnan(fs_obs[var])).where(fs_obs.LONGITUDE>=fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/fs_obs['TRANS'].where(~isnan(fs_obs[var])).where(fs_obs.LONGITUDE>=fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs_bso[var]=(bso_obs[var]*bso_obs['TRANS']).where(~isnan(bso_obs[var])).sum(dim='DEPTH').sum(dim='LONGITUDE')/bso_obs['TRANS'].where(~isnan(bso_obs[var])).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs[var]=AWN_obs_fs[var]*fs_obs['TRANS'].where(fs_obs.LONGITUDE>=fslim_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/(-AWN_obs['TRANS'])+AWN_obs_bso[var]*bso_obs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')/(-AWN_obs['TRANS'])


def plot_NorESM_WMprops_clim():
    f,axx=subplots(4,1,figsize=(6,20),sharex=True)
    for ii,var in enumerate(['TRANS','PSAL','PTMP','PDEN']):
            PWS[var].groupby('TIME.month').mean('TIME').plot(label='PWS NorESM',ax=axx[ii],linewidth=3)
            AWS[var].groupby('TIME.month').mean('TIME').plot(label='AWS NorESM',ax=axx[ii],linewidth=3)
            DWS[var].groupby('TIME.month').mean('TIME').plot(label='DWS NorESM',ax=axx[ii],linewidth=3)
            AWN[var].groupby('TIME.month').mean('TIME').plot(label='AWN NorESM',ax=axx[ii],linewidth=3)
            PWN[var].groupby('TIME.month').mean('TIME').plot(label='PWN NorESM',ax=axx[ii],linewidth=3)

            PWS_obs[var].groupby('TIME.month').mean('TIME').plot(label='PWS Obs',ax=axx[ii],linewidth=3,linestyle='--',color='C0')
            AWS_obs[var].groupby('TIME.month').mean('TIME').plot(label='AWS Obs',ax=axx[ii],linewidth=3,linestyle='--',color='C1')
            DWS_obs[var].groupby('TIME.month').mean('TIME').plot(label='DWS Obs',ax=axx[ii],linewidth=3,linestyle='--',color='C2')
            AWN_obs[var].groupby('TIME.month').mean('TIME').plot(label='AWN Obs',ax=axx[ii],linewidth=3,linestyle='--',color='C3')
            PWN_obs[var].groupby('TIME.month').mean('TIME').plot(label='PWN Obs',ax=axx[ii],linewidth=3,linestyle='--',color='C4')
            axx[ii].set_ylabel(var)
            axx[ii].set_xlabel('month')
    axx[0].legend(ncol=2,loc=(1.05,0.1))
    savefig(figdir+'NorESM_WM_obscomp_clim.png',bbox_inches='tight')

plot_NorESM_WMprops_clim()

plot_NorESM_WMprops_clim()


############################################################################################
################  PLOT AND COMPARE WATER MASS PARTITIONING   #############################################


salvec=linspace(31,36,103)
tmpvec=linspace(-3,16,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)
SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
SA_vec_1000=gsw.SA_from_SP(salvec,1e3*ones(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
sigma1mat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3
        sigma1mat[jj,ii]=gsw.sigma1(SA_vec[ii],CT_vec[jj])



def plot_TS_bylayer():
    figure(figsize=(5,4))
    scatter(PWS['PSAL'].mean(dim='TIME'),PWS['PTMP'].mean(dim='TIME'),s=PWS['TRANS'].mean(dim='TIME')**2*4,zorder=50,linewidth=3,label='PWS')
    scatter(AWS['PSAL'].mean(dim='TIME'),AWS['PTMP'].mean(dim='TIME'),s=AWS['TRANS'].mean(dim='TIME')**2*4,zorder=50,linewidth=3,label='AWS')
    scatter(DWS['PSAL'].mean(dim='TIME'),DWS['PTMP'].mean(dim='TIME'),s=DWS['TRANS'].mean(dim='TIME')**2*4,zorder=50,linewidth=3,label='DWS')
    scatter(AWN['PSAL'].mean(dim='TIME'),AWN['PTMP'].mean(dim='TIME'),s=AWN['TRANS'].mean(dim='TIME')**2*4,zorder=50,linewidth=3,label='AWN')
    scatter(PWN['PSAL'].mean(dim='TIME'),PWN['PTMP'].mean(dim='TIME'),s=PWN['TRANS'].mean(dim='TIME')**2*4,zorder=51,linewidth=3,label='PWN')
    contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=5)
    xlabel('salinity')
    ylabel('pot.temperature [$^\circ$C]')
    xlim(34,36)
    ylim(-2,12)
    lgnd=legend(loc=(1.05,0.2),ncol=2)
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [40]
    title('Transport-weighted water mass properties in NorESM')
    savefig(figdir+'NorESM_WMTS.png',bbox_inches='tight')

plot_TS_bylayer()



# ################################################################################################################################
# ################################################################################################################################
# ###########################################      SAVE WM props    #######################################################
# ################################################################################################################################
# ################################################################################################################################

WM_obs={}
for ii,xrw in enumerate([PWS_obs,AWS_obs,DWS_obs,PWN_obs,AWN_obs]):
    WM_obs[ii]=xr.concat([xrw['PDEN'],xrw['PSAL'],xrw['PTMP'],xrw['TRANS']],pd.Index(['PDEN','PSAL','PTMP','TRANS'],name='var'))

WMall_obs=xr.concat([WM_obs[ww] for ww in WM_obs],pd.Index(['PWS','AWS','DWS','PWN','AWN'],name='WM'))

WMall_obs=WMall_obs.to_dataset('var')
# WMall_obs.to_netcdf(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc','w')

WM={}
for ii,xrw in enumerate([PWS,AWS,DWS,AWN,PWN]):
    WM[ii]=xr.concat([xrw['PDEN'],xrw['PSAL'],xrw['PTMP'],xrw['TRANS']],pd.Index(['PDEN','PSAL','PTMP','TRANS'],name='var'))

WMall=xr.concat([WM[ww] for ww in WM],pd.Index(['PWS','AWS','DWS','AWN','PWN'],name='WM'))
WMall=WMall.to_dataset('var')

# WMall.to_netcdf(datadir+'NorESM/NorESM_WMs_1912.nc','w')

WMall['TRANS'].sel(WM='PWS').mean(dim='TIME')+WMall['TRANS'].sel(WM='AWS').mean(dim='TIME')+WMall['TRANS'].sel(WM='DWS').mean(dim='TIME')

(WMall['TRANS'].sel(WM='PWN').mean(dim='TIME')+WMall['TRANS'].sel(WM='AWN').mean(dim='TIME'))

################################################################################################################################
################################################################################################################################
###########################################      PLOT TS SPACE    #######################################################
################################################################################################################################
################################################################################################################################



coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}

WM_nor=WMall
WM_obs=WMall_obs

# WM_lit=xr.Dataset({'PSAL':('WM',[35.2,35,33.8,33.8,35]),'PTMP':('WM',[0,0,0,0,4]),
#                    'TRANS':('WM',[17.2,-12.7,-3.8,8,-8]),'COLOR':('WM',['red','grey','royalblue','purple','orange'])},
#                    coords={'WM':['AWS', 'DWS', 'PWS', 'PWN', 'AWN'],})



def plot_TS_bylayer():
    f,[ax1,ax2]=subplots(1,2,figsize=(11,5),sharex=True,sharey=True)
    ii=0
    for tt in WM_nor.TIME:
        for wm in WM_nor.WM:
            ax2.scatter(WM_nor['PSAL'].sel(WM=wm).sel(TIME=tt).values,WM_nor['PTMP'].sel(WM=wm).sel(TIME=tt).values,
            s=WM_nor['TRANS'].sel(WM=wm).sel(TIME=tt).values**2*4,zorder=50,linewidth=3,label=str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
        if ii==0:
            lgnd=ax2.legend(loc=(1.05,0.2))
        ii+=1
    for wm in WM_nor.WM:
        ax2.scatter(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,
                s=WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label='NorESM : '+str(wm.values),color='k',zorder=100)
        if ('AWS' in wm) | ('PWN' in wm):
            ax2.plot(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
        else:
            ax2.plot(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
        # else:
        #     ax2.plot(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=10,color='k')
    for wm in WM_obs.WM:
        for tt in WM_obs.TIME:
            ax1.scatter(WM_obs['PSAL'].sel(WM=wm).sel(TIME=tt).values,WM_obs['PTMP'].sel(WM=wm).sel(TIME=tt).values,
            s=WM_obs['TRANS'].sel(WM=wm).sel(TIME=tt).values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
        ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,
                s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label=str(wm.values),color='k',zorder=100)
        if ('AWS' in wm) | ('PWN' in wm):
            ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
        # else:
        #     ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=10,color='k',zorder=51)

    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [40]
    for axx in [ax1,ax2]:
        axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
        axx.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    xlim(34,36)
    ylim(-2,12)
    f.suptitle('Transport-weighted water mass properties\n',fontsize=16)
    ax1.set_title('Observations',fontsize=14)
    ax2.set_title('NorESM',fontsize=14)
    savefig(figdir+'TS_comp.png',bbox_inches='tight')
    savefig(figdir+'TS_comp.pdf',bbox_inches='tight')

plot_TS_bylayer()
M96={}
M96['WMS']=['AW_IS','AW_BS','AW_FS','AE_PW','AE_AAW','AW_IW','NE_RAW','ED_PW','ED_AAW','ED_RAW','AW_IR','AW_EI','DW']
M96['PTMP']=[7,4,2,-1,0.5,0,1.5,-1,0.5,1,4,2.5,-1]
M96['PSAL']=[35.3,35.05,35,34.3,34.85,34.89,34.95,34.3,34.85,34.92,34,34.9,34.9]
M96['TRANS']=[6.8,-1.6,-3.4,1.5,1.5,0.17,1.1,-1.6,-2,-0.8,0.9,-0.7,-2]

sign=array(M96['TRANS'])/array([abs(tt) for tt in M96['TRANS']])

sign
sum(M96['TRANS'])

def plot_M96_TS():
    f,axx=subplots(1,1,figsize=(6,5),sharex=True,sharey=True)
    for ii in range(len(M96['WMS'])):
        axx.scatter(M96['PSAL'][ii],M96['PTMP'][ii],s=M96['TRANS'][ii]**2*50,zorder=50,linewidth=3,label=M96['WMS'][ii])
        if sign[ii]>0:
            axx.plot(M96['PSAL'][ii],M96['PTMP'][ii],'k+',zorder=51)
        else:
            axx.plot(M96['PSAL'][ii],M96['PTMP'][ii],'k_',zorder=51)
    axx.legend(loc=(1.05,0))
    axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    xlim(34,36)
    ylim(-2,12)
    f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)

    f.suptitle('Mauritzen 1996 (Pt.II)\nTransport-weighted water mass properties\n',fontsize=16)
    savefig(figdir+'TS_Mauritzen.png',bbox_inches='tight')
    savefig(figdir+'TS_Mauritzen.pdf',bbox_inches='tight')

plot_M96_TS()

############################################################################################
################  PLOT SEASONAL EVOLUTION OF NORESM WMS   #############################################
############################################################################################
# This was not v. helpful, 12 panels all looked v similar, hard to compare
# def plot_TS_NORESM_CLIM():
#     axx={}
#     f,[(axx[0],axx[1],axx[2]),(axx[3],axx[4],axx[5]),(axx[6],axx[7],axx[8]),(axx[9],axx[10],axx[11])]=subplots(4,3,figsize=(15,15),sharex=True,sharey=True)
#     for wm in WM_nor.WM:
#         WM_clim=WM_nor.groupby('TIME.month').mean('TIME').sel(WM=wm)
#         for ii in range(12):
#             axx[ii].scatter(WM_clim['PSAL'][ii].values,WM_clim['PTMP'][ii].values,
#                 s=WM_clim['TRANS'][ii].values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
#             if ('AWS' in wm) | ('PWN' in wm):
#                 axx[ii].plot(WM_clim['PSAL'][ii].values,WM_clim['PTMP'][ii].values,'+',markersize=12,color='k')
#             axx[ii].contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
#             axx[ii].contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=5)
#     lgnd=axx[5].legend(loc=(2.3,0.2))
#     for ii in range(5):
#         lgnd.legendHandles[ii]._sizes = [40]
#     f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
#     f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
#     xlim(34,36)
#     ylim(-3,12)
#     # f.suptitle('Transport-weighted water mass properties\n',fontsize=16)
#     # savefig(figdir+'TS_comp.png',bbox_inches='tight')
#     # savefig(figdir+'TS_comp.pdf',bbox_inches='tight')
#
#
# plot_TS_NORESM_CLIM()

def plot_TS_NORESM_ALL():
    f,axx=subplots(1,1,figsize=(10,8))
    for wm in WM_nor.WM:
        for tt in WM_nor.TIME:
            axx.scatter(WM_nor['PSAL'].sel(WM=wm).sel(TIME=tt).values,WM_nor['PTMP'].sel(WM=wm).sel(TIME=tt).values,
            s=WM_nor['TRANS'].sel(WM=wm).sel(TIME=tt).values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
        # if ('AWS' in wm) | ('PWN' in wm):
        #     axx.plot(WM_clim['PSAL'][3*ii:3*ii+3].mean(dim='month').values,WM_clim['PTMP'][3*ii:3*ii+3].mean(dim='month').values,'+',markersize=12,color='k',zorder=100)
        if tt==0:
            lgnd=axx.legend(loc=(1.05,0.2))
            for jj in range(5):
                lgnd.legendHandles[jj]._sizes = [40]
    axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    axx.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=5)

    f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    xlim(34,36)
    ylim(-3,12)
    f.suptitle('Seasaonal evolution of transport-weighted water mass properties in NorESM',fontsize=16)
    savefig(figdir+'TS_all_NorESM.png',bbox_inches='tight')
    savefig(figdir+'TS_all_NorESM.pdf',bbox_inches='tight')


plot_TS_NORESM_ALL()


def plot_TS_NORESM_CLIM():
    f,axx=subplots(1,1,figsize=(10,8))
    for ii in range(4):
        for wm in WM_nor.WM:
            WM_clim=WM_nor.groupby('TIME.month').mean('TIME').sel(WM=wm)
            axx.scatter(WM_clim['PSAL'][3*ii:3*ii+3].mean(dim='month').values,WM_clim['PTMP'][3*ii:3*ii+3].mean(dim='month').values,
                s=WM_clim['TRANS'][3*ii:3*ii+3].mean(dim='month').values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
            if ('AWS' in wm) | ('PWN' in wm):
                axx.plot(WM_clim['PSAL'][3*ii:3*ii+3].mean(dim='month').values,WM_clim['PTMP'][3*ii:3*ii+3].mean(dim='month').values,'+',markersize=12,color='k',zorder=100)
        if ii==0:
            lgnd=axx.legend(loc=(1.05,0.2))
            for jj in range(5):
                lgnd.legendHandles[jj]._sizes = [40]
    axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    axx.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=5)

    f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    xlim(34,36)
    ylim(-3,12)
    f.suptitle('Seasaonal evolution of transport-weighted water mass properties in NorESM',fontsize=16)
    savefig(figdir+'TS_seas_NorESM.png',bbox_inches='tight')
    savefig(figdir+'TS_seas_NorESM.pdf',bbox_inches='tight')


plot_TS_NORESM_CLIM()

def plot_XPORT_S_bylayer():
    f,[ax1,ax2]=subplots(1,2,figsize=(11,5),sharex=True,sharey=True)
    for wm in WM_nor.WM:
        ax2.scatter(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values,
                s=WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
        ax2.plot([WM_nor['PSAL'].sel(WM=wm).min(dim='TIME'),WM_nor['PSAL'].sel(WM=wm).max(dim='TIME')],[WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values,WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values],color=coldic[str(wm.values)],)
        ax2.plot([WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME'),WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME')],[WM_nor['TRANS'].sel(WM=wm).min(dim='TIME').values,WM_nor['TRANS'].sel(WM=wm).max(dim='TIME').values],color=coldic[str(wm.values)],)
    for wm in WM_obs.WM:
        ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values,
                s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,zorder=50,linewidth=3,label=str(wm.values),color=coldic[str(wm.values)])
        ax1.plot([WM_obs['PSAL'].sel(WM=wm).min(dim='TIME'),WM_obs['PSAL'].sel(WM=wm).max(dim='TIME')],[WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values,WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values],color=coldic[str(wm.values)],)
        ax1.plot([WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME'),WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME')],[WM_obs['TRANS'].sel(WM=wm).min(dim='TIME').values,WM_obs['TRANS'].sel(WM=wm).max(dim='TIME').values],color=coldic[str(wm.values)],)
    lgnd=ax1.legend(loc=(2.3,0.3))
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [40]
    f.text(0.05, 0.5, 'volume transport [Sv]', va='center', rotation='vertical',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    xlim(33.5,35.75)
    for axx in [ax1,ax2]:
        axx.axhline(0,color='k',linewidth=3)
    # ylim(-2,12)
    f.suptitle('Volume transport and salinity by water mass\n',fontsize=16)
    ax1.set_title('Observations',fontsize=14)
    ax2.set_title('NorESM',fontsize=14)
    savefig(figdir+'XPORT_S_comp.png',bbox_inches='tight')

plot_XPORT_S_bylayer()



# THIS WAS FOR DEMOING TO FEILI -- WOULD RE-EVALUATE NOW THAT I HAVE MORE OBS/ MAYBE ITS NOT THE POINT
# CODE COULD BE USEFUL.

Se={'AWS': 0.06, 'PWS': 0.8, 'DWS': 0.1, 'AWN': 0.2, 'PWN': 0.8, 'FW': 0, 'SI': 0}
Ue={'AWS': 2.0,'PWS': 2.0,'DWS': 2.0,'AWN': 2.0,'PWN': 2.0,'FW': 0.01,'SI': 0.01}



WM_obs['COLOR']=['royalblue','red','grey','purple','orange']

WM_obs
def barplot():
    f,[ax1,ax2,ax3]=subplots(3,1,figsize=(16,11))
    ax1.bar(WM_obs['WM'],WM_obs['TRANS'].mean(dim='TIME'),color=WM_obs['COLOR'].values)#,yerr=[Ue[kk] for kk in WM_lit.WM.values],capsize=10)
    ax1.axhline(0,color='k',linewidth=2)
    ax1.plot(WM_nor['WM'],WM_nor['TRANS'].mean(dim='TIME'),'d',color='k',markersize=20)
    ax2.bar(WM_obs['WM'],WM_obs['PSAL'].mean(dim='TIME'),color=WM_obs['COLOR'].values)#,yerr=[Se[kk] for kk in WM_lit.WM.values],capsize=10)
    ax2.set_ylim(32.5,36)
    ax2.plot(WM_nor['WM'],WM_nor['PSAL'].mean(dim='TIME'),'d',color='k',markersize=20)
    ax3.bar(WM_obs['WM'],(WM_obs['TRANS']*WM_obs['PSAL']).mean(dim='TIME'),color=WM_obs['COLOR'].values)#,yerr=[Ue[kk]*WM_lit['PSAL'].sel(WM=kk)+Se[kk]*WM_lit['TRANS'].sel(WM=kk) for kk in WM_lit.WM.values],capsize=10)
    ax3.plot(WM_nor['WM'],(WM_nor['PSAL']*WM_nor['TRANS']).mean(dim='TIME'),'d',color='k',markersize=20)
    ax3.axhline(0,color='k',linewidth=2)
    ax1.set_ylabel('Transport, U [Sv]')
    ax2.set_ylabel('Salinity, S')
    ax3.set_ylabel('Salinity budget term, U*S')
    savefig(figdir+'Barchart_synth.png',bbox_inches='tight')


barplot()

# THIS IS FOR NORESM ONLY PLOTTING -- KEEPING AROUND FOR POSTERITY
# ################################################################################################################################
# ################################################################################################################################
# ###########################################      PLOT MEAN SECTION    #######################################################
# ################################################################################################################################
# ################################################################################################################################
# sigmax=osnap.SIGMEAN.values
#
# sigmax
# oslim=-40
# fslim=6
#
#
# univec['tmp']=['pot. temperature', arange(-1,11,0.5),cm.RdYlBu_r,range(0,11,2),'[$^\\circ$C]']
# univec['sal']=['salinity',arange(34,35.6,0.1),cm.PiYG_r, array([34., 34.5 ,35,35.5,]),'']
# vmax=0.1
# univec['VELO']=['across track velocity',arange(-vmax,vmax+0.01,0.01),cm.RdBu_r,arange(-vmax,vmax+0.05,0.05),'[m/s]']
#
# # VELOCITY, SALINITY, TEMPERATURE AND DENSITY MEANS OVER PERIOD OF OVERLAP
# def plot_each(dat,xxvar,axx,var,var2,ymax):
#     filled=axx.pcolor(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),cmap=univec[var2][2])#,extend='both')#univec[var2][1],
#     axx.contour(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),levels=univec[var2][1][::2],colors='k')
#     axx.contour(xxvar,dat.DEPTH,dat['PDEN'].mean(dim='TIME'),levels=[sigmax],colors='k',linewidths=4) # add NorESM isopycnal of maximum overturning (for 2014-2016)
#     axx.set_facecolor('k')
#     axx.set_ylim(ymax,0)
#     return filled
#
# def plot_all_secs(var,var2,tit):
#     f,[ax1,ax2,ax3,ax4]=subplots(1,4,figsize=(20,4))
#     plot_each(fs,fs.LONGITUDE,ax1,var,var2,3000)
#     ax1.set_xlabel('Longitude [$^\circ$W]')
#     ax1.set_title('Fram Strait',fontsize=14)
#
#     plot_each(bso,bso.LATITUDE,ax2,var,var2,600)
#     ax2.set_xlabel('Latitude [$^\circ$N]')
#     ax2.set_title('Barents Sea Opening',fontsize=14)
#     fill=plot_each(osnap,osnap.LONGITUDE,ax3,var,var2,3500)
#     ax3.set_xlabel('Longitude [$^\circ$W]')
#     ax3.set_title('OSNAP East',fontsize=14)
#     fill=plot_each(ns,ns.LONGITUDE,ax4,var,var2,300)
#     ax4.set_xlabel('Longitude [$^\circ$W]')
#     ax4.set_title('North Sea',fontsize=14)
#     caxit=f.add_axes([0.93,0.1,0.01,0.8])
#     ax1.axvline(fslim,color='k',linewidth=4)
#     ax3.plot([oslim, oslim],[0,100],color='k',linewidth=4)
#     colorbar(fill,label=tit,cax=caxit)
#     savefig(figdir+'NorESM_sections_'+var2+'.png',bbox_inches='tight')
#
#
#
# plot_all_secs('PTMP','tmp','pot. temperature [$^\circ$C]')
#
# plot_all_secs('PSAL','sal','salinity')
#
# plot_all_secs('PDEN','pden','potential density [kg m$^{-3}$]')
#
# plot_all_secs('VELO','VELO','velocity [m s$^{-1}$]')
