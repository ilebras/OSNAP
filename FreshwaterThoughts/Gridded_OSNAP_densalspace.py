from aux_funcs import *

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))
eastind=dat.LONGITUDE>-45
westind=dat.LONGITUDE<-45
allind=dat.LONGITUDE<0
psi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Streamfunction_201408_201604_2018.nc')


salvec=linspace(31,36,103)
tmpvec=linspace(-3,16,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3


#
sref_west=34.85
sref_east=34.97
sref_all=34.92

dat['TRANS']=dat.VELO*dat.AREA/1e6
dat['FWT']=dat['TRANS']*(sref_all-dat['PSAL'])/sref_all*1e3

#
# sal_lev=arange(31,35.6,0.01)
# sal_xport=zeros((len(dat.TIME),len(dat.LONGITUDE),len(sal_lev)))
# sal_fwt=zeros((len(dat.TIME),len(dat.LONGITUDE),len(sal_lev)))
# for ii in range(len(sal_lev)-1):
#     sal_xport[:,:,ii]=(dat['TRANS']).where(dat.PSAL>=sal_lev[ii]).where(dat.PSAL<sal_lev[ii+1]).sum(dim='DEPTH')
#     sal_fwt[:,:,ii]=(dat['FWT']).where(dat.PSAL>=sal_lev[ii]).where(dat.PSAL<sal_lev[ii+1]).sum(dim='DEPTH')
#
# den_lev=psi.LEVEL.values
# den_fwt=NaN*ones((len(dat.TIME),len(dat.LONGITUDE),len(den_lev)))
# den_xport=NaN*ones((len(dat.TIME),len(dat.LONGITUDE),len(den_lev)))
# den_sal=NaN*ones((len(dat.TIME),len(dat.LONGITUDE),len(den_lev)))
# den_tmp=NaN*ones((len(dat.TIME),len(dat.LONGITUDE),len(den_lev)))
# den_den=NaN*ones((len(dat.TIME),len(dat.LONGITUDE),len(den_lev)))
# den_area=NaN*ones((len(dat.TIME),len(dat.LONGITUDE),len(den_lev)))
# for ii in range(len(den_lev)-1):
#     den_xport[:,:,ii]=(dat['TRANS']).where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1]).sum(dim='DEPTH')
#     den_fwt[:,:,ii]=dat['FWT'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1]).sum(dim='DEPTH')
#     den_sal[:,:,ii]=(dat['PSAL'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1])*dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1])).sum(dim='DEPTH')/dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1]).sum(dim='DEPTH')
#     den_tmp[:,:,ii]=(dat['PTMP'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1])*dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1])).sum(dim='DEPTH')/dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1]).sum(dim='DEPTH')
#     den_den[:,:,ii]=(dat['PDEN'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1])*dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1])).sum(dim='DEPTH')/dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1]).sum(dim='DEPTH')
#     den_area[:,:,ii]=dat['AREA'].where(dat.PDEN>=den_lev[ii]).where(dat.PDEN<den_lev[ii+1]).sum(dim='DEPTH').T
#
#
# den_space=xr.Dataset({'xport': (['date','LONGITUDE','den_lev'],  den_xport),'fwt': (['date','LONGITUDE','den_lev'],  den_fwt),'pden': (['date','LONGITUDE','den_lev'],  den_den),'sal': (['date','LONGITUDE','den_lev'],  den_sal),'tmp': (['date','LONGITUDE','den_lev'],  den_tmp),'area': (['date','LONGITUDE','den_lev'],  den_area),},coords={'date':dat.TIME.values,'LONGITUDE': dat.LONGITUDE.values,'den_lev': den_lev})
# sal_space=xr.Dataset({'xport': (['date','LONGITUDE','sal_lev'],  sal_xport),'fwt': (['date','LONGITUDE','sal_lev'],  sal_fwt)},coords={'date':dat.TIME.values,'LONGITUDE': dat.LONGITUDE.values,'sal_lev': sal_lev})

sal_cols = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(sal_cols,position=[0,0.6,0.69,0.9,1],bit=True)
slevs=array([34,34.2,34.4,34.6,34.8,34.85,34.9,34.95,35,35.05,35.1,35.15,35.2,35.25,35.3,35.35,35.4])

def plotsec_east():
    f,axx=subplots(2,1,sharex=True,sharey=True,figsize=(12,12))

    cxx=0.92
    cy=0.325
    cax1=f.add_axes([cxx,0.55,0.0175,cy])
    vel=axx[0].contourf(dat.LONGITUDE,dat.DEPTH,dat.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
    axx[0].contour(dat.LONGITUDE,dat.DEPTH,dat.VELO.mean(dim='TIME'),levels=[0],colors='k')

    # clabel(c2,fmt='%1.2f')
    axx[0].set_ylabel('depth [m]')
    colorbar(vel,ax=axx[0],cax=cax1,label='Velocity [m/s]')
    axx[0].fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
    # axx[0].text(-44,3200,'Greenland',color='white',fontsize=16)
    # axx[0].text(-34.5,3200,'Mid-Atlantic Ridge',color='white',fontsize=16)
    # axx[0].text(-11,3200,'Scotland',color='white',fontsize=16)
    # axx[1].text(-19,2800,'Mean 2014-2016',color='white',fontsize=15)

    cax2=f.add_axes([cxx,0.125,0.0175,cy])
    sal=axx[1].contourf(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),slevs,cmap=sal_cmap,extend='both')
    c1=axx[1].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[34.9],colors='grey')

    axx[1].set_ylabel('depth [m]')
    axx[1].fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
    clabel(c1,fmt='%1.1f')
    # c2=axx[1].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[sref_east],colors='k')


    axx[1].set_xlabel('Longitude')
    colorbar(sal,ax=axx[1],ticks=slevs[::2],cax=cax2,label='Salinity')
    axx[0].set_ylim(3500,0)
    savefig(figdir+'Freshwater/Salvecsec_all.png',bbox_inches='tight',dpi=300)
    axx[0].set_ylim(500,0)
    axx[0].set_xlim(-56,-50)
    savefig(figdir+'Freshwater/Salvecsec_all_zoom.png',bbox_inches='tight',dpi=300)
    axx[0].set_ylim(3500,0)
    axx[0].set_xlim(-44,-7)
    c2=axx[0].contour(dat.LONGITUDE,dat.DEPTH,dat.PDEN.mean(dim='TIME'),levels=[27.53],colors='k',linewidths=3)
    c2=axx[1].contour(dat.LONGITUDE,dat.DEPTH,dat.PDEN.mean(dim='TIME'),levels=[27.53],colors='k',linewidths=3)
    clabel(c2,fmt='%1.2f')
    [axx[ii].axvline(-40,linewidth=4,color='grey') for ii in [0,1]]
    axx[0].set_ylim(500,0)
    axx[0].set_xlim(-44,-30)
    savefig(figdir+'Freshwater/Salvecsec_labellayer_east_zoom.png',bbox_inches='tight',dpi=300)
    axx[0].set_ylim(3500,0)
    axx[0].set_xlim(-44,-7)
    axx[1].text(-23,400,'upper layer',fontsize=17,color='white')
    axx[1].text(-37,1500,'lower layer',fontsize=17)
    axx[0].set_title('OSNAP East')
    savefig(figdir+'Freshwater/Salvecsec_labellayer_east.png',bbox_inches='tight',dpi=300)


plotsec_east()


lonbnd=-44
denlim=27.53

layer={}
layer['upper']={}
layer['upper']['xport']=dat.TRANS.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim).sum(dim='DEPTH')
layer['upper']['sal']=(dat.PSAL.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim)*dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim)).sum(dim='DEPTH')/dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim).sum(dim='DEPTH')
layer['upper']['tmp']=(dat.PTMP.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim)*dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim)).sum(dim='DEPTH')/dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim).sum(dim='DEPTH')
layer['upper']['pden']=(dat.PDEN.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim)*dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim)).sum(dim='DEPTH')/dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')<denlim).sum(dim='DEPTH')
layer['lower']={}
layer['lower']['xport']=dat.TRANS.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim).sum(dim='DEPTH')
layer['lower']['sal']=(dat.PSAL.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim)*dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim)).sum(dim='DEPTH')/dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim).sum(dim='DEPTH')
layer['lower']['tmp']=(dat.PTMP.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim)*dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim)).sum(dim='DEPTH')/dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim).sum(dim='DEPTH')
layer['lower']['pden']=(dat.PDEN.mean(dim='TIME').where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim)*dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim)).sum(dim='DEPTH')/dat.AREA.where(dat.LONGITUDE>lonbnd).where(dat.PDEN.mean(dim='TIME')>=denlim).sum(dim='DEPTH')

NA_wbnd=-40

def twolayer_timemean():
    f,axx=subplots(4,1,sharex=True,figsize=(5,12))
    axx[0].plot(dat.LONGITUDE,layer['lower']['xport'].cumsum(dim='LONGITUDE'),label='Lower layer')
    axx[0].plot(dat.LONGITUDE[::-1],layer['upper']['xport'][::-1].cumsum(dim='LONGITUDE'),label='Upper layer')
    axx[0].set_ylabel('Cumulative transport [Sv]')
    axx[0].legend(loc=(1.05,0.2),fontsize=16)
    axx[0].set_title('Layer averaged properties',fontsize=16)

    axx[0].axhline(0,color='k')
    axx[1].plot(dat.LONGITUDE,layer['lower']['sal'])
    axx[1].plot(dat.LONGITUDE,layer['upper']['sal'])
    fs=20
    axx[1].text(-30,35.3,'AW',color='C1',fontsize=fs)
    axx[1].text(-45,35.3,'PW',color='C1',fontsize=fs)
    axx[1].set_ylabel('Salinity')
    axx[1].set_ylim(34.5,35.5)


    axx[2].plot(dat.LONGITUDE,layer['lower']['tmp'])
    axx[2].plot(dat.LONGITUDE,layer['upper']['tmp'])
    axx[2].set_ylabel('Theta [$^\circ$C]')

    axx[3].plot(dat.LONGITUDE,layer['lower']['pden'])
    axx[3].plot(dat.LONGITUDE,layer['upper']['pden'])
    axx[3].set_ylabel('$\sigma_0$ [kg/m$^3$]')
    axx[3].set_xlim(-45,-5)
    axx[3].axhline(27.53,color='k')
    axx[3].set_xlabel('Longitude [$^\circ$W]')

    for ii in range(4):
        # axx[ii].grid('on')
        axx[ii].axvline(NA_wbnd,color='C1')
    savefig(figdir+'Freshwater/Layer_properties.png',bbox_inches='tight',dpi=300)

twolayer_timemean()


## Get  transport weighted North Atlantic water T,S,
tweight={}
tweight['AW']={}
tweight['AW']['sal']=(layer['upper']['sal']*layer['upper']['xport']).where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')/layer['upper']['xport'].where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')
tweight['AW']['tmp']=(layer['upper']['tmp']*layer['upper']['xport']).where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')/layer['upper']['xport'].where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')
tweight['AW']['pden']=(layer['upper']['pden']*layer['upper']['xport']).where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')/layer['upper']['xport'].where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')
tweight['AW']['xport']=layer['upper']['xport'].where(dat.LONGITUDE>NA_wbnd).sum(dim='LONGITUDE')


tweight['UW']={}
tweight['UW']['sal']=(layer['upper']['sal']*layer['upper']['xport']).sum(dim='LONGITUDE')/layer['upper']['xport'].sum(dim='LONGITUDE')
tweight['UW']['tmp']=(layer['upper']['tmp']*layer['upper']['xport']).sum(dim='LONGITUDE')/layer['upper']['xport'].sum(dim='LONGITUDE')
tweight['UW']['pden']=(layer['upper']['pden']*layer['upper']['xport']).sum(dim='LONGITUDE')/layer['upper']['xport'].sum(dim='LONGITUDE')
tweight['UW']['xport']=layer['upper']['xport'].sum(dim='LONGITUDE')


## Transport weighted polar water
tweight['PW']={}
tweight['PW']['sal']=(layer['upper']['sal']*layer['upper']['xport']).where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')/layer['upper']['xport'].where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')
tweight['PW']['tmp']=(layer['upper']['tmp']*layer['upper']['xport']).where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')/layer['upper']['xport'].where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')
tweight['PW']['pden']=(layer['upper']['pden']*layer['upper']['xport']).where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')/layer['upper']['xport'].where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')
tweight['PW']['xport']=layer['upper']['xport'].where(dat.LONGITUDE<=NA_wbnd).sum(dim='LONGITUDE')

## Transport weighted lower limb (deep water)
tweight['DW']={}
tweight['DW']['sal']=(layer['lower']['sal']*layer['lower']['xport']).sum(dim='LONGITUDE')/layer['lower']['xport'].sum(dim='LONGITUDE')
tweight['DW']['tmp']=(layer['lower']['tmp']*layer['lower']['xport']).sum(dim='LONGITUDE')/layer['lower']['xport'].sum(dim='LONGITUDE')
tweight['DW']['pden']=(layer['lower']['pden']*layer['lower']['xport']).sum(dim='LONGITUDE')/layer['lower']['xport'].sum(dim='LONGITUDE')
tweight['DW']['xport']=layer['lower']['xport'].sum(dim='LONGITUDE')

wmvec=['AW','PW','DW','UW']

salind=salvec>34

def plot_TS_bylayer():
    figure(figsize=(9,5.5))
    scatter([tweight[wm]['sal'] for wm in wmvec],[tweight[wm]['tmp'] for wm in wmvec],s=[abs(tweight[wm]['xport'].values)**2*4 for wm in wmvec],c=['r','b','b','purple'],zorder=50,linewidth=3)
    den=contour(salvec,tmpvec,pdenmat,levels=arange(24.13,29,0.2),colors='grey',zorder=1)
    contour(salvec[salind],tmpvec,pdenmat[:,salind],levels=[denlim],colors='k',zorder=100,linewidth=3)
    text(34.5,4.75,'PW ($S_2,T_2$)',color='k',fontsize=18)
    text(35.075,9.35,'AW ($S_1,T_1$)',color='k',fontsize=18)
    text(35.35,10.4,'upper',color='k',fontsize=18)
    text(35.01,4.3,'DW ($S_3,T_3$)',color='k',fontsize=18)
    xlim(34.45,35.5)
    ylim(2,12)
    xlabel('Salinity')
    ylabel('Theta [$^\circ$C]')
    title('Overturning transport in T-S space: OSNAP EAST \n',fontsize=14)
    savefig(figdir+'Freshwater/TS_transportweighted_2layerONLY.png',bbox_inches='tight',dpi=300)

plot_TS_bylayer()

def plot_TS_bylayer():
    figure(figsize=(9,5))
    scatter(layer['upper']['sal'],layer['upper']['tmp'],s=layer['upper']['xport']**2*60,c=layer['upper']['xport'],vmin=-1.5,vmax=1.5,cmap=cm.RdBu_r,zorder=20)
    scatter(layer['lower']['sal'],layer['lower']['tmp'],s=layer['lower']['xport']**2*60,c=layer['lower']['xport'],vmin=-1.5,vmax=1.5,cmap=cm.RdBu_r,zorder=20)
    colorbar(label='Transport per grid point within layer [Sv]')
    scatter([tweight[wm]['sal'] for wm in wmvec],[tweight[wm]['tmp'] for wm in wmvec],s=[abs(tweight[wm]['xport'].values)**2*4 for wm in wmvec],edgecolors=['r','b','b','purple'],zorder=50,facecolors='none',linewidth=3)
    den=contour(salvec,tmpvec,pdenmat,levels=arange(24.13,29,0.2),colors='grey',zorder=1)
    contour(salvec[salind],tmpvec,pdenmat[:,salind],levels=[denlim],colors='k',zorder=100,linewidth=3)
    text(34.5,4.75,'PW ($S_2,T_2$)',color='k',fontsize=18)
    text(35.075,9.35,'AW ($S_1,T_1$)',color='k',fontsize=18)
    text(35.35,10.4,'upper',color='k',fontsize=18)
    text(35.05,4.25,'DW ($S_3,T_3$)',color='k',fontsize=18)
    xlim(34.45,35.5)
    ylim(2,12)
    xlabel('Salinity')
    ylabel('Theta [$^\circ$C]')
    title('Overturning transport in T-S space: OSNAP EAST \n',fontsize=14)
    savefig(figdir+'Freshwater/TS_transportweighted_2layer.png',bbox_inches='tight',dpi=300)

plot_TS_bylayer()
xport_corr=tweight['UW']['xport']+tweight['DW']['xport']

xport_corr

## Apply a crude mass balance for now

(tweight['AW']['xport']-xport_corr/2)

tweight['PW']['xport']

(tweight['DW']['xport']-xport_corr/2)

q_S=tweight['AW']['sal']*(tweight['AW']['xport']-xport_corr/2)+tweight['PW']['sal']*tweight['PW']['xport']+tweight['DW']['sal']*(tweight['DW']['xport']-xport_corr/2)



q_T=tweight['AW']['tmp']*(tweight['AW']['xport']-xport_corr/2)+tweight['PW']['tmp']*tweight['PW']['xport']+tweight['DW']['tmp']*(tweight['DW']['xport']-xport_corr/2)

q_S

q_S/(tweight['AW']['xport']-xport_corr/2)

q_S*0.029

q_T
q_T/(tweight['AW']['xport']-xport_corr/2)

q_T*4.1

q_S_test=tweight['DW']['sal']*tweight['DW']['xport']+tweight['UW']['sal']*tweight['UW']['xport']
q_T_test=tweight['DW']['tmp']*tweight['DW']['xport']+tweight['UW']['tmp']*tweight['UW']['xport']



for key in tweight:
    print(key,tweight[key]['xport'])


## testing whether FWT can be recovered from this bulk view:
FW_bulk=tweight['UW']['xport']*(tweight['UW']['sal']-sref_east)/sref_east+tweight['DW']['xport']*(tweight['DW']['sal']-sref_east)/sref_east


sal_cols = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(sal_cols,position=[0,0.6,0.69,0.9,1],bit=True)
slevs=array([34,34.2,34.4,34.6,34.8,34.85,34.9,34.95,35,35.05,35.1,35.15,35.2,35.25,35.3,35.35,35.4])

def sal_fwt_denspace():
    figure()
    contourf(dat.LONGITUDE,den_lev,den_space.sal.mean(dim='date').T,cmap=sal_cmap,levels=slevs,extend='both')
    colorbar()
    contour(dat.LONGITUDE,den_lev,den_space.sal.mean(dim='date').T,levels=[sref_all],colors='k')
    ylim(25,28)
    xlabel('Longitude')
    ylabel('$\sigma_0$ [kg/m$^3$]')
    title('Salinity in density-longitude space')
    savefig(figdir+'Freshwater/Denlon_sal.png',bbox_inches='tight')

    figure()
    pcolor(dat.LONGITUDE,den_lev,den_space.xport.mean(dim='date').T,cmap=cm.RdBu_r,vmin=-0.1,vmax=0.1)
    colorbar()
    ylim(25,28)
    xlabel('Longitude')
    ylabel('$\sigma_0$ [kg/m$^3$]')
    title('Volume transport in density-longitude space')
    savefig(figdir+'Freshwater/Denlon_xport.png',bbox_inches='tight')

    figure()
    pcolor(dat.LONGITUDE,den_lev,den_space.fwt.mean(dim='date').T,cmap=cm.RdBu_r,vmin=-0.5,vmax=0.5)
    colorbar()
    ylim(25,28)
    xlabel('Longitude')
    ylabel('$\sigma_0$ [kg/m$^3$]')
    title('Freshwater transport in density-longitude space')
    savefig(figdir+'Freshwater/Denlon_fwt.png',bbox_inches='tight')

sal_fwt_denspace()


def plot_den_psisal(minlon,maxlon,tit):
    f,axx=subplots(1,3,sharey=True,figsize=(12,8))
    axx[0].plot(den_space.xport.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum('LONGITUDE').T,den_lev,color='grey')
    axx[0].plot(den_space.xport.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).mean(dim='date').sum('LONGITUDE').T,den_lev,color='k',linewidth=3)
    axx[0].set_ylabel('$\sigma_0$ [kg/m$^3$]')
    axx[0].set_xlabel('Volume transport [Sv]')

    axx[1].plot(den_space.fwt.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum('LONGITUDE').T,den_lev,color='grey')
    axx[1].plot(den_space.fwt.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).mean(dim='date').sum('LONGITUDE').T,den_lev,color='k',linewidth=3)
    axx[1].set_ylim(28,24)
    axx[1].set_xlabel('Freshwater transport [mSv]')

    axx[1].set_title('Freshwater transport in density space: '+tit)

    axx[2].plot(den_space.fwt.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum('LONGITUDE').cumsum(dim='den_lev').T,den_lev,color='grey')
    axx[2].plot(den_space.fwt.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).mean(dim='date').sum('LONGITUDE').cumsum(dim='den_lev').T,den_lev,color='k',linewidth=3)
    axx[2].set_xlabel('Cumulative freshwater transport [mSv]')
    savefig(figdir+'Freshwater/FWT_denspace_'+tit+'.png',bbox_inches='tight')


plot_den_psisal(-70,0,'All')
plot_den_psisal(-70,-45,'West')
plot_den_psisal(-45,0,'East')


def plot_sal_psi(minlon,maxlon,tit):
    f,axx=subplots(1,3,sharey=True,figsize=(12,8))

    axx[0].plot(sal_space.xport.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum(dim='LONGITUDE').T,sal_lev,color='grey');
    axx[0].plot(sal_space.xport.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum(dim='LONGITUDE').mean(dim='date').T,sal_lev,color='k',linewidth=3);
    axx[0].axhline(sref_all,color='k',linestyle='--')
    axx[0].set_ylabel('Salinity')
    axx[0].set_xlabel('Volume transport per bin [Sv]')


    axx[1].plot(sal_space.xport.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum(dim='LONGITUDE').cumsum(dim='sal_lev').T,sal_lev,color='grey');
    axx[1].plot(sal_space.xport.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum(dim='LONGITUDE').cumsum(dim='sal_lev').mean(dim='date').T,sal_lev,color='k',linewidth=3);
    axx[1].axhline(sref_all,color='k',linestyle='--')
    axx[1].set_xlabel('Cumulative transport [Sv]')

    axx[1].set_title('Transports in salinity space: '+tit)

    axx[2].plot(sal_space.fwt.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum(dim='LONGITUDE').cumsum(dim='sal_lev').T,sal_lev,color='grey');
    axx[2].plot(sal_space.fwt.where(sal_space.lon>=minlon).where(sal_space.lon<=maxlon).sum(dim='LONGITUDE').cumsum(dim='sal_lev').mean(dim='date').T,sal_lev,color='k',linewidth=3);
    axx[2].axhline(sref_all,color='k',linestyle='--')
    axx[2].set_ylim(31,35.5)
    axx[2].set_xlabel('Cumulative freshwater transport [mSv]')

    savefig(figdir+'Freshwater/Psi_salspace_'+tit+'.png',bbox_inches='tight')


plot_sal_psi(-70,0,'All')
plot_sal_psi(-70,-45,'West')
plot_sal_psi(-45,0,'East')

#good: check's out, fwt is the same in every space...
# def check_fwt():
#     den_space.fwt.sum(dim='LONGITUDE').sum(dim='den_lev').plot()
#     sal_space.fwt.sum(dim='LONGITUDE').sum(dim='sal_lev').plot()
#     dat.FWT.sum(dim='LONGITUDE').sum(dim='DEPTH').plot()
#
#
# check_fwt()
