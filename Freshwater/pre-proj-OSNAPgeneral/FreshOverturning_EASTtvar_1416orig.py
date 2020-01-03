from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/'


#############################################
### Background T,S, density for plotting on TS diagrams
#############################################

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

#############################################
# Retiring this for now because I want to calculate it all meself so its clearer what I'm doing.
Feili_psi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Streamfunction_201408_201604_2018.nc')
denvec=Feili_psi.LEVEL.values
# ######### Get isopycnal of maximum overturning varying in time:
# MOC_east=ones(21)
# SIGmax_east=ones(21)
# for ii in range(21):
#     MOC_east[ii]=psi.T_EAST[ii,:].max()
#     SIGmax_east[ii]=psi.LEVEL[psi.T_EAST[ii,:].argmax()]
#############################################

#############################################
### Load data, establish some parameters
#############################################
# dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','r'))
dat=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.nc')


east.LATITUDE[dat.LONGITUDE.values>lonbnd].values

east.LONGITUDE[dat.LONGITUDE>lonbnd].values

east['TRANS']=east['VELO']*east['AREA']

east.TRANS.sum(dim='DEPTH').sum(dim='LONGITUDE')
vel_corr=float(east.TRANS.sum(dim='DEPTH').sum(dim='LONGITUDE').mean(dim='TIME')/east.AREA.sum(dim='DEPTH').sum(dim='LONGITUDE'))
vel_corr

east['VELO_corr']=east['VELO']-vel_corr

east['TRANS_corr']=east['VELO_corr']*east.AREA

east['TRANS_corr'].sum(dim='DEPTH').sum(dim='LONGITUDE')
dat
#############################################
### Calculate overturning in density space (and corresponding density of max overturning) from two transport estimates
#############################################
denstep=0.01

east=east.assign_coords(DENLEV=(denvec))

east

def getpsi(var,newname):
    psimat=NaN*ones((len(denvec),len(east.TIME)))
    for ii,dd in enumerate(denvec):
            psimat[ii,:]=east[var].where(east['PDEN']>dd-denstep/2).where(east['PDEN']<=dd+denstep/2).sum(dim='DEPTH').sum(dim='LONGITUDE').values
    east[newname]=((['DENLEV','TIME']),psimat)
    return east


east=getpsi('TRANS','PSI')
east=getpsi('TRANS_corr','PSI_corr')

east['MOC']=east.PSI.cumsum(dim='DENLEV').max(dim='DENLEV')/1e6
east['MOC_corr']=east.PSI_corr.cumsum(dim='DENLEV').max(dim='DENLEV')/1e6

east.MOC.plot()
east['MOC_corr'].plot()
east['SIGMAX']=(('TIME'),east.DENLEV[east.PSI.cumsum(dim='DENLEV').argmax(dim='DENLEV').values])
east['SIGMAX_corr']=(('TIME'),east.DENLEV[east.PSI_corr.cumsum(dim='DENLEV').argmax(dim='DENLEV').values])

east.SIGMAX.plot()
east['SIGMAX_corr'].plot()

#No real difference, seems like it doesn't matter... just use the _corr ones

#############################################
######### Upper and lower layer breakdown with sigma_max varying in time.
#############################################

upper=xr.Dataset()
lower=xr.Dataset()
for var in ['TRANS_corr','PTMP','PSAL','PDEN']:
    if var=='TRANS_corr':
        var2='TRANS'
    else:
        var2=var
    upper[var2]=east[var].where(east.PDEN<east.SIGMAX_corr)
    lower[var2]=east[var].where(east.PDEN>=east.SIGMAX_corr)



def plot_londep():
    f,[(ax1,ax2),(ax3,ax4),(ax5,ax6)]=subplots(3,2,figsize=(20,6),sharex=True)
    ax1.plot(east.LONGITUDE,upper.TRANS.sum(dim='DEPTH').T)
    ax1.set_title('Upper Layer')
    ax2.plot(east.LONGITUDE,lower.TRANS.sum(dim='DEPTH').T)
    ax2.set_title('Lower Layer')
    ax2.set_xlim(-44,-5)
    ax3.plot(east.LONGITUDE,upper.TRANS.sum(dim='DEPTH').cumsum(dim='LONGITUDE').T)
    ax4.plot(east.LONGITUDE,lower.TRANS.sum(dim='DEPTH').cumsum(dim='LONGITUDE').T)
    ax5.plot(east.LONGITUDE[::-1],upper.TRANS.sum(dim='DEPTH')[:,::-1].cumsum(dim='LONGITUDE').T)
    ax6.plot(east.LONGITUDE[::-1],lower.TRANS.sum(dim='DEPTH')[:,::-1].cumsum(dim='LONGITUDE').T)

plot_londep()

def get_uplow_const(const):
    upper_mean=xr.Dataset()
    lower_mean=xr.Dataset()
    for var in ['TRANS_corr','PTMP','PSAL','PDEN']:
        if var=='TRANS_corr':
            var2='TRANS'
        else:
            var2=var
        upper_mean[var2]=east[var].where(east.PDEN<const)
        lower_mean[var2]=east[var].where(east.PDEN>=const)
    return upper_mean,lower_mean

upper_eastm,lower_eastm=get_uplow_const(27.53)
upper_fullm,lower_fullm=get_uplow_const(27.66)

def transweighted_props(xray):
    xray_out={}
    xray_out['TRANS']=xray['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')/1e6
    for var in ['PSAL','PTMP','PDEN']:
        xray_out[var]=(xray[var]*xray['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
    return xray_out

UW={}
AW={}
PW={}
DW={}

NA_wbnd=-41

def get_WM(vers,upvar,lowvar):
    UW[vers]=transweighted_props(upvar)
    AW[vers]=transweighted_props(upvar.where(dat.LONGITUDE>NA_wbnd))
    PW[vers]=transweighted_props(upvar.where(dat.LONGITUDE<=NA_wbnd))
    DW[vers]=transweighted_props(lowvar)
    return UW,AW,PW,DW

UW,AW,PW,DW=get_WM('tvar',upper,lower)
UW,AW,PW,DW=get_WM('eastm',upper_eastm,lower_eastm)
UW,AW,PW,DW=get_WM('fullm',upper_fullm,lower_fullm)

def plot_WMprops(vers,tit):
    f,axx=subplots(4,1,figsize=(8,20))
    for ii,var in enumerate(['TRANS','PDEN','PSAL','PTMP']):
        AW[vers][var].plot(color='red',ax=axx[ii],label='Atlantic Water')
        PW[vers][var].plot(color='C0',ax=axx[ii],label='Polar Water')
        # UW[vers][var].plot(color='purple',ax=axx[ii])
        DW[vers][var].plot(color='k',ax=axx[ii],label='Deep Water')
        axx[ii].set_ylabel(var)
        axx[ii].set_xlabel('')
    axx[0].legend()
    axx[0].set_title(tit)
    savefig(figdir+'WMdecomp_tvar/OSNAPEast_decomp/PropDecomp_'+vers+'.png',bbox_inches='tight')

plot_WMprops('tvar','Time varying isopycnal of maximum overturning')

plot_WMprops('fullm','Using OSNAP-wide mean isopycnal of maximum overturning')

plot_WMprops('eastm','Using OSNAP East mean isopycnal of maximum overturning')




def plot_WMprops_comp(versvec,linevec):
    lsi=2
    f,axx=subplots(4,1,figsize=(8,20))
    for ii,var in enumerate(['TRANS','PDEN','PSAL','PTMP']):
        for ll,vers in enumerate(versvec):
            AW[vers][var].plot(color='red',ax=axx[ii],linestyle=linevec[ll],linewidth=lsi)
            PW[vers][var].plot(color='C0',ax=axx[ii],linestyle=linevec[ll],linewidth=lsi)
            # UW[vers][var].plot(color='purple',ax=axx[ii],linestyle=linevec[ll],linewidth=lsi)
            DW[vers][var].plot(color='k',ax=axx[ii],linestyle=linevec[ll],linewidth=lsi)
            axx[ii].set_ylabel(var)
            axx[ii].set_xlabel('')

PW['fullm']['PSAL'].plot()
PW['eastm']['PSAL'].plot()

for dic in [AW, DW, PW]:
    print(dic['eastm']['TRANS'].mean().values)
    print(dic['eastm']['PSAL'].mean().values)

PW['fullm']['TRANS'].plot()
PW['eastm']['TRANS'].plot()

plot_WMprops_comp(['tvar','eastm',],['--','-','-.'])


(-AW['fullm']['TRANS']).plot()
(PW['fullm']['TRANS']).plot()
(DW['fullm']['TRANS']).plot()

#############################################
### Look at TS props with density bounds, choose a salinity bound(?)
#############################################

lonbigmat=tile(east.LONGITUDE.values,[len(east.DEPTH),1])

shape(east.PSAL)

lonmat=tile(east.LONGITUDE.values,[len(east.DEPTH),1])
lonbigmat=tile(lonmat,[len(east.TIME),1,1])
shape(lonbigmat)


def plot_TS_withdefs():
        # figure()
        # scatter(east.PSAL.values.flatten(),east.PTMP.values.flatten(),c=lonbigmat.flatten(),alpha=0.3)
        # colorbar(label='longitude')
        # contour(salvec,tmpvec,pdenmat,levels=denvec[::30],colors='k',zorder=3)
        # contour(salvec,tmpvec,pdenmat,levels=[27.53],colors='r',linewidths=3,zorder=4) #A bit skeptical of this as these are all densities calculated as though the properties were at surface
        # contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='g',linewidths=3,zorder=4)
        # xlim(32,36)
        # ylim(-2,15)
        # [axvline(kk,color='k') for kk in [34.8]]
        # xlabel('salinity')
        # ylabel('pot. temperature [$^\circ$C]')
        # savefig(figdir+'WMdecomp_tvar/TS/TSmean_wbnds.png',bbox_inches='tight')
        # xlim(34.5,35.1)
        # ylim(0,6)
        # savefig(figdir+'WMdecomp_tvar/TS/TSmean_wbnds_zoom.png',bbox_inches='tight')

        for tt in range(21):
            figure()
            scatter(east.PSAL[tt,:,:],east.PTMP.values[tt,:,:],c=lonmat,alpha=0.3)
            colorbar()
            contour(salvec,tmpvec,pdenmat,levels=denvec[::30],colors='k',zorder=3)
            contour(salvec,tmpvec,pdenmat,levels=[float(east.SIGMAX_corr[tt])],colors='orange',zorder=3,linewidths=3)
            contour(salvec,tmpvec,pdenmat,levels=[27.53],colors='r',linewidths=3,zorder=4) #A bit skeptical of this as these are all densities calculated as though the properties were at surface
            contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='g',linewidths=3,zorder=4)
            xlim(32,36)
            ylim(-2,15)
            [axvline(kk,color='k') for kk in [34.8]]
            title(str(east.TIME[tt].values)[:10])
            xlabel('salinity')
            ylabel('pot. temperature [$^\circ$C]')
            savefig(figdir+'WMdecomp_tvar/TS/TS_wbnds_'+str(tt)+'.png',bbox_inches='tight')
            xlim(34.5,35.1)
            ylim(0,6)
            savefig(figdir+'WMdecomp_tvar/TS/TS_wbnds_zoom_'+str(tt)+'.png',bbox_inches='tight')


plot_TS_withdefs()

sal_cols = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(sal_cols,position=[0,0.6,0.69,0.9,1],bit=True)
slevs=array([34,34.2,34.4,34.6,34.8,34.85,34.9,34.95,35,35.05,35.1,35.15,35.2,35.25,35.3,35.35,35.4])

for tt in range(21):
    figure()
    contourf(east.LONGITUDE,east.DEPTH,east.VELO[tt,:,:],20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
    colorbar(label='cross-track velocity [m/s]')
    contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],levels=[float(east.SIGMAX_corr[tt])],linewidths=3,colors='orange')
    contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],levels=[27.53],linewidths=3,colors='r')
    contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],levels=[27.66],linewidths=3,colors='g')
    ylim(2000,-100)
    title(str(east.TIME[tt].values)[:10]+', $\sigma_\Theta$ = '+str(float(east.SIGMAX_corr[tt]))[:5])
    xlim(-44,-5)
    xlabel('longitude')
    ylabel('depth [m]')
    savefig(figdir+'WMdecomp_tvar/sections/vel_wsigmas_'+str(tt)+'.png',bbox_inches='tight')

for tt in range(21):
    figure(figsize=(20,10))
    contourf(east.LONGITUDE,east.DEPTH,east.PSAL[tt,:,:],cmap=sal_cmap,levels=slevs,extend='both')
    colorbar(label='salinity')
    contour(east.LONGITUDE,east.DEPTH,east.PSAL[tt,:,:],[34.8],colors='g',linewidths=3)
    contour(east.LONGITUDE,east.DEPTH,east.PSAL[tt,:,:],[34.9],colors='k',linewidths=3)
    contour(east.LONGITUDE,east.DEPTH,east.PSAL[tt,:,:],[34.95],colors='purple',linewidths=3)
    contour(east.LONGITUDE,east.DEPTH,east.PSAL[tt,:,:],[35],colors='grey',linewidths=3)
    # contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],levels=[float(east.SIGMAX_corr[tt])],linewidths=3,colors='orange')
    # contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],levels=[27.53],linewidths=3,colors='r')
    # contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],levels=[27.66],linewidths=3,colors='g')
    ylim(2000,-100)
    title(str(east.TIME[tt].values)[:10]+', $\sigma_\Theta$ = '+str(float(east.SIGMAX_corr[tt]))[:5])
    xlim(-44,-5)
    xlabel('longitude')
    ylabel('depth [m]')
    savefig(figdir+'WMdecomp_tvar/sections/sal_wsals_'+str(tt)+'.png',bbox_inches='tight')

####### Try splitting up by salinity instead:
# AW = saltier than 35
# PW = fresher than 34.8 (maybe split off the Scotland shelf)
# DW  = everything in between

vers='salbnd'
AW[vers]=transweighted_props(east.where(east.PSAL>35))
PW[vers]=transweighted_props(east.where(east.PSAL<=34.8))
DW[vers]=transweighted_props(east.where(east.PSAL>34.8).where(east.PSAL<=35))

plot_WMprops('salbnd')
