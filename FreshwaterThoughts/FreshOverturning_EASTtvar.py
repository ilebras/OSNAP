from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_gridded/'

#############################################
### Load data, establish some parameters
#############################################

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))
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


######### Get isopycnal of maximum overturning varying in time:
MOC_east=ones(21)
SIGmax_east=ones(21)
for ii in range(21):
    MOC_east[ii]=psi.T_EAST[ii,:].max()
    SIGmax_east[ii]=psi.LEVEL[psi.T_EAST[ii,:].argmax()]



for tt in range(21):
    figure()
    contourf(east.LONGITUDE,east.DEPTH,east.VELO[tt,:,:],20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
    contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],20,levels=[SIGmax_east[tt]],linewidths=3,colors='k')
    contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],20,levels=[27.53],linewidths=3,colors='r')
    ylim(2000,-100)
    title(str(east.TIME[tt].values)[:10]+', $\sigma_\Theta$ = '+str(SIGmax_east[tt])[:5])
    xlim(-44,-5)


lonbnd=-44
######### Upper and lower layer breakdown varying in time:

SIGmax_east



east=dat.where(dat.LONGITUDE>lonbnd)

upper=NaN*east
lower=NaN*east
for tt in range(21):
    for var in ['TRANS','PTMP','PSAL','PDEN']:
        upper[var][tt,:,:]=east[var][tt,:,:].where(east.PDEN[tt,:,:]<SIGmax_east[tt])
        lower[var][tt,:,:]=east[var][tt,:,:].where(east.PDEN[tt,:,:]>=SIGmax_east[tt])

upper.TRANS.sum(dim='DEPTH')[:,::-1]

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

upper.TRANS.sum(dim='DEPTH').plot()
lower.TRANS.sum(dim='DEPTH').plot()

NA_wbnd=-40


def transweighted_props(xray):
    xray_out={}
    xray_out['TRANS']=xray['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
    for var in ['PSAL','PTMP','PDEN']:
        xray_out[var]=(xray[var]*xray['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')

    return xray_out


UW=transweighted_props(upper)
AW=transweighted_props(upper.where(dat.LONGITUDE>NA_wbnd))
PW=transweighted_props(upper.where(dat.LONGITUDE<=NA_wbnd))
DW=transweighted_props(lower)



wmvec=[AW,PW,DW,UW]
plot(psi.T_EAST[0,:],psi.LEVEL)

def plot_TS_bylayer():
    for tt in range(2):
        f,[ax3,ax2,ax1]=subplots(1,3,figsize=(16,5))
        ax3.plot(psi.T_EAST[tt,:],psi.LEVEL,'k')
        ax3.set_ylim(28,26)
        ax3.set_xlim(-5,20)
        ax3.set_ylabel('$\sigma_\Theta$')
        ax3.set_xlabel('Transport [Sv]')
        vel=ax1.contourf(east.LONGITUDE,east.DEPTH,east.VELO[tt,:,:],20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
        ax1.contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],20,levels=[SIGmax_east[tt]],linewidths=3,colors='k')
        ax1.contour(east.LONGITUDE,east.DEPTH,east.PDEN[tt,:,:],20,levels=[27.53],linewidths=3,colors='r')
        ax1.set_ylim(2000,-100)
        ax1.set_xlim(-44,-5)
        ax1.set_xlabel('depth [m]')
        ax1.set_ylabel('longitude')
        ax1.axvline(-40,color='grey')
        colorbar(vel,label='velocity [m/s]')

        ax2.scatter([wm['PSAL'][tt] for wm in wmvec],[wm['PTMP'][tt] for wm in wmvec],s=[abs(wm['TRANS'][tt].values)**2*4 for wm in wmvec],c=['r','c','b','purple'],zorder=50,linewidth=3)
        den=ax2.contour(salvec,tmpvec,pdenmat,levels=arange(24.13,29,0.2),colors='grey',zorder=1)
        ax2.contour(salvec[salind],tmpvec,pdenmat[:,salind],levels=[SIGmax_east[tt]],colors='k',zorder=100,linewidth=3)

        ax2.set_xlim(33.75,35.75)
        ax2.set_ylim(0,12)
        ax2.set_xlabel('Salinity')
        ax2.set_ylabel('Theta [$^\circ$C]')

        f.suptitle(str(east.TIME[tt].values)[:10]+', MOC = '+str(MOC_east[tt])[:5]+', $\sigma_\Theta$ = '+str(SIGmax_east[tt])[:5])
        savefig(figdir+'Freshwater/east_tvar/LON_TS_decomp_'+str(tt)+'.png',bbox_inches='tight',dpi=300)

plot_TS_bylayer()

def comp_var(var):
    figure(figsize=(6,5))
    AW[var].plot(color='red')
    DW[var].plot(color='darkblue')
    PW[var].plot(color='cyan')

comp_var('TRANS')

comp_var('PDEN')
plot(dat.TIME,SIGmax_east,'k')

comp_var('PTMP')

comp_var('PSAL')
