from pylab import *
import xarray as xr
from aux_funcs import *

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

eastind=dat.LONGITUDE>-45
westind=dat.LONGITUDE<-45
allind=dat.LONGITUDE<0
CFind=(dat.LONGITUDE>-45)&(dat.LONGITUDE<=-40)

sref_west=34.85
sref_east=34.97
sref_all=34.92

dat['TRANS']=dat.VELO*dat.AREA/1e6

etrans_corr=dat['TRANS'][:,:,eastind].sum(dim='DEPTH').sum(dim='LONGITUDE')*1e6/dat['AREA'][:,eastind].sum(dim='DEPTH').sum(dim='LONGITUDE')
dat['TRANS_east']=(dat.VELO[:,:,eastind]-etrans_corr)*dat.AREA[:,eastind]/1e6

dat['FWT']=dat['TRANS']*(sref_all-dat['PSAL'])/sref_all*1e3
dat['FWT_east']=(dat.VELO[:,:,eastind]-etrans_corr)*dat.AREA[:,eastind]*(sref_east-dat['PSAL'][:,:,eastind])/sref_east/1e3

# (dat['TRANS']*1e3).sum(dim='DEPTH').sum(dim='LONGITUDE').plot(label='all')
# (dat['TRANS_east']*1e3).sum(dim='DEPTH').sum(dim='LONGITUDE').plot(label='east')
# legend()
#
#
# (dat['FWT_east']).sum(dim='DEPTH').sum(dim='LONGITUDE').plot(label='all')
# (dat['FWT_east'][:,:,CFind]).sum(dim='DEPTH').sum(dim='LONGITUDE').plot(label='east')
# legend()

# (((dat['FWT_east'][:,:,CFind]).sum(dim='DEPTH').sum(dim='LONGITUDE'))/((dat['FWT_east']).sum(dim='DEPTH').sum(dim='LONGITUDE'))).plot()
East_FWT=dat['FWT_east'].sum(dim='DEPTH').sum(dim='LONGITUDE').mean()
CF_FWT=(dat['FWT_east'][:,:,CFind]).sum(dim='DEPTH').sum(dim='LONGITUDE').mean()

sal_cols = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(sal_cols,position=[0,0.6,0.69,0.9,1],bit=True)
slevs=array([34,34.2,34.4,34.6,34.8,34.85,34.9,34.95,35,35.05,35.1,35.15,35.2,35.25,35.3,35.35,35.4])

def plotsal_east():
    f,axx=subplots(3,1,sharex=True,figsize=(12,9),gridspec_kw = {'height_ratios':[1,1.5,1.5]})

    axx[0].plot(dat.LONGITUDE,dat['FWT_east'].sum(dim='DEPTH').mean(dim='TIME').cumsum(dim='LONGITUDE'),color='k',linewidth=2)
    # axx[0].plot(CFlon[-1],-96,'ro')
    # axx[0].text(-42,-120,'direct estimate',color='red',fontsize=14)
    axx[0].plot(-40,CF_FWT,'ko')
    axx[0].text(-41,CF_FWT+25,str(int(CF_FWT))+' mSv',color='k',fontsize='16')
    axx[0].plot(-7.1,East_FWT,'ko')
    axx[0].text(-6.9,-170,str(int(East_FWT))+' mSv',color='k',fontsize='16')
    # axx[0].text(-35,-140,'East Greenland Current System carries half the OSNAP EAST freshwater transport',color='black')
    axx[0].set_ylabel('[mSv]')
    # axx[0].set_title('OSNAP EAST cumulative freshwater transport= $\int_{x_{w}} ^x v\ A \ (\overline{S}-S)\ /\ \overline{S}\ dx$ \n August 2014-April 2016 mean\n',fontsize=25)

    cxx=0.92
    cy=0.25
    cax1=f.add_axes([cxx,0.425,0.02,cy])
    vel=axx[1].contourf(dat.LONGITUDE,dat.DEPTH,dat.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
    axx[1].contour(dat.LONGITUDE,dat.DEPTH,dat.VELO.mean(dim='TIME'),levels=[0],colors='k')
    axx[1].set_ylim(3500,0)
    axx[1].set_ylabel('depth [m]')
    colorbar(vel,ax=axx[1],cax=cax1,label='Velocity [m/s]')
    axx[1].fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
    axx[1].text(-44,3200,'Greenland',color='white',fontsize=16)
    axx[1].text(-34.5,3200,'Mid-Atlantic Ridge',color='white',fontsize=16)
    axx[1].text(-11,3200,'Scotland',color='white',fontsize=16)
    # axx[1].text(-19,2800,'Mean 2014-2016',color='white',fontsize=15)

    cax2=f.add_axes([cxx,0.125,0.02,cy])
    sal=axx[2].contourf(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),slevs,cmap=sal_cmap,extend='both')
    c1=axx[2].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[34.9],colors='grey')
    axx[2].set_ylabel('depth [m]')
    axx[2].fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
    clabel(c1,fmt='%1.1f')
    c2=axx[2].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[sref_east],colors='k')
    clabel(c2,fmt='%1.2f')
    axx[2].set_ylim(3500,0)
    axx[2].set_xlabel('Longitude')
    colorbar(sal,ax=axx[2],ticks=slevs[::2],cax=cax2,label='Salinity')
    axx[0].set_xlim(-44,-7)
    savefig(figdir+'Freshwater/FWT_salvelsec_east.png',bbox_inches='tight')

plotsal_east()





def plotsal_all():
    f,axx=subplots(3,1,sharex=True,figsize=(15,9),gridspec_kw = {'height_ratios':[1,1.5,1.5]})

    axx[0].plot(dat.LONGITUDE,dat['FWT'].sum(dim='DEPTH').mean(dim='TIME').cumsum(dim='LONGITUDE'),color='k',linewidth=2)
    # axx[0].plot(CFlon[-1],-96,'ro')
    # axx[0].text(-42,-120,'direct estimate',color='red',fontsize=14)
    # axx[0].plot(-40,CF_FWT,'ko')
    # axx[0].text(-41,CF_FWT+25,str(int(CF_FWT))+' mSv',color='k',fontsize='16')
    # axx[0].plot(-7.1,East_FWT,'ko')
    # axx[0].text(-6.9,-170,str(int(East_FWT))+' mSv',color='k',fontsize='16')
    # axx[0].text(-35,-140,'East Greenland Current System carries half the OSNAP EAST freshwater transport',color='black')
    axx[0].set_ylabel('[mSv]')
    axx[0].grid('on')
    # axx[0].set_title('OSNAP EAST cumulative freshwater transport= $\int_{x_{w}} ^x v\ A \ (\overline{S}-S)\ /\ \overline{S}\ dx$ \n August 2014-April 2016 mean\n',fontsize=25)

    cxx=0.92
    cy=0.25
    cax1=f.add_axes([cxx,0.425,0.02,cy])
    vel=axx[1].contourf(dat.LONGITUDE,dat.DEPTH,dat.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
    axx[1].contour(dat.LONGITUDE,dat.DEPTH,dat.VELO.mean(dim='TIME'),levels=[0],colors='k')
    axx[1].set_ylim(3500,0)
    axx[1].set_ylabel('depth [m]')
    colorbar(vel,ax=axx[1],cax=cax1,label='Velocity [m/s]')
    axx[1].fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
    # axx[1].text(-44,3200,'Greenland',color='white',fontsize=16)
    # axx[1].text(-34.5,3200,'Mid-Atlantic Ridge',color='white',fontsize=16)
    # axx[1].text(-11,3200,'Scotland',color='white',fontsize=16)
    # axx[1].text(-19,2800,'Mean 2014-2016',color='white',fontsize=15)

    cax2=f.add_axes([cxx,0.125,0.02,cy])
    sal=axx[2].contourf(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),slevs,cmap=sal_cmap,extend='both')
    c1=axx[2].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[34.9],colors='grey')
    axx[2].set_ylabel('depth [m]')
    axx[2].fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
    clabel(c1,fmt='%1.1f')
    c2=axx[2].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[sref_all],colors='k')
    clabel(c2,fmt='%1.2f')
    axx[2].set_ylim(3500,0)
    axx[2].set_xlabel('Longitude')
    colorbar(sal,ax=axx[2],ticks=slevs[::2],cax=cax2,label='Salinity')
    savefig(figdir+'Freshwater/FWT_salvelsec_all.png',bbox_inches='tight')

plotsal_all()

def plotdensec_all():
        f,axx=subplots(1,1,figsize=(25,5))
        den=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat.PDEN.mean(dim='TIME'),univec['pden'][1],cmap=univec['pden'][2],extend='both')
        axx.contour(dat.LONGITUDE,dat.DEPTH,dat.PDEN.mean(dim='TIME'),levels=univec['pden'][1][::2],colors='k')
        axx.set_ylabel('depth [m]')
        axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
        axx.set_ylim(3750,0)
        axx.set_xlabel('Longitude')
        colorbar(den,label='$\sigma_0$ [kg/m$^{3}$]')
        title('21 month mean potential density across OSNAP')
        savefig(figdir+'Freshwater/pden/Pdensec_mean.png',bbox_inches='tight')

        for ii in range(21):

                f,axx=subplots(1,1,figsize=(25,15))
                den=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat.PDEN[ii,:,:],univec['pden'][1],cmap=univec['pden'][2],extend='both')
                axx.contour(dat.LONGITUDE,dat.DEPTH,dat.PDEN[ii,:,:],levels=univec['pden'][1][::2],colors='k')
                axx.set_ylabel('depth [m]')
                axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
                axx.set_ylim(3750,0)
                axx.set_xlabel('Longitude')
                colorbar(den,label='$\sigma_0$ [kg/m$^{3}$]')
                title('Potential density: '+str(dat.TIME[ii].values)[:10])
                savefig(figdir+'Freshwater/pden/Pdensec_'+str(ii).zfill(2)+'.png',bbox_inches='tight',dpi=300)

plotdensec_all()
univec['pden'][1]
diff(univec['pden'][1][::2])


trans_cols = [(37,52,148),(67,162,202) ,(255,255,255),(254,178,76),(240,59,32)]
trans_cmap = make_cmap(trans_cols,position=[0,0.45,0.5,0.55,1],bit=True)

d1=27.67
d2=27.73
d3=27.77

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


def TS_transmn():
    figure(figsize=(18,3))
    subplot(1,3,1)
    contour(salvec,tmpvec,pdenmat,levels=arange(27.06,28,0.2),colors='grey',zorder=1)
    ts=hexbin(dat.PSAL[:,:,eastind].mean(dim='TIME').values.flatten(),dat.PTMP[:,:,eastind].mean(dim='TIME').values.flatten(),
           C=(dat.VELO[:,:,eastind].mean(dim='TIME')*dat.AREA).values.flatten()/1e6,cmap=trans_cmap,vmin=-0.05,vmax=0.05,zorder=2)
    contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    colorbar(ts)
    axvline(sref_east,color='k')
    xlim(32,35.5)
    ylim(-2,12)
    title('EAST')
    ylabel('potential temperature')
    subplot(1,3,2)
    contour(salvec,tmpvec,pdenmat,levels=arange(27.06,28,0.2),colors='grey',zorder=1)
    ts=hexbin(dat.PSAL[:,:,westind].mean(dim='TIME').values.flatten(),dat.PTMP[:,:,westind].mean(dim='TIME').values.flatten(),
           C=(dat.VELO[:,:,westind].mean(dim='TIME')*dat.AREA).values.flatten()/1e6,cmap=trans_cmap,vmin=-0.05,vmax=0.05,zorder=2)
    contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    colorbar(ts)
    axvline(sref_east,color='k')
    xlim(32,35.5)
    ylim(-2,12)
    xlabel('Salinity')
    title('WEST')
    subplot(1,3,3)
    contour(salvec,tmpvec,pdenmat,levels=arange(27.06,28,0.2),colors='grey',zorder=1)
    ts=hexbin(dat.PSAL.mean(dim='TIME').values.flatten(),dat.PTMP.mean(dim='TIME').values.flatten(),
           C=(dat.VELO.mean(dim='TIME')*dat.AREA).values.flatten()/1e6,cmap=trans_cmap,vmin=-0.05,vmax=0.05,zorder=2)
    contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    # contour(salvec,tmpvec,pdenmat2,levels=[27.66],colors='r')
    colorbar(ts)
    axvline(sref_east,color='k')
    xlim(32,35.5)
    ylim(-2,12)
    title('ALL')
    savefig(figdir+'Freshwater/Xport_TS_timemean.png',bbox_inches='tight')
    # suptitle('Time-mean transport in TS space \n\n ',fontsize=16)

TS_transmn()


def TS_trans(ind,axx):
    axx.contour(salvec,tmpvec,pdenmat,levels=arange(24.06,28,0.4),colors='grey')
    ts=axx.hexbin(dat.PSAL[:,:,ind].values.flatten(),dat.PTMP[:,:,ind].values.flatten(),
           C=(dat.TRANS[:,:,ind]).values.flatten(),cmap=trans_cmap,gridsize=50,vmin=-0.075,vmax=0.075,)
    axx.contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    axx.axvline(sref_all,color='k')
    return ts


def TS_transplot():
    f,(ax1,ax2,ax3)=subplots(1,3,figsize=(14,3.5),sharex=True,sharey=True)
    ts=TS_trans(eastind,ax1)
    ax1.set_title('East')
    ts=TS_trans(westind,ax2)
    ax2.set_title('West')
    ts=TS_trans(allind,ax3)
    ax3.set_title('All')
    cbax=f.add_axes([0.95,0.1,0.02,0.8])
    colorbar(ts,cax=cbax,label='Transport [Sv]')
    ax1.set_xlim(31,36)
    ax1.set_ylim(-3,15)
    f.text(0.5, -0.05, 'Salinity', ha='center',fontsize=12)
    f.text(0.05, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=12)
    savefig(figdir+'Freshwater/Xport_TS_alltime.png',bbox_inches='tight')

TS_transplot()


def TS_trans_onetime(ind,axx,tchoose):
    axx.contour(salvec,tmpvec,pdenmat,levels=arange(24.06,28,0.4),colors='grey')
    ts=axx.hexbin(dat.PSAL[tchoose,:,ind].values.flatten(),dat.PTMP[tchoose,:,ind].values.flatten(),
           C=(dat.TRANS[tchoose,:,ind]).values.flatten(),cmap=trans_cmap,gridsize=50,vmin=-0.075,vmax=0.075,)
    axx.contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    axx.axvline(sref_all,color='k')
    return ts


def TS_transplot_monthly():
    for ii in range(21):
        f,(ax1,ax2,ax3)=subplots(1,3,figsize=(14,3.5),sharex=True,sharey=True)
        ts=TS_trans_onetime(eastind,ax1,ii)
        ax1.set_title('East')
        ts=TS_trans_onetime(westind,ax2,ii)
        ax2.set_title('West')
        ts=TS_trans_onetime(allind,ax3,ii)
        ax3.set_title('All')
        cbax=f.add_axes([0.95,0.1,0.02,0.8])
        colorbar(ts,cax=cbax,label='Transport [Sv]')
        ax1.set_xlim(31,36)
        ax1.set_ylim(-3,15)
        f.text(0.5, -0.05, 'Salinity', ha='center',fontsize=12)
        f.text(0.05, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=12)
        title(str(dat.TIME[ii].values)[:10])
        savefig(figdir+'Freshwater/TS/Xport_TS_'+str(ii).zfill(2)+'.png',bbox_inches='tight')

# TS_transplot_monthly()


def TS_fwt(ind,axx):
    axx.contour(salvec,tmpvec,pdenmat,levels=arange(24.06,28,0.4),colors='grey')
    ts=axx.hexbin(dat.PSAL[:,:,ind].values.flatten(),dat.PTMP[:,:,ind].values.flatten(),
           C=(dat.FWT[:,:,ind]).values.flatten(),cmap=trans_cmap,gridsize=50,vmin=-1,vmax=1,)
    axx.contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    axx.axvline(sref_all,color='k')
    return ts

def TS_fwtplot():
    f,(ax1,ax2,ax3)=subplots(1,3,figsize=(14,3.5),sharex=True,sharey=True)
    ts=TS_fwt(eastind,ax1)
    ax1.set_title('East')
    ts=TS_fwt(westind,ax2)
    ax2.set_title('West')
    ts=TS_fwt(allind,ax3)
    ax3.set_title('All')
    cbax=f.add_axes([0.95,0.1,0.02,0.8])
    colorbar(ts,cax=cbax,label='Freshwater transport [mSv]')
    ax1.set_xlim(31,36)
    ax1.set_ylim(-3,15)
    f.text(0.5, -0.05, 'Salinity', ha='center',fontsize=12)
    f.text(0.05, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=12)
    savefig(figdir+'Freshwater/FWT_TS_alltime.png',bbox_inches='tight')

TS_fwtplot()
