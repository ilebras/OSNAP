from pylab import *
import xarray as xr
from aux_funcs import *

TS=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Gridded_TS_201408_201604_2018.nc')
V=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Gridded_V_201408_201604_2018.nc')

dat=xr.merge([TS,V])

# xport=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Transports_201408_201604_2018.nc')
# psi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Streamfunction_201408_201604_2018.nc')

# bathy=io.loadmat(datadir+'Shipboard/jr302_osnap_sim_bathfromPenny.mat')
# bathlon=bathy['osnapnogap'][0][0][1][0]
# bathdepth=bathy['osnapnogap'][0][0][2][0]

bathy=io.loadmat(datadir+'OSNAPbathy_Fli.mat')

eastind=dat.LONGITUDE>-45
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
plot(bathy['lon'].flatten(),bathy['bathy'].flatten())


plot(dat.LONGITUDE,dat['FWT_east'].sum(dim='DEPTH').mean(dim='TIME').cumsum(dim='LONGITUDE'))

def plotsal():
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
    axx[1].fill_between(bathy['lon'].flatten(),-bathy['bathy'].flatten(),[4000]*len(bathy['lon'].flatten()),color='k')
    axx[1].text(-44,3200,'Greenland',color='white',fontsize=16)
    axx[1].text(-34.5,3200,'Mid-Atlantic Ridge',color='white',fontsize=16)
    axx[1].text(-11,3200,'Scotland',color='white',fontsize=16)
    # axx[1].text(-19,2800,'Mean 2014-2016',color='white',fontsize=15)

    cax2=f.add_axes([cxx,0.125,0.02,cy])
    sal=axx[2].contourf(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),slevs,cmap=sal_cmap,extend='both')
    c1=axx[2].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[34.9],colors='grey')
    axx[2].set_ylabel('depth [m]')
    axx[2].fill_between(bathy['lon'].flatten(),-bathy['bathy'].flatten(),[4000]*len(bathy['lon'].flatten()),color='k')
    clabel(c1,fmt='%1.1f')
    c2=axx[2].contour(dat.LONGITUDE,dat.DEPTH,dat.PSAL.mean(dim='TIME'),levels=[sref_east],colors='k')
    clabel(c2,fmt='%1.2f')
    axx[2].set_ylim(3500,0)
    axx[2].set_xlabel('Longitude')
    colorbar(sal,ax=axx[2],ticks=slevs[::2],cax=cax2,label='Salinity')
    axx[0].set_xlim(-44,-7)
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/OSNAP_FWsections.pdf',bbox_inches='tight')

plotsal()

trans_cols = [(37,52,148),(67,162,202) ,(255,255,255),(254,178,76),(240,59,32)]
trans_cmap = make_cmap(trans_cols,position=[0,0.45,0.5,0.55,1],bit=True)

def TS_transmn():
    contour(salvec,tmpvec,pdenmat,levels=arange(27.06,28,0.2),colors='grey',zorder=1)
    ts=hexbin(dat.PSAL[:,:,eastind].mean(dim='TIME').values.flatten(),dat.TEMP[:,:,eastind].mean(dim='TIME').values.flatten(),
           C=(dat.VELO[:,:,eastind].mean(dim='TIME')*dat.AREA).values.flatten()/1e6,cmap=trans_cmap,vmin=-0.05,vmax=0.05,zorder=2)
    contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    colorbar(ts)
    axvline(sref_east,color='k')
    xlim(33.5,35.5)
    ylim(0,12)
TS_transmn()

salvec=linspace(33,36,100)
tmpvec=linspace(-2,16,100)
salmat,tmpmat=meshgrid(salvec,tmpvec)

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])

SAdat=gsw.SA_from_SP(dat.PSAL,0,-45,57)

dat['PTMP']=dat.TEMP.copy()
for ii in range(len(dat.TIME)):
    for jj in range(len(dat.LONGITUDE)):
            dat['PTMP'][ii,:,jj]=gsw.pt0_from_t(SAdat[ii,:,jj],dat.TEMP[ii,:,jj],gsw.p_from_z(-dat.DEPTH,57))



def TS_trans(ind,axx):
    axx.contour(salvec,tmpvec,pdenmat,levels=arange(26.06,28,0.2),colors='grey')
    ts=axx.hexbin(dat.PSAL[:,:,ind].values.flatten(),dat.PTMP[:,:,ind].values.flatten(),
           C=(dat.TRANS[:,:,ind]).values.flatten(),cmap=trans_cmap,gridsize=50,vmin=-0.075,vmax=0.075,)
    axx.contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    axx.axvline(sref_east,color='k')
    return ts


def TS_transplot():
    f,(ax1,ax2)=subplots(1,2,figsize=(9,3),sharex=True,sharey=True)
    ts=TS_trans(eastind,ax1)
    ax1.set_title('OSNAP East')
    ts=TS_trans(CFind,ax2)
    ax2.set_title('East Greenland array')
    cbax=f.add_axes([0.95,0.1,0.02,0.8])
    colorbar(ts,cax=cbax,label='Transport [Sv]')
    ax1.set_xlim(33,36)
    ax1.set_ylim(-1,15)
    f.text(0.5, -0.05, 'Salinity', ha='center',fontsize=12)
    f.text(0.05, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=12)
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/OSNAP_TS_trans.pdf',bbox_inches='tight')


TS_transplot()
def TS_transplot1():
    f,ax1=subplots(1,1,figsize=(4.5,3),sharex=True,sharey=True)
    ts=TS_trans(eastind,ax1)
    colorbar(ts,label='Transport [Sv]')
    # ax1.set_title('OSNAP East')
    xlabel('Salinity')
    ylabel('Potential temperature [$^\circ$C]')
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/OSNAP_TS_trans1.pdf',bbox_inches='tight')
    ax1.plot([34.4,34.4],[2.5,6.1],'r-',linewidth=3)
    ax1.plot([35.4,35.4],[9.1,11],'r-',linewidth=3)
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/OSNAP_TS_trans1_wline.pdf',bbox_inches='tight')

TS_transplot1()

def TS_fwt(ind,axx):
    axx.contour(salvec,tmpvec,pdenmat,levels=arange(26.06,28,0.2),colors='grey')
    ts=axx.hexbin(dat.PSAL[:,:,ind].values.flatten(),dat.PTMP[:,:,ind].values.flatten(),
           C=(dat.FWT_east[:,:,ind]).values.flatten(),cmap=trans_cmap,gridsize=50,vmin=-1,vmax=1,)
    axx.contour(salvec,tmpvec,pdenmat,levels=[27.66],colors='k')
    axx.axvline(sref_east,color='k')
    return ts

def TS_fwtplot():
    f,(ax1,ax2)=subplots(1,2,figsize=(9,3),sharex=True,sharey=True)
    ts=TS_fwt(eastind,ax1)
    ax1.set_title('OSNAP East')
    ts=TS_fwt(CFind,ax2)
    ax2.set_title('East Greenland Current system')
    cbax=f.add_axes([0.95,0.1,0.02,0.8])
    colorbar(ts,cax=cbax,label='Freshwater transport [mSv]')
    ax1.set_xlim(33,36)
    ax1.set_ylim(-1,15)
    f.text(0.5, -0.05, 'Salinity', ha='center',fontsize=12)
    f.text(0.05, 0.5, 'Potential temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=12)
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/OSNAP_TS_fwt.pdf',bbox_inches='tight')
TS_fwtplot()

xxx
#
#

# CTdat=gsw.CT_from_pt(SAdat,dat.PTMP)
#
# for ii in range(len(dat.TIME)):
#     for jj in range(len(dat.LONGITUDE)):
#             dat.PDEN[ii,:,jj]=gsw.sigma0(SAdat[ii,:,jj],CTdat[ii,:,jj])

#
# denvec=arange(24.7,28.1,0.1)
# middenvec=denvec[:-1]+diff(denvec)/2
# varvec=['TRANS','PSAL','FWT','FWT_east']
# denmats={}
# varvec=['PSAL']
# for var in varvec:
#     denmats[var]=zeros((len(dat.TIME),len(middenvec),len(dat.LONGITUDE)))
#     for ll in range(len(dat.LONGITUDE)):
#         for dd in range(len(middenvec)):
#             if var=='PSAL':
#                 denmats[var][:,dd,ll]=dat[var][:,:,ll].where(dat.PDEN[:,:,ll]>denvec[dd]).where(dat.PDEN[:,:,ll]<=denvec[dd+1]).mean(dim='DEPTH');
#             else:
#                 denmats[var][:,dd,ll]=dat[var][:,:,ll].where(dat.PDEN[:,:,ll]>denvec[dd]).where(dat.PDEN[:,:,ll]<=denvec[dd+1]).sum(dim='DEPTH');
#
# def makedendat():
#     dendat=xr.Dataset({'PSAL': (['TIME','PDEN','LONGITUDE'],  denmats['PSAL']),
#                        'FWT': (['TIME','PDEN','LONGITUDE'],  denmats['FWT']/1e6),
#                        'FWT_east': (['TIME','PDEN','LONGITUDE'],  denmats['FWT_east']/1e6),
#                        'TRANS': (['TIME','PDEN','LONGITUDE'],  denmats['TRANS']/1e6)},
#                        coords={'PDEN': middenvec,
#                                'LONGITUDE': dat.LONGITUDE.values,
#                                'TIME': dat.TIME.values})
#
#     return dendat
#
# dendat=makedendat()

# dendat.PSAL.mean(dim='TIME').plot()
#
# dendat.FWT.mean(dim='TIME').plot()
#
# plot(dendat.FWT.mean(dim='TIME').sum(dim='LONGITUDE'),dendat.PDEN);
# axhline(27.35)
#
# figure(figsize=(12,3))
# dendat.TRANS.sum(dim='LONGITUDE').cumsum(dim='PDEN').max(dim='PDEN').plot()
#
# dvec=[NaN]*len(dat.LONGITUDE)
# for ii in range(len(dvec)):
#     if sum(~isnan(dat.AREA[:,ii]))>1:
#         dvec[ii]=max(dat.DEPTH[~isnan(dat.AREA[:,ii])])
#
# plot(dvec)
#
# dx=dat['AREA'].sum(dim='DEPTH')/dvec
#
# dendat.FWT.sum(dim='LONGITUDE').plot()
#
# (dendat['FWT']*dx).mean(dim='LONGITUDE').sum(dim='PDEN')
# dendat
#
# dendat['PSAL'].mean(dim='LONGITUDE').plot()
#
# (dendat['TRANS'].mean(dim='LONGITUDE')*dendat['PSAL'].mean(dim='LONGITUDE')).sum(dim='PDEN').plot()
# def plotFWTdecomp(fvar):
#     figure(figsize=(12,3))
#     (dendat['TRANS'].mean(dim='LONGITUDE')*(dendat['PSAL'].mean(dim='LONGITUDE'))/sref_all).sum(dim='PDEN').plot(label='Overturning')
#     # (dat[fvar]/1e6).sum(dim='LONGITUDE').sum(dim='DEPTH').plot(label='Total')
#     legend()
#
# plotFWTdecomp('FWT')
#
# plotFWTdecomp('FWT_east')
#
# dat
#
#
#
# contour(salvec,tmpvec,pdenmat,levels=arange(26.06,28,0.2),colors='grey')
#
# plot(dendat['FWT'].sum('LONGITUDE').mean(dim='TIME'),dendat.PDEN)
# [axhline(xx,color='grey') for xx in arange(26.06,28,0.2)]
# axhline(27.66,color='k')
# axvline(0,color='k')
# ylim(28,24)
