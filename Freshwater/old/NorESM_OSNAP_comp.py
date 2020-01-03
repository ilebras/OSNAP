from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'

################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################


#### load osnap data and cut out the eastern portion
lonbnd=-44
dat=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.nc')
osnap=dat.sel(LONGITUDE=slice(lonbnd,0))
osnap['LATITUDE']=dat['LATITUDE'].values[dat.LONGITUDE>lonbnd]

#### load noresm and force formatting to be the same
def load_noresm():
    noresm=xr.open_dataset(glob.glob(datadir+'NorESM/*OSNAP*.nc')[0])
    noresm=noresm.rename({'time': 'TIME','depth':'DEPTH','cell':'LONGITUDE'})
    noresm=noresm.rename({'vvel': 'VELO','salt':'PSAL','temp':'PTMP'})
    noresm=noresm.assign_coords(LONGITUDE=(noresm.lon.values))
    noresm=noresm.assign_coords(LATITUDE=(noresm.lat.values))
    noresm['vertices_lon']=noresm['vertices_lon']-360
    noresm=noresm.drop('lat').drop('lon')
    startyear=2010
    startmonth=1
    endyear=2018
    endmonth=12
    noresm['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    noresm=noresm.sel(TIME=slice(osnap.TIME[0],osnap.TIME[-1]))
    return noresm



noresm=load_noresm()


def add_PDEN(xray):
    if 'PRES' in list(xray.data_vars):
            PRES_out=xray['PRES']
    else:
            PRES_out=gsw.p_from_z(-xray['DEPTH'],60)
    SA_out=NaN*xray['PSAL'].copy()
    for ii,pp in enumerate(PRES_out):
        for jj,ll in enumerate(xray.LONGITUDE):
            SA_out[:,ii,jj]=gsw.SA_from_SP(xray['PSAL'][:,ii,jj],pp,ll,xray.LATITUDE[jj])
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['PDEN']=(('TIME','DEPTH','LONGITUDE'),PD_out)

add_PDEN(noresm)

################################################################################################################################
################################################################################################################################
###########################################      PLOT MEAN SECTIONS    #######################################################
################################################################################################################################
################################################################################################################################

# VELOCITY, SALINITY, TEMPERATURE AND DENSITY MEANS OVER PERIOD OF OVERLAP
def plot_each(dat,axx,var,var2,eachtit):
        filled=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat[var].mean(dim='TIME'),univec[var2][1],cmap=univec[var2][2],extend='both')
        axx.contour(dat.LONGITUDE,dat.DEPTH,dat[var].mean(dim='TIME'),levels=univec[var2][1][::2],colors='k')
        axx.set_xlabel('Longitude [$^\circ$W]')
        axx.set_facecolor('k')
        axx.set_ylim(3750,0)
        axx.set_xlim(-44,-8)
        axx.set_title(eachtit,fontsize=16)
        return filled

def comp_secs(var,var2,varlab):
    f,ax=subplots(1,2,figsize=(16,6),sharey=True)
    subplots_adjust(wspace=0.05)
    plot_each(osnap,ax[0],var,var2,'OSNAP Observations')
    fill=plot_each(noresm,ax[1],var,var2,'NorESM')
    ax[0].set_ylabel('depth [m]')
    cbax=f.add_axes([0.925,0.15,0.015,0.7])
    colorbar(fill,cax=cbax,label=varlab,ticks=univec[var2][3])
    savefig(figdir+'CompMeanSec_'+var2+'.png',bbox_inches='tight')

comp_secs('PDEN','pden','potential density [kg m$^{-3}$]')

vmax=0.3
univec['VELO']=['across track velocity',arange(-vmax,vmax+0.01,0.025),cm.RdBu_r,arange(-vmax,vmax+0.05,0.1),'[m/s]']
comp_secs('VELO','VELO','velocity [m s$^{-1}$]')

univec['sal']=['salinity',arange(34,35.6,0.1),cm.PiYG_r, array([34., 34.5 ,35,35.5,]),'']


comp_secs('PSAL','sal','salinity')

univec['tmp']=['pot. temperature', arange(-1,11,0.5),cm.RdYlBu_r,range(0,11,2),'[$^\\circ$C]']

comp_secs('PTMP','tmp','pot. temperature [$^\circ$C]')


################################################################################################################################
################################################################################################################################
###########################################      PLOT OVERALL TS SPACE    #######################################################
################################################################################################################################
################################################################################################################################

def TS_hex_plot(dat,titi):
    figure()
    hexbin(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),cmap='hot_r',mincnt=1,bins='log')
    colorbar(label='# of Observations')
    xlabel('salinity')
    ylabel('pot.temperature')
    xlim(34,36)
    ylim(-1,14)
    title(titi)
    savefig(figdir+'TShex_'+titi+'.png',bbox_inches='tight')

TS_hex_plot(osnap,'OSNAP')
TS_hex_plot(noresm,'NorESM')


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


def TS_comp_plot():
    figure()
    dat=osnap
    plot(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),'o',label='OSNAP')
    dat=noresm
    plot(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),'o',label='NorESM')

    xlabel('salinity')
    ylabel('pot.temperature')
    xlim(34,36)
    ylim(-1,14)
    contour(salvec,tmpvec,pdenmat,colors='k',levels=arange(25,29,0.2),zorder=5)
    legend(fontsize=16)
    savefig(figdir+'TScomp_all.png',bbox_inches='tight')


TS_comp_plot()

################################################################################################################################
################################################################################################################################
################   INVESTIGATE THE  STREAMFUNCTIONS AND WATER MASS PARTITIONING   #############################################
################################################################################################################################
################################################################################################################################


# depth_bounds:

nor_dist=array([57382.5094665826, 57367.6635038391, 57353.5409692236, 57340.1361679362, 57327.4434986111, 57315.4573352359, 57304.1721310823, 57293.5824492824,
57957.9063952388, 57948.4302559112, 57939.6484924917, 58619.6339361550, 58612.0555633080, 58605.1744353929, 59301.2369172294, 59295.6417676100, 59290.7490131978,
60003.1944856841, 59999.6566352716, 59996.8270689983, 59994.7052065002, 59993.2908518308, 59992.5837196858, 59992.5837196858, 59993.2908518308,59994.7052065002,
59996.8270689983, 59999.6566352716, 60003.1944856841, 60007.4409760742, 60012.3967978497, 60018.0626167572, 60024.4390999926, 60031.5270237768, 60039.3271940887, 60047.8404139799])

nor_depthdiff=diff([0, 2.50000000000000, 7.50000000000000, 12.5000000000000, 17.5000000000000, 22.5000000000000, 27.5000000000000, 35, 45, 56.2000000000000, 68.8000000000000, 81.2000000000000, 93.8000000000000, 106.200000000000, 118.800000000000, 131.200000000000, 143.800000000000, 162.500000000000, 187.500000000000, 212.500000000000, 237.500000000000, 262.500000000000, 287.500000000000, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1562.50000000000, 1687.50000000000, 1812.50000000000, 1937.50000000000, 2125, 2375, 2625, 2875, 3125, 3375, 3625, 3875, 4125, 4375, 4625, 4875, 5125, 5375, 5625, 5875, 6125, 6375, 6625, 8000])
noresm['AREA']=(['DEPTH','LONGITUDE'],tile(nor_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(nor_dist),1]).T)


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
noresm=getpsi(noresm)

osnap.SIGMEAN

noresm.SIGMEAN


def comp_MOC():
    figure(figsize=(7,3))
    osnap.MOCMEAN.plot(color='k',linestyle='--',label='time-varying $\sigma_{max}$')
    osnap.MOCMEAN.plot(color='k',linestyle='-',label='constant (mean) $\sigma_{max}$')
    osnap.MOC.plot(color='C0',linestyle='--',label='',linewidth=2)
    noresm.MOC.plot(color='C1',linestyle='--',label='',linewidth=2)
    osnap.MOCMEAN.plot(label='NorESM',color='C0',linewidth=2)
    noresm.MOCMEAN.plot(label='OSNAP',color='C1',linewidth=2)
    legend(loc=(1.05,0),fontsize=12)
    xlabel('')
    ylabel('Oveturning transport [Sv]')
    savefig(figdir+'MOCcomp.png',bbox_inches='tight')

comp_MOC()

def plot_comp_psi():
    figure()
    plot(osnap.PSI,osnap.DENLEV,'C0',alpha=0.2,label='')
    plot(noresm.PSI,noresm.DENLEV,'C1',alpha=0.2,label='')
    plot(osnap.PSI.mean(dim='TIME'),osnap.DENLEV,linewidth=3,label='OSNAP')
    plot(noresm.PSI.mean(dim='TIME'),noresm.DENLEV,linewidth=3,label='NorESM')
    xlabel('Streamfunction [Sv]')
    ylabel('pot. density anomaly [kg m$^{-3}$]')
    legend(fontsize=16)
    ylim(28.3,26)
    savefig(figdir+'PsiComp.png',bbox_inches='tight')

plot_comp_psi()

############################################################################################
################   WATER MASS PARTITIONING   #############################################

AW={}
PW={}
DW={}
for prod in ['osnap','noresm']:
    AW[prod]={}
    PW[prod]={}
    DW[prod]={}
    if prod=='osnap':
        xray=osnap
    elif prod=='noresm':
        xray=noresm
    for type in ['mean','tvar']:
        AW[prod][type]={}
        PW[prod][type]={}
        DW[prod][type]={}
        if type=='mean':
            sigmax=xray.SIGMEAN
        elif type=='tvar':
            sigmax=xray.SIGMAX
        for var in ['PSAL','PTMP','PDEN','TRANS']:
            if var=='TRANS':
                DW[prod][type][var]=xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
                PW[prod][type][var]=xray['TRANS'].where(xray.LONGITUDE<-41).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
                AW[prod][type][var]=xray['TRANS'].where(xray.LONGITUDE>=-41).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')

            else:
                DW[prod][type][var]=(xray[var]*xray['TRANS']).where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
                PW[prod][type][var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE<-41).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE<-41).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')
                AW[prod][type][var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=-41).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=-41).where(xray.PDEN<sigmax).sum(dim='DEPTH').sum(dim='LONGITUDE')


for var in ['PSAL','PTMP','PDEN','TRANS']:
    f,axx=subplots(1,2,figsize=(12,3),sharex=True,sharey=True)
    for ii,prod in enumerate(['osnap','noresm']):
        DW[prod]['mean'][var].plot(linestyle='-',color='k',label='DW',ax=axx[ii])
        AW[prod]['mean'][var].plot(linestyle='-',color='r',label='AW',ax=axx[ii])
        PW[prod]['mean'][var].plot(linestyle='-',color='b',label='PW',ax=axx[ii])
        DW[prod]['tvar'][var].plot(linestyle='--',color='k',label='',ax=axx[ii])
        AW[prod]['tvar'][var].plot(linestyle='--',color='r',label='',ax=axx[ii])
        PW[prod]['tvar'][var].plot(linestyle='--',color='b',label='',ax=axx[ii])
    axx[0].set_ylabel(var)
    axx[0].set_xlabel('')
    axx[0].set_xlabel('')
    axx[1].set_ylabel('')
    axx[0].set_title('OSNAP',fontsize=16)
    axx[1].set_title('NorESM',fontsize=16)
    axx[1].legend(loc=(1.05,0.2))
    savefig(figdir+'WMComp_tvar_'+var+'.png',bbox_inches='tight')

wmvec=[AW,PW,DW]

def plot_TS_bylayer():
    figure(figsize=(5,4))
    for ii,prod in enumerate(['osnap','noresm']):
        scatter(AW[prod]['mean']['PSAL'].mean(dim='TIME'),AW[prod]['mean']['PTMP'].mean(dim='TIME'),s=[abs(AW[prod]['mean']['TRANS'].mean(dim='TIME'))**2*4 for wm in wmvec],c=['r'],zorder=50,linewidth=3,alpha=1-0.5*ii,label='AW '+prod)
        scatter(PW[prod]['mean']['PSAL'].mean(dim='TIME'),PW[prod]['mean']['PTMP'].mean(dim='TIME'),s=[abs(PW[prod]['mean']['TRANS'].mean(dim='TIME'))**2*4 for wm in wmvec],c=['b'],zorder=50,linewidth=3,alpha=1-0.5*ii,label='PW '+prod)
        scatter(DW[prod]['mean']['PSAL'].mean(dim='TIME'),DW[prod]['mean']['PTMP'].mean(dim='TIME'),s=[abs(DW[prod]['mean']['TRANS'].mean(dim='TIME'))**2*4 for wm in wmvec],c=['k'],zorder=50,linewidth=3,alpha=1-0.5*ii,label='DW '+prod)
    contour(salvec,tmpvec,pdenmat,colors='k',levels=[osnap.SIGMEAN.values],zorder=5,alpha=1)
    contour(salvec,tmpvec,pdenmat,colors='k',levels=[noresm.SIGMEAN.values],zorder=5,alpha=0.5)
    xlabel('salinity')
    ylabel('pot.temperature')
    xlim(34,36)
    ylim(0,12)
    lgnd=legend(loc=(1.05,0.2),ncol=2)
    for ii in range(6):
        lgnd.legendHandles[ii]._sizes = [40]
    title('Transport-weighted water mass properties')
    savefig(figdir+'TS_WMcomp.png',bbox_inches='tight')

plot_TS_bylayer()

prod='osnap'
val='PSAL'
osnap.SIGMEAN
AW[prod]['mean'][val].min(dim='TIME')
AW[prod]['mean'][val].max(dim='TIME')
PW[prod]['mean'][val].min(dim='TIME')
PW[prod]['mean'][val].max(dim='TIME')
DW[prod]['mean'][val].min(dim='TIME')
DW[prod]['mean'][val].max(dim='TIME')
