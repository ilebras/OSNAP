from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Summarize_1912/'

################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################

osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')

osnap_obs=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full_wMOC.nc')
fs_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_fs_xray_1912.nc')
bso_obs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_bso_xray_1912.nc')

WM_nor=xr.open_dataset(datadir+'NorESM/NorESM_WMs_1912.nc')
WM_obs=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')


################################################################################################################################
################################################################################################################################
########################################      PLOT MEAN SECTION PANELS    #######################################################
################################################################################################################################
################################################################################################################################
sigmax=osnap.SIGMEAN.values

fslim=6

vmax=0.3
univec['VELO']=['across track velocity',arange(-vmax,vmax+0.05,0.05),cm.RdBu_r,arange(-vmax,vmax+0.1,0.1),'[m/s]']
univec['PTMP']=['pot. temperature', arange(-1,11,0.5),cm.RdYlBu_r,range(0,11,2),'[$^\\circ$C]']
univec['PSAL']=['salinity',arange(34,35.6,0.1),cm.PiYG_r, array([34., 34.5 ,35,35.5,]),'']

prop_key=['VELO','PSAL','PTMP']
label_key=['velocity [m s$^{-1}$]','salinity','pot. temperature [$^\circ$C]']

# VELOCITY, SALINITY, TEMPERATURE AND DENSITY MEANS OVER PERIOD OF OVERLAP
def plot_each(dat,xxvar,axx,var,ymax,sigi):
    filled=axx.contourf(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),univec[var][1],cmap=univec[var][2],extend='both')
    axx.contour(xxvar,dat.DEPTH,dat[var].mean(dim='TIME'),levels=univec[var][1][::2],colors='k')
    # if 'PDEN' in dat:
    axx.contour(xxvar,dat.DEPTH,dat['PDEN'].mean(dim='TIME'),levels=[sigi],colors='k',linewidths=4) # add NorESM isopycnal of maximum overturning (for 2014-2016)
    axx.set_facecolor('k')
    axx.set_ylim(ymax,0)
    return filled


def comp_obs_model(nor,obs,ymax,tit,xlab,stit):
    f,axx=subplots(3,2,figsize=(12,10),sharex=True,sharey=True)
    cout={}
    for ii in range(3):
        if 'Barents' in tit:
            nor_x=nor.LATITUDE
            if ii ==0:
                obs_x=obs.MLAT
            else:
                obs_x=obs.LATITUDE
        else:
            nor_x=nor.LONGITUDE
            if 'Fram' in tit:
                if ii==0:
                    obs_x=obs.MLON
                else:
                    obs_x=obs.LONGITUDE
            else:
                obs_x=obs.LONGITUDE
        plot_each(obs,obs_x,axx[ii,0],prop_key[ii],ymax,sigmax_obs)
        cout[ii]=plot_each(nor,nor_x,axx[ii,1],prop_key[ii],ymax,sigmax)
    fs=14
    axx[0,0].set_title('Observations',fontsize=fs)
    axx[0,1].set_title('NorESM model',fontsize=fs)
    f.text(0.05, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fs)
    f.text(0.5, 0.05, xlab, ha='center',fontsize=fs)
    suptitle(tit,fontsize=fs+2)
    caxit={}
    xc=0.95
    cwi=0.01
    cle=0.2
    caxit[0]=f.add_axes([xc,0.67,cwi,cle])
    caxit[1]=f.add_axes([xc,0.4,cwi,cle])
    caxit[2]=f.add_axes([xc,0.13,cwi,cle])
    for ii in range(3):
        colorbar(cout[ii],label=label_key[ii],cax=caxit[ii])
    savefig(figdir+'CompSect_Obs_NorESM_'+stit+'.png',bbox_inches='tight')

comp_obs_model(fs,fs_obs,3000,'Fram Strait','Longitude [$^\circ$W]','fs')

comp_obs_model(bso,bso_obs,500,'Barents Sea Opening','Latitude [$^\circ$N]','bso')

comp_obs_model(osnap,osnap_obs,3200,'OSNAP East','Longitude [$^\circ$W]','osnap')


################################################################################################################################
################################################################################################################################
###########################################      PLOT TS SPACE    #######################################################
################################################################################################################################
################################################################################################################################

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



coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}
coldic.keys()

WM_lit=xr.Dataset({'PSAL':('WM',[35.2,35,33.8,33.8,35]),'PTMP':('WM',[0,0,0,0,4]),
                   'TRANS':('WM',[17.2,-12.7,-3.8,8,-8]),'COLOR':('WM',['red','grey','royalblue','purple','orange'])},
                   coords={'WM':['AWS', 'DWS', 'PWS', 'PWN', 'AWN'],})

def plot_TS_bylayer():
    f,[ax1,ax2]=subplots(1,2,figsize=(11,5),sharex=True,sharey=True)
    for wm in WM_nor.WM:
        ax2.scatter(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,
                s=WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
    for wm in WM_obs.WM:
        ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,
                s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,zorder=50,linewidth=3,label=str(wm.values),color=coldic[str(wm.values)])
    for wm in ['AWN','PWN']:
        ax1.scatter(WM_lit['PSAL'].sel(WM=wm).values,WM_lit['PTMP'].sel(WM=wm).values,
                s=WM_lit['TRANS'].sel(WM=wm).values**2*4,zorder=50,linewidth=3,label=wm,color=coldic[wm])
    lgnd=ax1.legend(loc=(2.3,0.2))
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [40]
    for axx in [ax1,ax2]:
        axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
        axx.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=5)
    f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    xlim(33.5,36)
    ylim(-2,12)
    f.suptitle('Transport-weighted water mass properties',fontsize=16)
    ax1.set_title('Observations',fontsize=14)
    ax2.set_title('NorESM',fontsize=14)
    savefig(figdir+'TS_comp.png',bbox_inches='tight')

plot_TS_bylayer()

Se={'AWS': 0.06, 'PWS': 0.8, 'DWS': 0.1, 'AWN': 0.2, 'PWN': 0.8, 'FW': 0, 'SI': 0}
Ue={'AWS': 2.0,'PWS': 2.0,'DWS': 2.0,'AWN': 2.0,'PWN': 2.0,'FW': 0.01,'SI': 0.01}

def barplot():
    f,[ax1,ax2,ax3]=subplots(3,1,figsize=(16,11))
    ax1.bar(WM_lit['WM'],WM_lit['TRANS'],color=WM_lit['COLOR'].values,yerr=[Ue[kk] for kk in WM_lit.WM.values],capsize=10)
    ax1.plot(WM_nor['WM'],WM_nor['TRANS'].mean(dim='TIME'),'d',color='k',markersize=20)
    ax2.bar(WM_lit['WM'],WM_lit['PSAL'],color=WM_lit['COLOR'].values,yerr=[Se[kk] for kk in WM_lit.WM.values],capsize=10)
    ax2.set_ylim(32.5,36)
    ax2.plot(WM_nor['WM'],WM_nor['PSAL'].mean(dim='TIME'),'d',color='k',markersize=20)
    ax3.bar(WM_lit['WM'],WM_lit['TRANS']*WM_lit['PSAL'],color=WM_lit['COLOR'].values,yerr=[Ue[kk]*WM_lit['PSAL'].sel(WM=kk)+Se[kk]*WM_lit['TRANS'].sel(WM=kk) for kk in WM_lit.WM.values],capsize=10)
    ax3.plot(WM_nor['WM'],WM_nor['PSAL'].mean(dim='TIME')*WM_nor['TRANS'].mean(dim='TIME'),'d',color='k',markersize=20)
    ax1.set_ylabel('Transport, U [Sv]')
    ax2.set_ylabel('Salinity, S')
    ax3.set_ylabel('Salinity budget term, U*S')
    savefig(figdir+'SU_comp.png',bbox_inches='tight')


barplot()
