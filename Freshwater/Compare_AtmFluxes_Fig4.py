from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/airsea/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'

era5=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/ERA5_ilebras/ERA5_from2000_ephf_cutint_ilebras.nc')
ncep1=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/NCEP_CFS/NCEP_pre2011_ephf_cutint_ilebras.nc')
ncep2=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/NCEP_CFSv2/NCEP_ephf_cutint_ilebras.nc')
ncep=xr.concat([ncep1,ncep2],dim='time')
noresm=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_18yrs_2004.nc')
hf1=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_200001-200912.nc')
hf2=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_201001-201812.nc')
noresm_hf=xr.concat([hf1,hf2],dim='time')
jra55=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/JRA55/JRA55_2000_2013_ephf_cutint_ilebras.nc')
noresm['ep_int']=-noresm.evap-noresm.liqprec-noresm.solprec
noresm['hf_int']=('TIME',-noresm_hf.NORDIC_hflx.values)
noresm=noresm.rename({'TIME':'time'})
era5['ep_int']=-era5.e_int-era5.tp_int

era5['hf_int']=-era5.sshf_int-era5.slhf_int-era5.ssr_int-era5.str_int
noresm['mltfrz']=-1*abs(noresm['mltfrz'])
noresm['runoff']=-1*abs(noresm['runoff'])


def complot(var,ylab):
    labvec=['ERA 5','NCEP CFS','NCEP CFSv2','JRA55 (NorESM)','JRA55']
    figure(figsize=(15,3))
    for ii,xray in enumerate([era5,ncep1,ncep2,noresm,jra55]):
        (-xray[var]).plot(label=labvec[ii],linewidth=2)
    xlim(datetime.datetime(2000,1,1),datetime.datetime(2019,1,1))
    legend(loc=(1.01,0.3))
    ylabel(ylab)
    savefig(figdir+'ReanalysisComp_'+var+'.png',bbox_inches='tight')


complot('ep_int','Precipitation-Evaporation [Sv]')

complot('hf_int','Air-Sea Heat Flux [TW]')

coldic={'ERA5':'limegreen','JRA55':'#8c510a','NCEP':'#f768a1','JRA55 (NorESM)':'#8c510a'}


def all_airsea_fig():
    f,axx=subplots(2,2,figsize=(10,6),sharex=True)
    ylp={'ep_int':-0.03,'mltfrz':0.02,'runoff':0.06,'hf_int':-600}
    for var in ['ep_int','mltfrz','hf_int','runoff']:
        if var=='ep_int':
            prodvec=[noresm,era5,ncep]
            labvec=['JRA55 (NorESM)','ERA5','NCEP']
            jj=0
            kk=0
        elif var=='hf_int':
            prodvec=[era5,ncep,jra55,noresm]
            labvec=['ERA5','NCEP','JRA55','JRA55 (NorESM)']
            jj=1
            kk=1
        else:
            prodvec=[noresm]
            labvec=['JRA55 (NorESM)']
        if var=='mltfrz':
            jj=0
            kk=1
        elif var=='runoff':
            jj=1
            kk=0
        for ii,xray in enumerate(prodvec):
            if 'JRA55'==labvec[ii]:
                lsty='--'
            else:
                lsty='-'
            seas_mean=(-xray[var]).sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').mean(dim='time')
            seas_std=(-xray[var]).sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').std(dim='time')
            seas_mean.plot(label=labvec[ii],linewidth=3,color=coldic[labvec[ii]],linestyle=lsty,ax=axx[jj,kk],zorder=200)
            axx[jj,kk].fill_between(range(1,13),seas_mean-seas_std,seas_mean+seas_std,alpha=0.3,color=coldic[labvec[ii]],zorder=100)
            if 'hf_int' not in var:
                axx[jj,kk].plot(6.5,mean(seas_mean),color=coldic[labvec[ii]],marker='d',zorder=3,markersize=10)
            else:
                axx[jj,kk].plot((ii+11.5)/2,mean(seas_mean),color=coldic[labvec[ii]],marker='d',zorder=3,markersize=10)
            if labvec[ii]=='JRA55 (NorESM)':
                if var=='hf_int':
                    axx[jj,kk].text(5.5,ylp[var],str(int(mean(seas_mean).values))+' TW',color=coldic[labvec[ii]],fontsize=13,fontweight='bold')
                else:
                    axx[jj,kk].text(5.5,ylp[var],str(int(mean(seas_mean*1e3).values))+'mSv',color=coldic[labvec[ii]],fontsize=13,fontweight='bold')
    axx[0,0].set_xlabel('')
    axx[0,0].set_ylim(-0.05,0.25)
    axx[0,1].set_ylim(-0.05,0.25)
    axx[1,0].set_ylim(-0.05,0.25)
    axx[0,0].set_yticks(arange(0,0.3,0.1))
    axx[0,1].set_yticks(arange(0,0.3,0.1))
    axx[1,0].set_yticks(arange(0,0.3,0.1))
    axx[1,1].set_yticks(arange(600,-900,-300))
    fsi=13
    axx[0,0].set_title('Precipitation-Evaporation [Sv]',fontsize=fsi)
    axx[0,1].set_title('Sea ice [Sv]',fontsize=fsi)
    axx[1,0].set_title('Runoff [Sv]',fontsize=fsi)
    axx[1,1].set_title('Heat flux [TW]',fontsize=fsi)
    for qq in range(2):
        for rr in range(2):
            axx[qq,rr].set_xlabel('')
            axx[qq,rr].set_ylabel('')
    f.text(0.5, 0.05, 'Month', ha='center',fontsize=fsi)
    axx[0,0].legend(loc=(0.3,-1.75),fontsize=fsi,ncol=3)
    xlim(1,12)
    savefig(figdir_paper+'AirSeaF_FWsrc.png',bbox_inches='tight')
    savefig(figdir_paper+'AirSeaF_FWsrc.pdf',bbox_inches='tight')


all_airsea_fig()
(noresm.mltfrz+noresm.ep_int).groupby('time.month').mean(dim='time').plot()

era5['hf_int'].sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').mean(dim='time').mean(dim='month')
jra55['hf_int'].sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').mean(dim='time').mean(dim='month')
noresm['hf_int'].sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').mean(dim='time').mean(dim='month')
noresm['mltfrz'].sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').mean(dim='time').mean(dim='month')*1000

noresm['ep_int'].sel(time=slice('2000-1-1','2019-1-1')).groupby('time.month').mean(dim='time').mean(dim='month')*1000
