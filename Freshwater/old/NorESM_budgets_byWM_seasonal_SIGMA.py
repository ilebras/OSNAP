from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'
################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################
osnap={}
fs={}
bso={}
ns={}
osnap['dpth']=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_18yrs_2004.nc')
osnap['sig']=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_18yrs_2004_sigma.nc')
fs['sig']=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_18yrs_2004_sigma.nc')
fs['dpth']=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_18yrs_2004.nc')
bso['sig']=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_18yrs_2004_sigma.nc')
bso['dpth']=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_18yrs_2004.nc')
ns['sig']=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_18yrs_2004_sigma.nc')
ns['dpth']=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_18yrs_2004.nc')
so=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_18yrs_2004.nc')

vts1=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_volumetransports_200001-200912.nc')
vts2=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_volumetransports_201001-201812.nc')
vts=xr.concat([vts1,vts2],dim='time')
startyear=2000
startmonth=1
endyear=2018
endmonth=12
vts['time']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])

def comp_trans(var,var2):
    figure(figsize=(20,4))
    var['dpth'].TRANS.sum(dim='LONGITUDE').sum(dim='DEPTH').plot(label='depth')
    var['sig'].TRANS.sum(dim='LONGITUDE').sum(dim='DZ').plot(label='sigma')
    vts[var2].plot(label='diagnostic')
    ylabel('Transport [Sv]')
    title(var2[7:])
    legend()


comp_trans(osnap,'net_vt_OSNAP')
comp_trans(fs,'net_vt_FS')
comp_trans(bso,'net_vt_BSO')
comp_trans(ns,'net_vt_NS')

hf1=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_200001-200912.nc')
hf2=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_201001-201812.nc')
hf=xr.concat([hf1,hf2],dim='time')
hf['time']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])

WM=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc')


WM_obs=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')

vts['tot']=(-vts['net_vt_FS']-vts['net_vt_BSO']+vts['net_vt_OSNAP'])+vts['net_vt_NS']

vts.tot.mean()


coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}
################################################################################################################################
################################################################################################################################
###########################################      VOL BALANCE     #######################################################
################################################################################################################################
################################################################################################################################

def plot_volfluxes_seas():
    f,axx=subplots(2,3,figsize=(14,6),sharex=True)
    for gg,wm in enumerate(WM.WM.values):
        if gg>2:
            jj=1
            gg=gg-3
        else:
            jj=0
        trans_std=WM['TRANS'].sel(WM=wm).groupby('TIME.month').std(dim='TIME')
        trans_mean=WM['TRANS'].sel(WM=wm).groupby('TIME.month').mean(dim='TIME')
        axx[jj,gg].fill_between(range(1,13),trans_mean-trans_std,trans_mean+trans_std,color='grey',alpha=0.4)
        trans_mean.plot(color='k',linewidth=4,ax=axx[jj,gg],linestyle='--')
        WM_obs['TRANS'].sel(WM=wm).groupby('TIME.month').mean(dim='TIME').plot(color=coldic[wm],linewidth=4,ax=axx[jj,gg])
        axx[jj,gg].set_title(wm)

    axx[0,0].set_ylim(8,22)
    axx[0,2].set_ylim(-22,-8)
    axx[1,0].set_ylim(-14,0)
    axx[1,1].set_ylim(0,14)


    vts_std=(vts['tot']).groupby('time.month').std(dim='time')
    vts_mean=(vts['tot']).groupby('time.month').mean(dim='time')
    axx[1,2].fill_between(range(1,13),vts_mean-vts_std,vts_mean+vts_std,color='grey',alpha=0.4)
    vts_mean.plot(color='k',linewidth=4,ax=axx[1,2],linestyle='--')
    # axhline(mean(vts_mean),color='green',linewidth=3)
    so_mean=so['FW+SI'].groupby('TIME.month').mean(dim='TIME')
    so_std=so['FW+SI'].groupby('TIME.month').std(dim='TIME')
    axx[1,2].fill_between(range(1,13),so_mean-so_std,so_mean+so_std,color='limegreen',alpha=0.4)
    so_mean.plot(color='limegreen',linewidth=3)

    for ggg in axx:
        for axi in ggg:
            axi.set_ylabel('')
            axi.set_xlabel('')
            axi.set_xlim(1,12)
    axx[1,2].set_title('Lateral volume transport divergence')
    fsz=14
    f.text(0.075, 0.5, 'Transport [Sv]', va='center',fontsize=fsz,rotation='vertical')
    f.text(0.5, 0.05, 'Month', ha='center',fontsize=fsz)
    savefig(figdir+'Budget_volume_seas.png',bbox_inches='tight')
    savefig(figdir+'Budget_volume_seas.pdf',bbox_inches='tight')



plot_volfluxes_seas()

################################################################################################################################
################################################################################################################################
###########################################    HEAT AND SALT STORAGE     #######################################################
################################################################################################################################
################################################################################################################################
cp=3850
rhow=1000
tera=10**12

# #HEAT YEARLY Note, yearly storage doesn't work out because I don't have the correctly gridded heat/salt fluxes
# htd_mean=(-cp*rhow*1e6*WM['TRANS']*WM['PTMP']/tera).sum(dim='WM').groupby('TIME.year').mean(dim='TIME')
# htd_mean.plot(label='-Lateral heat transport divergence',ax=axx[1,1],color='k',linewidth=2)
# hf_mean=hf['NORDIC_hflx'].groupby('time.year').mean(dim='time')
# hf_mean.plot(label='Air-Sea heat flux',color='limegreen',linewidth=2,ax=axx[1,1])

def plot_saltheat_storage():
    f,axx=subplots(1,2,figsize=(13,4),sharex='row')
    fsz=14
    #SALT MONTHLY
    std_mean=(WM['TRANS']*WM['PSAL']).sum(dim='WM').groupby('TIME.month').mean(dim='TIME')
    std_std=(WM['TRANS']*WM['PSAL']).sum(dim='WM').groupby('TIME.month').std(dim='TIME')
    axx[0].fill_between(range(1,13),std_mean-std_std,std_mean+std_std,color='grey',alpha=0.4)
    std_mean.plot(label='-Lateral salt transport divergence',ax=axx[0],color='k',linewidth=4)
    axx[0].set_ylabel('Salt transport, S x [Sv]')

    #HEAT MONTHLY
    htd_mean=(-cp*rhow*1e6*WM['TRANS']*WM['PTMP']/tera).sum(dim='WM').groupby('TIME.month').mean(dim='TIME')
    htd_std=(-cp*rhow*1e6*WM['TRANS']*WM['PTMP']/tera).sum(dim='WM').groupby('TIME.month').std(dim='TIME')
    axx[1].fill_between(range(1,13),htd_mean-htd_std,htd_mean+htd_std,color='grey',alpha=0.4)
    htd_mean.plot(label='-Lateral heat transport divergence',ax=axx[1],color='k',linewidth=4)
    hf_mean=hf['NORDIC_hflx'].groupby('time.month').mean(dim='time')
    hf_std=hf['NORDIC_hflx'].groupby('time.month').std(dim='time')
    axx[1].fill_between(range(1,13),hf_mean-hf_std,hf_mean+hf_std,color='limegreen',alpha=0.4)
    axx[1].plot(range(1,13),hf_mean, label='Air-Sea heat flux',color='limegreen',linewidth=4)
    axx[1].set_ylabel('Heat flux [TW]')
    for axi in axx:
        axi.set_xlabel('')
    f.text(0.5, 0, 'Month', ha='center',fontsize=fsz)
    savefig(figdir+'SaltHeat_storage.png',bbox_inches='tight')
    savefig(figdir+'SaltHeat_storage.pdf',bbox_inches='tight')


plot_saltheat_storage()

################################################################################################################################
################################################################################################################################
###########################################      HEAT BALANCE     #######################################################
################################################################################################################################
################################################################################################################################



def plot_tmpfluxes():
    f,axx=subplots(1,2,figsize=(14,4),sharex=True)
    for wm in WM.WM:
        tt_mean=(WM['TRANS']*WM['PTMP']).sel(WM=wm).groupby('TIME.month').mean(dim='TIME')
        tt_std=(WM['TRANS']*WM['PTMP']).sel(WM=wm).groupby('TIME.month').std(dim='TIME')
        axx[0].fill_between(range(1,13),tt_mean-tt_std,tt_mean+tt_std,color=coldic[str(wm.values)],alpha=0.4)
        tt_mean.plot(label=wm.values,ax=axx[0],color=coldic[str(wm.values)],linewidth=4)
    # wm='PWS'
    # (WM['TRANS']*WM['PTMP']).sel(WM=wm).groupby('TIME.month').mean(dim='TIME').plot(label='',ax=axx[0],color=coldic[wm],linewidth=4)
    htd_mean=(-cp*rhow*1e6*WM['TRANS']*WM['PTMP']/tera).sum(dim='WM').groupby('TIME.month').mean(dim='TIME')
    htd_std=(-cp*rhow*1e6*WM['TRANS']*WM['PTMP']/tera).sum(dim='WM').groupby('TIME.month').std(dim='TIME')
    axx[1].fill_between(range(1,13),htd_mean-htd_std,htd_mean+htd_std,color='grey',alpha=0.4)
    htd_mean.plot(label='-Lateral heat transport divergence',ax=axx[1],color='k',linewidth=4)
    hf_mean=hf['NORDIC_hflx'].groupby('time.month').mean(dim='time')
    hf_std=hf['NORDIC_hflx'].groupby('time.month').std(dim='time')
    axx[1].fill_between(range(1,13),hf_mean-hf_std,hf_mean+hf_std,color='limegreen',alpha=0.4)
    axx[1].plot(range(1,13),hf_mean, label='Air-Sea heat flux',color='limegreen',linewidth=4)
    for axi in axx:
        axi.set_xlabel('')
    axx[0].set_ylabel('Temp. transport [$^\circ$C x Sv]')
    axx[1].set_ylabel('Heat flux [TW]')

    axx[0].legend(loc=(0,1.05),ncol=5)
    axx[1].legend()
    axx[0].set_title('')
    xlim(1,12)
    f.text(0.5, 0, 'Month', ha='center',fontsize=14)

    savefig(figdir+'Budget_heat_seas.png',bbox_inches='tight')
    savefig(figdir+'Budget_heat_seas.pdf',bbox_inches='tight')

plot_tmpfluxes()


################################################################################################################################
################################################################################################################################
###########################################      SALT BALANCE     #######################################################
################################################################################################################################
################################################################################################################################
def plot_saltfluxes_seas():
    f,axx=subplots(2,3,figsize=(16,7),sharex=True)
    for gg,wm in enumerate(WM.WM.values):
        if gg>2:
            jj=1
            gg=gg-3
        else:
            jj=0
        for ii in range(19):
            axx[jj,gg].plot(range(1,13),(WM['TRANS']*WM['PSAL']).sel(WM=wm)[12*ii:12*ii+12],'-',color='grey',alpha=0.5)
        (WM['TRANS']*WM['PSAL']).sel(WM=wm).groupby('TIME.month').mean(dim='TIME').plot(color='k',linewidth=4,ax=axx[jj,gg],linestyle='--')
        (WM_obs['TRANS']*WM_obs['PSAL']).sel(WM=wm).groupby('TIME.month').mean(dim='TIME').plot(color=coldic[wm],linewidth=4,ax=axx[jj,gg])
        axx[jj,gg].set_title(wm)

    # axx[0,0].set_ylim(8,22)
    # axx[0,2].set_ylim(-22,-8)
    # axx[1,0].set_ylim(-14,0)
    # axx[1,1].set_ylim(0,14)

    # for ii in range(9):
    #     axx[1,2].plot(range(1,13),WM['TRANS'].sum(dim='WM')[12*ii:12*ii+12],'-',color='grey',alpha=0.5)
    #     axx[1,2].plot(range(1,13),vts['tot'][12*ii:12*ii+12],'-',color='green',alpha=0.5)
    #     # axx[1,2].axhline(WM['TRANS'].sum(dim='WM')[12*ii:12*ii+12].mean(dim='TIME'),color='k',alpha=0.75)
    # WM['TRANS'].sum(dim='WM').groupby('TIME.month').mean(dim='TIME').plot(color='k',linewidth=4,ax=axx[1,2],linestyle='--')
    # # axx[1,2].axhline(vts['tot'].groupby('time.month').mean(dim='time').mean(dim='month'),color='green',linewidth=4)
    # (vts['tot']).groupby('time.month').mean(dim='time').plot(color='green',linewidth=4,ax=axx[1,2])

    for ggg in axx:
        for axi in ggg:
            axi.set_ylabel('')
            axi.set_xlabel('')
            axi.set_xlim(1,12)
    axx[1,2].set_title('Lateral volume transport divergence')
    fsz=14
    f.text(0.075, 0.5, 'Transport [Sv]', va='center',fontsize=fsz,rotation='vertical')
    f.text(0.5, 0.05, 'Month', ha='center',fontsize=fsz)
    # savefig(figdir+'Budget_salinity_seas.png',bbox_inches='tight')
    # savefig(figdir+'Budget_salinity_seas.pdf',bbox_inches='tight')


plot_saltfluxes_seas()


def plot_saltfluxes():
    f,axx=subplots(2,1,figsize=(8,8),sharex=True)
    for wm in WM.WM:
        (WM['TRANS']*WM['PSAL']).sel(WM=wm).plot(label=wm.values,ax=axx[0],color=coldic[str(wm.values)])
    (WM['TRANS']*WM['PSAL']).sum(dim='WM').plot(label='Lateral salinity transport convergence',ax=axx[1],color='k')
    for axi in axx:
        axi.set_ylabel('Salt transport [S x Sv]')
        axi.set_xlabel('')
    axx[0].legend(loc=(1.01,0.2))
    axx[1].legend()
    axx[0].set_title('Salinity budget in NorESM')
    axx[0].set_xlim(datetime.datetime(2010,1,1),datetime.datetime(2019,1,1))
    # savefig(figdir+'Budget_NorESM_salinity.png',bbox_inches='tight')
    # savefig(figdir+'Budget_NorESM_salinity.pdf',bbox_inches='tight')

plot_saltfluxes()
(WM['TRANS']).sum(dim='WM').plot(label='Lateral volume transport convergence',color='k',figsize=(8,4))


corrcoef((WM['TRANS']*WM['PSAL']).sum(dim='WM'),(WM['TRANS']).sum(dim='WM'))
