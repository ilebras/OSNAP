from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'
################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################

osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_1912.nc')
so=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_1912.nc')
hf=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_201001-201812.nc')

WM=xr.open_dataset(datadir+'NorESM/NorESM_WMs_1912.nc')

vts=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_volumetransports_201001-201812.nc')
startyear=2010
startmonth=1
endyear=2018
endmonth=12
vts['time']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}
################################################################################################################################
################################################################################################################################
###########################################      VOL BALANCE     #######################################################
################################################################################################################################
################################################################################################################################


def plot_volfluxes():
    f,axx=subplots(2,1,figsize=(11,6),sharex=True)
    for wm in WM.WM:
        WM['TRANS'].sel(WM=wm.values).plot(label=wm.values,ax=axx[0],color=coldic[str(wm.values)])
    WM['TRANS'].sum(dim='WM').plot(label='Lateral volume transport convergence',ax=axx[1],color='k')
    so['FW+SI'].plot(label='Freshwater volume sources',ax=axx[1],color='c')
    for axi in axx:
        axi.set_ylabel('Transport [Sv]')
        axi.set_xlabel('')
    axx[0].legend(loc=(1.01,0.2))
    axx[1].legend()
    axx[0].set_title('Volume budget in NorESM')
    axx[0].set_xlim(datetime.datetime(2010,1,1),datetime.datetime(2019,1,1))
    savefig(figdir+'Budget_NorESM_volume.png',bbox_inches='tight')
    savefig(figdir+'Budget_NorESM_volume.pdf',bbox_inches='tight')

plot_volfluxes()

################################################################################################################################
################################################################################################################################
###########################################      SALT BALANCE     #######################################################
################################################################################################################################
################################################################################################################################
def plot_saltfluxes():
    f,axx=subplots(2,1,figsize=(11,6),sharex=True)
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
    savefig(figdir+'Budget_NorESM_salinity.png',bbox_inches='tight')
    savefig(figdir+'Budget_NorESM_salinity.pdf',bbox_inches='tight')

plot_saltfluxes()

################################################################################################################################
################################################################################################################################
###########################################      HEAT BALANCE     #######################################################
################################################################################################################################
################################################################################################################################
cp=3850
rhow=1000
tera=10**12

(hf['NORDIC_hflx'].mean(dim='time')*tera/cp/rhow/1e6)

def plot_tmpfluxes():
    f,axx=subplots(2,1,figsize=(11,6),sharex=True)
    for wm in WM.WM:
        (WM['TRANS']*WM['PTMP']).sel(WM=wm).plot(label=wm.values,ax=axx[0],color=coldic[str(wm.values)])
    (-WM['TRANS']*WM['PTMP']).sum(dim='WM').plot(label='-Lateral temperature transport convergence',ax=axx[1],color='k')
    axx[1].plot(osnap['TIME'],hf['NORDIC_hflx']*tera/cp/rhow/1e6, label='Heat flux [$^\circ$C x Sv]',color='limegreen')
    for axi in axx:
        axi.set_ylabel('Temp. transport [$^\circ$C x Sv]')
        axi.set_xlabel('')
    axx[0].legend(loc=(1.01,0.2))
    axx[1].legend()
    axx[0].set_title('Heat budget in NorESM')
    axx[0].set_xlim(datetime.datetime(2010,1,1),datetime.datetime(2019,1,1))
    savefig(figdir+'Budget_NorESM_temp.png',bbox_inches='tight')
    savefig(figdir+'Budget_NorESM_temp.pdf',bbox_inches='tight')

plot_tmpfluxes()
