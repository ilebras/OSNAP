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

################################################################################################################################
################################################################################################################################
###########################################      VOL BALANCE     #######################################################
################################################################################################################################
################################################################################################################################
labstr=['OSNAP','North Seas','Fram Strait', 'Barents Sea Opening']
def plot_volfluxes():
    f,axx=subplots(3,1,figsize=(11,8),sharex=True)
    for ii,sec in enumerate([osnap,ns,fs,bso]):
        (sec['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(label=labstr[ii],ax=axx[0])
    (osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='Southern boundary')
    (fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='Northern boundary')
    (osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='Volume convergence')
    so['FW+SI'].plot(label='Freshwater volume sources',ax=axx[2])
    so['U_storage'].plot(label='Volume storage',ax=axx[2])
    for jj in range(3):
        axx[jj].set_ylabel('Transport [Sv]')
        axx[jj].legend(loc=(1.01,0.2))
        axx[jj].set_xlabel('')
    savefig(figdir+'FW_massbal.png',bbox_inches='tight')


plot_volfluxes()

################################################################################################################################
################################################################################################################################
###########################################      SALT BALANCE     #######################################################
################################################################################################################################
################################################################################################################################
def plot_saltfluxes():
    f,axx=subplots(2,1,figsize=(11,6),sharex=True)
    for ii,sec in enumerate([osnap,ns,fs,bso]):
        ((sec['TRANS']*sec['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(label=labstr[ii],ax=axx[0])
    ((osnap['TRANS']*osnap['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')+(ns['TRANS']*ns['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')
    -(fs['TRANS']*fs['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')-(bso['TRANS']*bso['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(label='Salt transport convergence',ax=axx[1])
    (-so['SU_storage']/10).plot(label='Salt volume storage (/10)',ax=axx[1])
    for jj in range(2):
        axx[jj].set_ylabel('Salt transport [S x Sv]')
        axx[jj].legend(loc=(1.01,0.2))
        axx[jj].set_xlabel('')
    savefig(figdir+'FW_saltbal.png',bbox_inches='tight')

plot_saltfluxes()
