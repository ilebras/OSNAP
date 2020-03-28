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


vts=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_volumetransports_201001-201812.nc')
startyear=2010
startmonth=1
endyear=2018
endmonth=12
vts['time']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])


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
        # net_vt_OSNAP.plot(c
    (osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='Southern boundary')
    (fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='Northern boundary')
    (osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='Volume convergence')
    so['FW+SI'].plot(label='Freshwater volume sources',ax=axx[2])
    so['U_storage'].plot(label='Volume storage',ax=axx[2])
    for jj in range(3):
        axx[jj].set_ylabel('Transport [Sv]')
        axx[jj].legend(loc=(1.01,0.2))
        axx[jj].set_xlabel('')
    savefig(figdir+'FW_massbal_fromfields.png',bbox_inches='tight')

plot_volfluxes()

def plot_volfluxes_corr():
    f,axx=subplots(3,1,figsize=(11,8),sharex=True)
    # for ii,sec in enumerate([osnap,ns,fs,bso]):
    #     (sec['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(label=labstr[ii],ax=axx[0])
    vts.net_vt_OSNAP.plot(label=labstr[0],ax=axx[0])#,linestyle='--',color='C0')
    vts.net_vt_NS.plot(label=labstr[1],ax=axx[0])#,linestyle='--',color='C1')
    vts.net_vt_FS.plot(label=labstr[2],ax=axx[0])#,linestyle='--',color='C2')
    vts.net_vt_BSO.plot(label=labstr[3],ax=axx[0])#,linestyle='--',color='C3')
    (vts.net_vt_OSNAP+vts.net_vt_NS).plot(ax=axx[1],label='Southern boundary')
    (vts.net_vt_FS+vts.net_vt_BSO).plot(ax=axx[1],label='Northern boundary')
    (vts.net_vt_OSNAP+vts.net_vt_NS-vts.net_vt_FS-vts.net_vt_BSO).plot(ax=axx[1],label='Volume convergence')
    so['FW+SI'].plot(label='Freshwater volume sources',ax=axx[2])
    # (vts.net_vt_OSNAP+vts.net_vt_NS-vts.net_vt_FS-vts.net_vt_BSO).plot(ax=axx[2],label='Volume convergence')
    (so['U_storage']).plot(label='Volume storage',ax=axx[2])
    for jj in range(3):
        axx[jj].set_ylabel('Transport [Sv]')
        axx[jj].legend(loc=(1.01,0.2))
        axx[jj].set_xlabel('')
    savefig(figdir+'FW_massbal_fromdiag.png',bbox_inches='tight')

plot_volfluxes_corr()


def comp_calc():
    f,axx=subplots(3,1,figsize=(11,8),sharex=True)
    for ii,sec in enumerate([osnap]):
        (sec['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(label='fields',ax=axx[0])
    vts.net_vt_OSNAP.plot(label='diagnostics',ax=axx[0])
    (fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[1],label='fields')
    (vts.net_vt_FS+vts.net_vt_BSO).plot(ax=axx[1],label='diagnostics')
    (osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(ax=axx[2],label='fields')
    (vts.net_vt_OSNAP+vts.net_vt_NS-vts.net_vt_FS-vts.net_vt_BSO).plot(ax=axx[2],label='diagnostics')
    axx[0].set_ylabel('OSNAP Transport [Sv]')
    axx[1].set_ylabel('FS+BSO Transport [Sv]')
    axx[2].set_ylabel('Volume convergence [Sv]')
    axx[2].axhline(0,color='k')
    axx[1].legend(loc=(1.01,0.2))
    for jj in range(3):
        axx[jj].set_xlabel('')

ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE').mean(dim='TIME')

(osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).mean()

(vts.net_vt_OSNAP+vts.net_vt_NS-vts.net_vt_FS-vts.net_vt_BSO).mean()
so['FW+SI'].mean()
((osnap['TRANS']*osnap['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')+(ns['TRANS']*ns['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')-(fs['TRANS']*fs['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')-(bso['TRANS']*bso['PSAL']).sum(dim='DEPTH').sum(dim='LONGITUDE')).mean(dim='TIME')

comp_calc()
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

################################################################################################################################
################################################################################################################################
###########################################      HEAT BALANCE     #######################################################
################################################################################################################################
################################################################################################################################
cp=3850
rhow=1000
tera=10**12

hf['NORDIC_hflx'].plot()

Tconv=((osnap['TRANS']*osnap['PTMP']).sum(dim='DEPTH').sum(dim='LONGITUDE')+(ns['TRANS']*ns['PTMP']).sum(dim='DEPTH').sum(dim='LONGITUDE')-(fs['TRANS']*fs['PTMP']).sum(dim='DEPTH').sum(dim='LONGITUDE')-(bso['TRANS']*bso['PTMP']).sum(dim='DEPTH').sum(dim='LONGITUDE'))
def plot_tmpfluxes():
    f,axx=subplots(2,1,figsize=(11,6),sharex=True)
    for ii,sec in enumerate([osnap,ns,fs,bso]):
        ((sec['TRANS']*sec['PTMP']).sum(dim='DEPTH').sum(dim='LONGITUDE')).plot(label=labstr[ii],ax=axx[0])
    (-Tconv).plot(label='- Temperature transport convergence',ax=axx[1])
    axx[1].plot(osnap['TIME'],hf['NORDIC_hflx']*tera/cp/rhow/1e6, label='Heat flux [$^\circ$C x Sv]')
    for jj in range(2):
        axx[jj].legend(loc=(1.01,0.2))
        axx[jj].set_xlabel('')
    f.text(0.05, 0.5, 'Temperature transport [$^\circ$C x Sv]', va='center', rotation='vertical',fontsize=14)
    savefig(figdir+'FW_tmpbal.png',bbox_inches='tight')

plot_tmpfluxes()
