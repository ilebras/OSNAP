from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'
################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################
so=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_18yrs_2004.nc')

so['FW+SI'].mean()

startyear=2000
startmonth=1
endyear=2018
endmonth=12

def load_allyears(varlab):
    tmp1=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_200001-200912.nc')[0])
    tmp2=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_201001-201812.nc')[0])
    varout=xr.concat([tmp1,tmp2],dim='time')
    varout=varout.rename({'time':'TIME'})
    varout['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return varout

hf=load_allyears('heatloss')

vts=load_allyears('volumetransports')
sts=load_allyears('salinitytransports')
tts=load_allyears('temperaturetransports')
hts=load_allyears('heattransports')
hf=load_allyears('heatloss')

cp=3850
rhow=1000
################################################################################################################################
###########################################      Check out budgets    #######################################################
################################################################################################################################
WM=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004_sigma.nc')
WM_dpth=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc')
#volume:
latvol=vts.net_vt_OSNAP-vts.net_vt_BSO-vts.net_vt_FS#+vts.net_vt_NS
latvol_WM=WM['TRANS'].sum(dim='WM')
latvol_dpth=WM_dpth['TRANS'].sum(dim='WM')

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM_testing/'


def voltrans_plot():
    figure(figsize=(12,3))
    latvol.plot(label='Convergence from diagnostics = '+str(int(latvol.mean()*1e3))+'mSv')
    latvol_WM.plot(label='Convergence from fields (sigma) = '+str(int(latvol_WM.mean()*1e3))+'mSv')
    latvol_dpth.plot(label='Convergence from fields (depth) = '+str(int(latvol_dpth.mean()*1e3))+'mSv')
    so['FW+SI'].plot(label='FW sources = '+str(int(so['FW+SI'].mean()*1e3))+'mSv')
    legend(loc=(1.01,0.2))
    title('Volume budget')
    savefig(figdir+'volbud_comp.png',bbox_inches='tight')
    savefig(figdir+'volbud_comp.pdf',bbox_inches='tight')

voltrans_plot()

#salinity transports
latsal=sts.net_st_OSNAP-sts.net_st_BSO-sts.net_st_FS#+sts.net_st_NS
latsal_WM=latvol_WM.copy()
latsal_dpth=latvol_WM.copy()

for ii in range(len(latsal_WM.TIME)):
    latsal_WM[ii]=(WM['TRANS'].isel(TIME=ii)*WM['PSAL'].isel(TIME=ii)).sum(dim='WM')
    latsal_dpth[ii]=(WM_dpth['TRANS'].isel(TIME=ii)*WM_dpth['PSAL'].isel(TIME=ii)).sum(dim='WM')


def saltrans_plot():
    figure(figsize=(12,3))
    latsal.plot(label='Convergence from diagnostics = '+str(int(latsal.mean()))+' S x Sv')
    latsal_WM.plot(label='Convergence from fields (sigma) = '+str(int(latsal_WM.mean()))+' S x Sv')
    latsal_dpth.plot(label='Convergence from fields (depth) = '+str(int(latsal_dpth.mean()))+' S x Sv')
    legend(loc=(1.01,0.2))
    title('Salinity budget')
    savefig(figdir+'salbud_comp.png',bbox_inches='tight')
    savefig(figdir+'salbud_comp.pdf',bbox_inches='tight')

saltrans_plot()

#temperature transports
# lattmp=tts.net_ht_OSNAP-tts.net_ht_BSO-tts.net_ht_FS
latheat=(hts.net_ht_OSNAP-hts.net_ht_BSO-hts.net_ht_FS)

lattmp_WM=latvol_WM.copy()
lattmp_dpth=latvol_WM.copy()
for ii in range(len(latsal_WM.TIME)):
    lattmp_WM[ii]=((WM['TRANS'].isel(TIME=ii)*WM['PTMP'].isel(TIME=ii)).sum(dim='WM'))*cp*rhow/1e6
    lattmp_dpth[ii]=((WM_dpth['TRANS'].isel(TIME=ii)*WM_dpth['PTMP'].isel(TIME=ii)).sum(dim='WM'))*cp*rhow/1e6

def heattrans_plot():
    figure(figsize=(12,3))
    (latheat/1e12).plot(label='Convergence from diagnostics = '+str(int(latheat.mean()/1e12))+' TW')
    (lattmp_WM).plot(label='Convergence from fields (sigma) = '+str(int(lattmp_WM.mean()))+' TW')
    (lattmp_dpth).plot(label='Convergence from fields (depth) = '+str(int(lattmp_dpth.mean()))+' TW')
    (-hf.NORDIC_hflx).plot(label='Surface heat flux = '+str(int(hf.NORDIC_hflx.mean()))+' TW')
    legend(loc=(1.01,0.2))
    title('Heat budget')
    savefig(figdir+'heatbud_comp.png',bbox_inches='tight')
    savefig(figdir+'heatbud_comp.pdf',bbox_inches='tight')

heattrans_plot()

lattmp.mean()

lattmp.plot()
lattmp_WM.plot()
hf_unit=(hf.NORDIC_hflx/cp/rhow*1e6)
hf_unit.plot()
lattmp.mean()
lattmp_WM.mean()
hf.NORDIC_hflx.mean()
(hf_unit).mean()
WM_OSNAP=WM.sel(WM=slice('AWS','DWS'))
lattmp_WM_OSNAP=lattmp_WM.copy()
for ii in range(len(latsal_WM.TIME)):
    lattmp_WM_OSNAP[ii]=(WM_OSNAP['TRANS'].isel(TIME=ii)*WM_OSNAP['PTMP'].isel(TIME=ii)).sum(dim='WM')

tts.net_ht_OSNAP.plot()
(tts.net_ht_OSNAP+tts.net_ht_NS).plot()
lattmp_WM_OSNAP.plot()

WM_NORTH=WM.sel(WM=slice('AWN','PWN'))
lattmp_WM_NORTH=lattmp_WM.copy()
for ii in range(len(latsal_WM.TIME)):
    lattmp_WM_NORTH[ii]=(WM_NORTH['TRANS'].isel(TIME=ii)*WM_NORTH['PTMP'].isel(TIME=ii)).sum(dim='WM')

(-tts.net_ht_BSO-tts.net_ht_FS).plot()
(lattmp_WM_NORTH).plot()
