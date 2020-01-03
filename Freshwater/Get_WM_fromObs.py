from firstfuncs_1618 import *

osnap_obs=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full_wMOC.nc')
fs=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_bso_xray_1912.nc')


AWS_obs={}
PWS_obs={}
DWS_obs={}
AWN_obs={}
PWN_obs={}

xray=osnap_obs
sigmax_obs=osnap_obs.SIGMEAN.values
oslim_obs=-41
for var in ['PSAL','PTMP','PDEN','TRANS']:
    if var=='TRANS':
        DWS_obs[var]=xray['TRANS'].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS_obs[var]=xray['TRANS'].where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS_obs[var]=xray['TRANS'].where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
    else:
        DWS_obs[var]=(xray[var]*xray['TRANS']).where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.PDEN>=sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        PWS_obs[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE<oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWS_obs[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=oslim_obs).where(xray.PDEN<sigmax_obs).sum(dim='DEPTH').sum(dim='LONGITUDE')


#split up FS into PW and AW using definition in Tsubouchi et al. 2018, which is boundary between "EGC" + "Middle", 2W
fslim=-2

fs

for var in ['PSAL','PTMP','PDEN','TRANS']:
    if var=='TRANS':
        PWN_obs[var]=-fs['TRANS'].where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs[var]=-fs['TRANS'].where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')
    else:
        PWN_obs[var]=(fs[var]*fs['TRANS']).where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')/fs['TRANS'].where(fs.LONGITUDE<fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')
        AWN_obs[var]=((fs[var]*fs['TRANS']).where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')+(bso[var]*bso['TRANS']).sum(dim='DEPTH').sum(dim='LONGITUDE'))/(fs['TRANS'].where(fs.LONGITUDE>=fslim).sum(dim='DEPTH').sum(dim='LONGITUDE')+bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE'))





WM_obs={}
for ii,xrw in enumerate([PWS_obs,AWS_obs,DWS_obs,PWN_obs,AWN_obs]):
    WM_obs[ii]=xr.concat([xrw['PDEN'],xrw['PSAL'],xrw['PTMP'],xrw['TRANS']],pd.Index(['PDEN','PSAL','PTMP','TRANS'],name='var'))

WMall_obs=xr.concat([WM_obs[ww] for ww in WM_obs],pd.Index(['PWS','AWS','DWS'],name='WM'))

WMall_obs=WMall_obs.to_dataset('var')
WMall_obs.to_netcdf(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc','w')
