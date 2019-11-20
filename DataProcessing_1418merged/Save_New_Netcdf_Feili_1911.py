from firstfuncs_1618 import *

import datetime as datetime

#####################################################################################
######## Go through moorings one by one and create new files as needed #############
#####################################################################################
#####################################################################################

def toTimestamp(d):
  return calendar.timegm(d.timetuple())

def to_datetime_andordinal(date):
    """
    Converts a numpy datetime64 object to a python datetime object
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    dtimeobj=datetime.datetime.utcfromtimestamp(timestamp)
    return datetime.datetime.toordinal(dtimeobj)

def rm_coords(elvario):
    for kk in elvario:
        if 'coordinates' in elvario[kk].attrs:
            del elvario[kk].attrs['coordinates']


#####################################################################################
####################################### CF1 #########################################
#####################################################################################

moornum=1

dlist=glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*')

cf1=xr.open_dataset(dlist[0])
for dd in dlist[1:]:
    if 'Nov2019' not in dd:
        tmp=xr.open_dataset(dd)
        cf1=xr.concat([cf1,tmp],dim='DEPTH')

mtime=cf1.resample(TIME='1M').mean(dim='TIME').TIME
time_month=[to_datetime_andordinal(ddd) for ddd in mtime.values]
time_all=[to_datetime_andordinal(ddd) for ddd in cf1.TIME.values]

sal_mdiff=(cf1.PSAL.sel(DEPTH=50)-cf1.PSAL.sel(DEPTH=100)).resample(TIME='1M').mean(dim='TIME').values

sal_mdiff_fill_100=linspace(sal_mdiff[8],sal_mdiff[0],5)
sal_mdiff_int=hstack((sal_mdiff[:9],sal_mdiff_fill_100[1:-1],sal_mdiff[:9],sal_mdiff_fill_100[1:-1],sal_mdiff[:2]))
f100=interp1d(time_month,sal_mdiff_int,bounds_error=False)
sal_mdiff_fulltime=f100(time_all)

psal_r50=cf1.PSAL.sel(DEPTH=100)+sal_mdiff_fulltime

cf1_50_corr=cf1.sel(DEPTH=50).copy()
cf1_50_corr['PRES'][isnan(cf1['PSAL'].sel(DEPTH=50))]=mean(cf1['PRES'].sel(DEPTH=50)).values
cf1_50_corr['PSAL_QC'][isnan(cf1.sel(DEPTH=50)['PSAL'])]=8
cf1_50_corr['PSAL'][isnan(cf1.sel(DEPTH=50).PSAL)]=psal_r50[isnan(cf1.sel(DEPTH=50).PSAL)].values

cf1_50_corr['TEMP_QC'][isnan(cf1.sel(DEPTH=50)['TEMP'])]=8
cf1_50_corr['TEMP'][isnan(cf1['TEMP'].sel(DEPTH=50))]=cf1.TEMP.sel(DEPTH=100)[isnan(cf1['TEMP'].sel(DEPTH=50))]+mean(cf1.TEMP.sel(DEPTH=50)-cf1.TEMP.sel(DEPTH=100))
#
cf1_50_corr.attrs['comment']=cf1.attrs['comment']+'; missing data filled in by extrapolation from below instruments (flagged as 8).'

rm_coords(cf1_50_corr)
cf1_50_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF1_MCTD_2016_14713_50m_Nov2019.nc','w',format='netCDF4')

####################################################################################
###################################### CF2 #########################################
####################################################################################
moornum=2

dlist=sort(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*'))
cf2=xr.open_dataset(dlist[0])
for dd in dlist[1:]:
    if 'Nov2019' not in dd:
        tmp=xr.open_dataset(dd)
        cf2=xr.concat([cf2,tmp],dim='DEPTH')


mtime=cf2.resample(TIME='1M').mean(dim='TIME').TIME
time_month=[to_datetime_andordinal(ddd) for ddd in mtime.values]
time_all=[to_datetime_andordinal(ddd) for ddd in cf2.TIME.values]

cf2_100_corr=cf2.sel(DEPTH=100).copy()

cf2_100_corr['PRES'][isnan(cf2['PSAL'].sel(DEPTH=100))]=mean(cf2['PRES'].sel(DEPTH=100)).values
cf2_100_corr['PSAL_QC'][isnan(cf2.sel(DEPTH=100)['PSAL'])]=8
cf2_100_corr['PSAL'][isnan(cf2.sel(DEPTH=100).PSAL)]=cf2.PSAL.sel(DEPTH=200)[isnan(cf2['PSAL'].sel(DEPTH=100))]+mean(cf2.PSAL.sel(DEPTH=100)-cf2.PSAL.sel(DEPTH=200))
cf2_100_corr['TEMP_QC'][isnan(cf2.sel(DEPTH=100)['TEMP'])]=8
cf2_100_corr['TEMP'][isnan(cf2['TEMP'].sel(DEPTH=100))]=cf2.TEMP.sel(DEPTH=200)[isnan(cf2['TEMP'].sel(DEPTH=100))]+mean(cf2.TEMP.sel(DEPTH=100)-cf2.TEMP.sel(DEPTH=200))
# cf2_100_corr.TEMP.plot()
#delete "coordinates" attribute from each... (ts just "time", which is in there a few times)
rm_coords(cf2_100_corr)
cf2_100_corr.attrs['comment']=cf1.attrs['comment']+'; missing data filled in by extrapolation from below instruments (flagged as 8).'
cf2_100_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF2_MCTD_2016_14631_100m_Nov2019.nc','w',format='netCDF4')

cf2_50_corr=cf2.sel(DEPTH=100).copy()
cf2_50_corr['DEPTH']=50
cf2_50_corr['PRES']=50*cf2['PRES'].sel(DEPTH=100)/cf2['PRES'].sel(DEPTH=100)
cf2_50_corr['PSAL_QC']=8*cf2['PRES'].sel(DEPTH=100)/cf2['PRES'].sel(DEPTH=100)
cf2_50_corr['PSAL']=cf2.PSAL.sel(DEPTH=100)-0.35968952 # Copied from CF2_merged_reconstruction.py
cf2_50_corr['TEMP_QC']=8*cf2['PRES'].sel(DEPTH=100)/cf2['PRES'].sel(DEPTH=100)
cf2_50_corr['TEMP']=cf2.TEMP.sel(DEPTH=100)-0.8 # Approx value copied from CF2_merged_reconstruction.py
cf2_50_corr.attrs['instrument']='instrument lost: this is a reconstruction from 100m instrument, based on difference between 50m and 100m instrument during the first deployment'
cf2_50_corr.attrs['comment']=cf1.attrs['comment']+'; missing data filled in by extrapolation from below instruments (flagged as 8).'
rm_coords(cf2_50_corr)


cf2_50_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF2_MCTD_2016_recon_50m_Nov2019.nc','w',format='netCDF4')

#####################################################################################
####################################### CF4 #########################################
#####################################################################################
moornum=4

cf4_350_corr=xr.open_dataset(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*350m*p.nc')[0])

cf4=cf4_350_corr.copy() #just lazy coding
cf4_350_corr['PSAL']=NaN*cf4.PSAL
cf4_350_corr['PSAL_QC']=5*cf4.PSAL/cf4.PSAL
cf4_350_corr.attrs['comment']=cf4.attrs['comment']+'; salinity was found to have uncorrectable drifts: removed and flagged here.'

cf4_350_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF4_MCTD_2016_14642_350m_adcp_Nov2019.nc','w',format='netCDF4')


#####################################################################################
####################################### CF6 #########################################
#####################################################################################
moornum=6
cf6_1500=xr.open_dataset(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*1500m.nc')[0])
cf6_bot_corr=xr.open_dataset(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*BOTTOM.nc')[0])

cf6_bot_corr['PSAL']=(['TIME'],cf6_1500['PSAL'].values)
cf6_bot_corr['PSAL_QC']=8*cf6_bot_corr.PSAL/cf6_bot_corr.PSAL
cf6_bot_corr.attrs['comment']=cf4.attrs['comment']+'; salinity was found to have uncorrectable drifts: replaced with 1500m salinity -- use caution when interpreting.'


cf6_bot_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF6_MCTD_2016_5932_BOTTOM_Nov2019.nc','w',format='netCDF4')

#####################################################################################
####################################### CF7 #########################################
#####################################################################################
moornum=7

cf7_1500=xr.open_dataset(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*1500m.nc')[0])
cf7_bot_corr=xr.open_dataset(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*BOTTOM.nc')[0])

cf7_bot_corr['PSAL']=(['TIME'],cf7_1500['PSAL'].values)
cf7_bot_corr['PSAL_QC']=8*cf7_bot_corr.PSAL/cf7_bot_corr.PSAL
cf7_bot_corr.attrs['comment']=cf4.attrs['comment']+'; salinity was found to have uncorrectable drifts: replaced with 1500m salinity -- use caution when interpreting.'

cf7_bot_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF7_MCTD_2016_3588_BOTTOM_Nov2019.nc','w',format='netCDF4')

cf7_100_corr=xr.open_dataset(glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*100m.nc')[0])

cf7_100_corr['PSAL']=NaN*cf7_100_corr.PSAL
cf7_100_corr['PSAL_QC']=5*cf7_100_corr.TEMP/cf7_100_corr.TEMP
cf7_100_corr.attrs['comment']=cf4.attrs['comment']+'; salinity was found to have uncorrectable drifts: removed and flagged here.'

cf7_100_corr.to_netcdf(datadir+'OSNAP2018recovery/mcat_nc/CF7_MCTD_2016_14614_100m_Nov2019.nc','w',format='netCDF4')


#####################################################################################
# Finally, load all altered files, change some of the descriptive attributes and re-save
corrdatlist=glob.glob(datadir+'OSNAP2018recovery/mcat_nc/*Nov2019.nc')
for kk in corrdatlist:
    tmp=xr.open_dataset(kk)
    tmp.attrs['date_created']='2019-11-13'
    tmp.attrs['product_version']='v2'
    tmp.attrs['history']='Processed by Jamie Holte, Scripps Institution of Oceanography, March 2019; Additional processing by Isabela Le Bras, Nov 2019'
    tmp.to_netcdf(kk[:-3]+'v2'+'.nc','w',format='netCDF4')
