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

def load_source_sto():
    int=xr.open_dataset(glob.glob(datadir+'NorESM/*NordicSeas*new.nc')[0])
    int=int.rename({'time': 'TIME'})
    for xx in int:
        int=int.rename({xx : xx[7:]})
    startyear=2010
    startmonth=1
    endyear=2018
    endmonth=12
    int['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return int

int=load_source_sto()

for xx in int:
    figure()
    (int[xx]/1e6).plot()

figdir

(int['saltstorage']/int['volume']).plot()

#######################################################################################
##############    NEW UNITS    ###########################################
###########################################################################
######################################################################################

# convert freshwater sources to Sv...

rho_w=1000

units=int.copy()
for kk in ['iceberg','liqprec','mltfrz','runoff','solprec','evap']:
    if kk=='evap':
        units[kk]=-int[kk]/rho_w/1e6
    else:
        units[kk]=int[kk]/rho_w/1e6

units['FW+SI']=units['iceberg']+units['liqprec']+units['solprec']+units['runoff']+units['mltfrz']+units['evap']
units=units.drop('saltstorage').drop('volume')

def plot_NorESM_FW():
    figure(figsize=(8,3))
    for xx in ['iceberg','liqprec','mltfrz','runoff','solprec','evap','FW+SI']:
        (units[xx]).plot(label=xx)
    ylabel('Transport [Sv]')
    legend(loc=(1.05,0.2))
    savefig(figdir+'FW_sources.png',bbox_inches='tight')

plot_NorESM_FW()
Time_between=int.TIME[:-1]+diff(int.TIME)/2

units=units.assign_coords(TIME_I=Time_between.values)

diff(int.TIME)
plot(diff(int.volume))

plot(diff(int.saltstorage))
plot(diff(int.TIME))
# Get storage terms
units['SU_storage_I']=(['TIME_I'],diff(int.saltstorage)/[float(xx) for xx in diff(int.TIME)]*1e3) # time gives nanoseconds(1e-9), so multiply by 1e3 to get Sv.

units['U_storage_I']=(['TIME_I'],diff(int.volume)/[float(xx) for xx in diff(int.TIME)]*1e3)


units['U_storage']=(['TIME'],units['U_storage_I'].interp(TIME_I=units.TIME))
units['SU_storage']=(['TIME'],units['SU_storage_I'].interp(TIME_I=units.TIME))

units['U_storage_I'].plot()
units['U_storage'].plot()

units['SU_storage_I'].plot()
units['SU_storage'].plot()


units.to_netcdf(datadir+'NorESM/NorESM_source_storage_xray_1912.nc','w')
