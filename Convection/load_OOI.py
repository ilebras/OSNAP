from pylab import *
from netCDF4 import Dataset
import glob
import xarray as xr
import datetime as dt
%matplotlib inline

datadir='/home/isabela/Documents/projects/OSNAP/data/'

# #quick plot
# def checkout_dataset(LIST):
#     figure(figsize=(8,8))
#     for ff in LIST:
#         print(ff)
#         FLA=xr.open_dataset(ff)
#         salind=FLA.practical_salinity>5
#         subplot(211)
#         plot(FLA.practical_salinity[salind],FLA.ctdmo_seawater_temperature[salind],'o');
#         subplot(212)
#         plot(FLA.time,FLA.ctdmo_seawater_pressure);
#
# checkout_dataset(FLA_list)

def np64ToDatetime(DA):
    return array([datetime.datetime.utcfromtimestamp((dd-np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for dd in DA])

FLAlist_d1=glob.glob(datadir+'aux_data/OOI/FLA/*t0001*.nc')
FLAlist_d2=glob.glob(datadir+'aux_data/OOI/FLA/*t0002*.nc')

FLBlist_d1=glob.glob(datadir+'aux_data/OOI/FLB/*t0001*.nc')
FLBlist_d2=glob.glob(datadir+'aux_data/OOI/FLB/*t0002*FLMB*.nc')

SUMO_list=glob.glob(datadir+'aux_data/OOI/SUMO/*.nc')
# SUMO_list

#re-organize and save as pickle for use elsewhere
def loadset(datlist):
    for ii,dd in enumerate(datlist):
        dat=xr.open_dataset(dd)

        Stmp=dat['practical_salinity'].values
        Ttmp=dat['ctdmo_seawater_temperature'].values
        Ptmp=dat['practical_salinity'].values

        plevtmp=int(mean(dat.ctdmo_seawater_pressure/10))*10

        salind=dat.practical_salinity<5

        Stmp[salind]=NaN
        Ttmp[salind]=NaN
        Ptmp[salind]=NaN

        if ii==0:
            sal=Stmp
            tmp=Ttmp
            prs=Ptmp
            plev=plevtmp
        else:
            sal=vstack((sal,Stmp))
            tmp=vstack((tmp,Ttmp))
            prs=vstack((prs,Ptmp))
            plev=hstack((plev,plevtmp))

    # print(dat.time.values)
    date_in_dt=np64ToDatetime(dat.time.values)
    # print(date_in_dt)

    dset=xr.Dataset({'salinity': (['depth','date'], sal),
                     'temperature': (['depth','date'], tmp),
                     'pressure': (['depth','date'], prs)},
                     coords={'date': date_in_dt, 'depth': plev})

    dset_sorted=dset.sortby(dset.depth)

    return dset_sorted





FLA_d1=loadset(FLAlist_d1)
FLA_d2=loadset(FLAlist_d2)

FLB_d1=loadset(FLBlist_d1)
FLB_d2=loadset(FLBlist_d2)

import pickle

FLA=xr.merge([FLA_d1,FLA_d2])
FLB=xr.merge([FLB_d1,FLB_d2])
pickle.dump([FLA,FLB],open(datadir+'OSNAP2016recovery/pickles/OOI/FLMAB_xrays.pickle','wb'),protocol=2)
