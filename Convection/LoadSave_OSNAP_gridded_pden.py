from pylab import *
import xarray as xr
from aux_funcs import *

TS=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Gridded_TS_201408_201604_2018.nc')
V=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Gridded_V_201408_201604_2018.nc')

dat=xr.merge([TS,V])

bathy=io.loadmat(datadir+'OSNAPbathy_Fli.mat')

SAdat=gsw.SA_from_SP(dat.PSAL,0,-45,57)
dat['PTMP']=dat.TEMP.copy()
dat['PDEN']=dat.TEMP.copy()
for ii in range(len(dat.TIME)):
    for jj in range(len(dat.LONGITUDE)):
        SA=gsw.SA_from_SP(dat['PSAL'][ii,:,jj],gsw.p_from_z(-dat.DEPTH,dat.LATITUDE[jj]),dat.LONGITUDE[jj],dat.LATITUDE[jj])
        dat['PTMP'][ii,:,jj]=gsw.pt0_from_t(SA,dat.TEMP[ii,:,jj],gsw.p_from_z(-dat.DEPTH,dat.LATITUDE[jj]))
        CT=gsw.CT_from_pt(SA,dat['PTMP'][ii,:,jj])
        dat['PDEN'][ii,:,jj]=gsw.sigma0(SA,CT)

pickle.dump(dat,open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','wb'))
