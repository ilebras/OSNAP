from pylab import *
import xarray as xr
from aux_funcs import *

xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Gridded_TS_201408_201604_2018.nc')
xport=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Transports_201408_201604_2018.nc')

xport

xport.MOC_SIGMA.plot()
