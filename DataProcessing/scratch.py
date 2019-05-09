
from aux_funcs import *

ctdlist=hstack((glob.glob(datadir+'OSNAP2016recovery/NOC_M1/nocm1_01_2014/microcat/*.nc')))

ctdlist
xr.open_dataset(ctdlist[0])
