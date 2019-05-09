from aux_funcs import *

dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1803extrap.pickle','rb'))

dat.to_netcdf('../pickles/xarray/CFmoorings_notid_1804.nc')
