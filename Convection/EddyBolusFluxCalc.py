from aux_funcs import *
dendat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/density_gridded_props_cf5-oom_from5m.nc')

WM={}
WM['upper']=dendat.thickness[:,p1:p2,:].sum(dim='den')
WM['deep']=dendat.thickness[:,p2:p3,:].sum(dim='den')

WM['outcrop_upper']=isnan(dendat.thickness[:,p1,:])
WM['outcrop_deep']=isnan(dendat.thickness[:,p2,:])

hbar={}
for kk in ['upper','deep']:
    hbar[kk]=(WM[kk][2,:]-WM[kk][0,:])/(dendat.distance[2]-dendat.distance[0])/1e3

# for kk in hbar:
#     hbar[kk][WM['outcrop_'+kk][2,:]==1]=NaN

for kk in hbar:
    figure(figsize=(12,3))
    hbar[kk].plot()
    axhline(0,color='k')

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap.pickle','rb'))

# get hprime and vprime for uIIW and dIIW from daily fields, just to see what they are like

vprime={}
vprime['upper']=dat['along track velocity'].where((dat['potential density']<=d2)&(dat['potential density']>d1)).mean(dim='depth')
vprime['deep']=dat['along track velocity'].where((dat['potential density']<=d3)&(dat['potential density']>d2)).mean(dim='depth')

hprime={}
hprime['upper']=dat['depth'].where((dat['potential density']<=d2)&(dat['potential density']>d1)).diff(dim='depth').sum(dim='depth')

for kk in range(4,8):
    (hprime['upper'][kk,:]*vprime['upper'][kk,:]).resample(date='1M').mean(dim='date').plot(figsize=(12,3))

hprime['upper'][4,:].plot(figsize=(12,3))
hprime['upper'][5,:].plot()
hprime['upper'][6,:].plot()
hprime['upper'][7,:].plot()

vprime['upper'][4,:].plot(figsize=(12,3))
axhline(0,color='k')

hprime['upper'][5,:].plot(figsize=(12,3))
vprime['upper'][5,:].plot(figsize=(12,3))
axhline(0,color='k')
vprime['upper'][6,:].plot(figsize=(12,3))
axhline(0,color='k')
