from aux_funcs import *

datadir

dat=io.loadmat(datadir+'aux_data/EN4/en4irm_surface.mat')

dat.keys()

shape(dat['surface_pden'])

shape(dat['date'])

shape(dat['lon'])

dat['depth']

help(datetime.timedelta)

mkdate=datetime.datetime(2002,1,15) + pd.to_timedelta(np.arange(192), 'M')


en4=xr.Dataset({'pden': (['date', 'lon', 'lat'],  dat['surface_pden']),},
                coords={'lat': dat['lat'].flatten(),
                        'lon': dat['lon'].flatten(),
                        'date': mkdate})

en4

en4.pden.mean(dim='date').T

en4.pden.sel(date=slice('2015-1-1','2015-12-31')).mean(dim='date').T.plot()

en4['date.month']

en4_1416=en4.pden.sel(date=slice('2014-1-1','2016-12-1'))
en4_1416[(en4_1416['date.month']>=1)&(en4_1416['date.month']<=4)].mean(dim='date').T.plot(vmin=27.3,vmax=27.77,cmap=cm.Greens)


en4.to_netcdf(datadir+'aux_data/EN4_surfden.nc','w',format='netCDF4')
