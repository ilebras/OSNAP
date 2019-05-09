from aux_funcs import *

nclist=glob.glob('/media/isabela/Seagate Backup Plus Drive/NARR/*.nc')

#For now, extract along ~60N from -45 to -40 ?

#Use initial_time0_hours to get time from each, seems easiest to use. Build array that is lon x time.
minlat=59.80
maxlat=60.1
minlon=-43.0
maxlon=-40.5

dat=Dataset(nclist[0])

locind=(dat['gridlat_221'][:]>minlat)&(dat['gridlat_221'][:]<=maxlat)&(dat['gridlon_221'][:]>minlon)&(dat['gridlon_221'][:]<=maxlon)

lon=dat['gridlon_221'][:][locind]
lat=dat['gridlat_221'][:][locind]

plot(CFlon,CFlat,'o-',label='OSNAP mooring positions')
plot(lon,lat,'o-',label='Wind speed positions')
legend()
savefig('../figures/wind/windpositions.png')

for na,dd in enumerate(sort(nclist)):
    dat=Dataset(dd)
    uwindtmp=dat['U_GRD_221_HTGL'][:][:,locind]
    vwindtmp=dat['V_GRD_221_HTGL'][:][:,locind]
    hourtmp=dat['initial_time0_hours'][:]

    if na==0:
        uwind=uwindtmp.copy()
        vwind=vwindtmp.copy()
        hours=hourtmp.copy()
    else:
        uwind=vstack((uwind,uwindtmp))
        vwind=vstack((vwind,vwindtmp))
        hours=hstack((hours,hourtmp))


date=[datetime.datetime(1800,1,1)+datetime.timedelta(hours=hh) for hh in hours]

figure(figsize=(14,3))
plot(date,uwind);
figure(figsize=(14,3))
plot(date,vwind);

distvec=cumsum(hstack((0,sw.dist(lat,lon)[0])))


#create a wind xarray that can be loaded elsewhere

winddat=xr.Dataset({'u wind speed': (['date', 'distance'],  uwind),
                'v wind speed': (['date','distance'],  vwind),
                'longitude':(['distance'],lon),
                'latitude':(['distance'],lat)},
                coords={'distance': distvec,
                        'date': date})

theta
cos(theta)
winddat['across track wind speed']=winddat['v wind speed']*sin(theta)+winddat['u wind speed']*cos(theta)

winddat['along track wind speed']=winddat['v wind speed']*cos(theta)-winddat['u wind speed']*sin(theta)

pickle.dump(winddat,open('../pickles/wind/NARR_linet_newtheta.pickle','wb'))
