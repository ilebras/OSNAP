from firstfuncs_1618 import *

dat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_2m.nc')


dat['across track velocity']=dat.UCUR*cos(theta)+dat.VCUR*sin(theta)
dat['along track velocity']=-dat.UCUR*sin(theta)+dat.VCUR*cos(theta)

dat=dat.rename({'PSAL':'salinity','PDEN':'potential density','PTMP': 'temperature'})
dat=dat.drop('UCUR').drop('VCUR').drop('ASAL').drop('CTMP')
dat=dat.transpose('distance','depth','date','lat','lon')

dat.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_rotnewfieldnames_2m.nc','w',format='netCDF4')

fcor=sw.f(60)

# 0. Add shelf extrapolation
#first, extrapolate density and salinity using horizontal gradient
shdist=hstack((-10, dat.distance))
shdat=dat.reindex(distance=shdist)
geovel=zeros((len(dat.depth),len(dat.date)))
for vv in ['salinity','potential density']:
        for tt,adate in enumerate(dat.date):
            s1=shdat[vv][1,:,tt]
            s2=shdat[vv][2,:,tt]

            x0=shdist[0]
            x1=shdist[1]
            x2=shdist[2]

            shdat[vv][0,:,tt]=s1-(s2-s1)*(x1-x0)/(x2-x1)

            if 'den' in vv:
                geoshear=(-diff(shdat[vv][:,:,tt],axis=0)[0,:].T*9.8/fcor/1027/diff(shdist)[0].T/1e3).T
                if nansum(geoshear)>0:
                    #ref to bottom
                    geovel[:,tt]=nancumsum(geoshear[::-1])[::-1]*10+dat['across track velocity'][0,:,tt][~isnan(dat['across track velocity'][0,:,tt])][-1].values
#nan out bottom bits of geovel
geovel[isnan(dat.salinity[0,:,100])]=nan



# 1. re-index to finer distance grid
newdist=sort(hstack((-10,-7.5,-5,-2.5,
        dat.distance[:-1]+diff(dat.distance)/4,
        dat.distance[:-1]+diff(dat.distance)/2,
        dat.distance[:-1]+3*diff(dat.distance)/4,
        dat.distance)))


bathatmoor=zeros(len(dat.distance))
#get deepest measurement point at each mooring
for dd in range(len(dat.distance)):
    bathatmoor[dd]=dat.depth[min([where(isnan(dat['salinity'].mean(dim='date').isel(distance=dd))==False)[0][-1],where(isnan(dat['across track velocity'].mean(dim='date').isel(distance=dd))==False)[0][-1]])]

bathatmoor[bathatmoor<200]=200

newdpth=sort(hstack(([bathatmoor[0]]*4,
        bathatmoor[:-1]+diff(bathatmoor)/4,
        bathatmoor[:-1]+diff(bathatmoor)/2,
        bathatmoor[:-1]+3*diff(bathatmoor)/4,
        bathatmoor)))

# pickle.dump([bathatmoor,hstack((newdpth[0],newdpth,newdpth[-1]))],open('../pickles/moordpths.pickle','wb'))

figure()
plot(dat.distance,dat.distance,'o')
plot(newdist,newdist,'x')

figure()
plot(bathdist,bathbath)
plot(newdist,newdpth,'*')
gca().invert_yaxis()


newgrid=shdat.reindex(distance=newdist)

#1.5: Add geostrophic velocity to "across track" vel
newgrid['across track velocity'][0,:,:]=geovel

newgrid['across track velocity'].mean(dim='date').plot()

#2. insert along-bathymetry interpolation to new grid and dont extend downwards -- this was creating a trap!
# Not populating all depths in bottom triangles, so this is currently incomplete.
xtenddown=newgrid.copy()
for vv in dat:
    for tt,adate in enumerate(dat.date):
        bot_y=array([float(dat[vv].isel(date=tt,distance=dd).sel(depth=bathatmoor[dd]).values) for dd in range(len(dat.distance))])
        nanind=~isnan(bot_y)
        if sum(nanind)>1:
            distf=interpolate.interp1d(dat.distance.values[nanind],bot_y[nanind],bounds_error=False)
            for dd in range(len(xtenddown.distance[4:])):
                xtenddown[vv][dd+4,int(argmin(abs(dat.depth-newdpth[dd+4]))),tt]=distf(newgrid.distance[dd+4])

xtenddown['salinity'].mean(dim='date').plot()

xtenddown['across track velocity'].mean(dim='date').plot()

# 3. regrid horizontally (linearly) at finer mesh for each time step
# dold,zold=meshgrid(xtenddown.distance,xtenddown.depth)
newgrid=xtenddown.interpolate_na(dim='distance')
# for vv in dat:
#     for tt,adate in enumerate(dat.date):
#         for zz in range(len(dat.depth)):
#             nonanind=~isnan(xtenddown[vv][:,zz,tt])
#             if sum(nonanind)>1:
#                 dinterp=interpolate.interp1d(newdist[nonanind],xtenddown[vv][:,zz,tt][nonanind],bounds_error=False)
#                 newgrid[vv][:,zz,tt]=dinterp(newdist)

newgrid['salinity'].mean(dim='date').plot()

newgrid['across track velocity'].mean(dim='date').plot()

# 4. now interpolate vertically once more to fill in bottom triangles...
newgrid=newgrid.interpolate_na(dim='depth')
# for vv in dat:
#     if vv[0]!='d':
#         for tt,adate in enumerate(dat.date):
#             for dd in range(len(newgrid.distance)):
#                 nonanind=~isnan(newgrid[vv][dd,:,tt])
#                 if sum(nonanind)>1:
#                     dinterp=interpolate.interp1d(newgrid.depth[nonanind],newgrid[vv][dd,:,tt][nonanind],bounds_error=False)
#                     newgrid[vv][dd,:,tt]=dinterp(newgrid.depth)


newgrid['salinity'].mean(dim='date').plot()

newgrid['across track velocity'].mean(dim='date').plot()


newgrid.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid_2m.nc','w',format='netCDF4')
