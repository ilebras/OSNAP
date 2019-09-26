from aux_funcs import *

dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1810JHIL.pickle','rb'))

dat

#make a dictionary of vars and save as a mat for Alex (masters student)

datadic={}
datadic['lon']=CFlon
datadic['lat']=CFlat
datadic['depth']=dat.depth.values[::5]
datadic['ptmp']=dat.temperature.values[:,::5,:]
datadic['sal']=dat.salinity.values[:,::5,:]
datadic['u_alongstream']=dat['across track velocity'].values[:,::5,:]
datadic['u_acrossstream']=dat['along track velocity'].values[:,::5,:]
datadic['date']=arange(datetime2matlabdn(dt.datetime(2014,8,17)),datetime2matlabdn(dt.datetime(2014,8,17))+712)

io.savemat(datadir+'OSNAP_CFmoor_1416_daily.mat',datadic)

XXXXXXXXXXXXXXXXXXXXX
################################################################################
######## Test run horizontal extrap of sal and den
################################################################################
#
# saltest=dat['potential density'].mean(dim='date')
#
# newsal=saltest.reindex(distance=hstack((-20,saltest.distance)))
#
# s1=newsal[1,:]
# s2=newsal[2,:]
#
# x0=newsal.distance[0]
# x1=newsal.distance[1]
# x2=newsal.distance[2]
#
# newsal[0,:]=s1-(s2-s1)*(x1-x0)/(x2-x1)
#
# fcor=fcor=gsw.f(60)
# geoshear=(-diff(newsal,axis=0).T*9.8/fcor/1028/diff(newsal.distance).T/1e3).T
#
# geovel=nancumsum(geoshear[:,::-1],axis=1)[:,::-1]*2
# geodist=newsal.distance[:-1]+diff(newsal.distance)/2
#
# newvel=dat['across track velocity'].mean(dim='date').reindex(distance=hstack((-10,saltest.distance)))
# newvel[0,:]=geovel[0,:]+newvel[1,:][~isnan(newvel[1,:])][-1].values
#
#
# figure()
# contourf(geodist,newsal.depth,geovel.T,vmin=-0.7,vmax=0)
# ylim([80,0])
# xlim([-10,90])
# colorbar()
# title('Geostrophic vel referenced to zero at the bottom')
# figure()
# contourf(dat.distance,newsal.depth,dat['across track velocity'].mean(dim='date').T,vmin=-0.7,vmax=0)
# ylim([80,0])
# xlim([-10,90])
# colorbar()
# title('Directly measured velocity')
# figure()
# contourf(hstack((-10,saltest.distance)),newsal.depth,newvel.T,vmin=-0.7,vmax=0)
# ylim([80,0])
# xlim([-10,90])
# colorbar()
# title('Geostrophic vel ref to bottom CF1 vel fills shelf')
#
#
# def qplot(saltest):
#     figure()
#     contourf(saltest.distance,saltest.depth,saltest.T,vmin=26,vmax=28)
#     ylim([80,0])
#     xlim([-20,90])
#     colorbar()
#
#
# qplot(saltest)
# qplot(newsal)

################################################################################
################################################################################
######## Regrid to finer horizotal grid and make nice
################################################################################
################################################################################


#first, fill in one day that is missing from M1

tprob=where(isnan(dat['salinity'].isel(distance=-1).mean(dim='depth'))==True)[0][0]
for vv in dat:
    if vv[0]!='d':
        dat[vv][-1,:,tprob]=(dat[vv].isel(date=tprob+1,distance=-1)+dat[vv].isel(date=tprob-1,distance=-1))/2

# 0. Add shelf extrapolation
#first, extrapolate density and salinity using horizontal gradient
shdist=hstack((-10, dat.distance))
shdat=dat.reindex(distance=shdist)
geovel=zeros((len(dat.depth),len(dat.date)))
for vv in ['potential density','salinity']:
        for tt,adate in enumerate(dat.date):
            s1=shdat[vv][1,:,tt]
            s2=shdat[vv][2,:,tt]

            x0=shdist[0]
            x1=shdist[1]
            x2=shdist[2]

            shdat[vv][0,:,tt]=s1-(s2-s1)*(x1-x0)/(x2-x1)

            if 'den' in vv:
                geoshear=(-diff(shdat[vv][:,:,tt],axis=0)[0,:].T*9.8/fcor/1027/diff(shdist)[0].T/1e3).T
                #ref to bottom
                geovel[:,tt]=nancumsum(geoshear[::-1])[::-1]*2+dat['across track velocity'][0,:,tt][~isnan(dat['across track velocity'][0,:,tt])][-1].values
                #ref to surface
                # geovel_inter=nancumsum(geoshear[::-1])[::-1]*2
                # geovel[:,tt]=geovel_inter+(dat['across track velocity'][0,0,tt].values-geovel_inter[0])
#nan out bottom bits of geovel
geovel[isnan(dat.salinity[0,:,0])]=nan

# 1. re-index to finer distance grid
newdist=sort(hstack((-10,-7.5,-5,-2.5,
        dat.distance[:-1]+diff(dat.distance)/4,
        dat.distance[:-1]+diff(dat.distance)/2,
        dat.distance[:-1]+3*diff(dat.distance)/4,
        dat.distance)))


bathatmoor=zeros(len(dat.distance))
#get deepest measurement point at each mooring
for dd in range(len(dat.distance)):
    bathatmoor[dd]=dat.depth[where(isnan(dat['salinity'].mean(dim='date').isel(distance=dd))==False)[0][-1]]

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

newgrid.distance[4]

#1.5: Add geostrophic velocity to "across track" vel
newgrid['across track velocity'][0,:,:]=geovel


newgrid['across track velocity'].mean(dim='date').plot()

#2. insert along-bathymetry interpolation to new grid and dont extend downwards -- this was creating a trap!
# Not populating all depths in bottom triangles, so this is currently incomplete.
xtenddown=newgrid.copy()
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            distf=interpolate.interp1d(dat.distance,
                  [dat[vv].isel(date=tt,distance=dd,
                  depth=int(argmin(abs(dat.depth-bathatmoor[dd])))) for dd in range(len(dat.distance))])
            for dd in range(len(xtenddown.distance[4:])):
                xtenddown[vv][dd+4,int(argmin(abs(dat.depth-newdpth[dd+4]))),tt]=distf(newgrid.distance[dd+4])
xtenddown['across track velocity'].mean(dim='date').plot()

# 3. regrid horizontally (linearly) at finer mesh for each time step
# dold,zold=meshgrid(xtenddown.distance,xtenddown.depth)
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            for zz in range(len(dat.depth)):
                nonanind=~isnan(xtenddown[vv][:,zz,tt])
                if sum(nonanind)>1:
                    dinterp=interpolate.interp1d(newdist[nonanind],xtenddown[vv][:,zz,tt][nonanind],bounds_error=False)
                    newgrid[vv][:,zz,tt]=dinterp(newdist)

newgrid['across track velocity'].mean(dim='date').plot()

# 4. now interpolate vertically once more to fill in bottom triangles...
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            for dd in range(len(newgrid.distance)):
                nonanind=~isnan(newgrid[vv][dd,:,tt])
                if sum(nonanind)>1:
                    dinterp=interpolate.interp1d(newgrid.depth[nonanind],newgrid[vv][dd,:,tt][nonanind],bounds_error=False)
                    newgrid[vv][dd,:,tt]=dinterp(newgrid.depth)



## Dont do this last step since I'm extrapolating better on shelf and its not really necessary on the other side.

# 4a. extend distance coord further to include 5 km on either side and repeat outermost measurements on either end of array (for plotting)
# distplussides=hstack((-5,newdist,newdist[-1]+5))
# plussides=newgrid.reindex(distance=distplussides)
# for vv in dat:

#     if vv[0]!='d':
#         for tt,adate in enumerate(dat.date):
#             plussides[vv][0,:,tt]=plussides[vv][1,:,tt]
#             plussides[vv][-1,:,tt]=plussides[vv][-2,:,tt]
# dat['salinity'].mean(dim='date').plot()
#
# plussides['salinity'].mean(dim='date').plot()

### Save this version of the xarray (for plotting)
pickle.dump(newgrid,open('../pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','wb'))


newgrid['across track velocity'].mean(dim='date').plot()

## Took this part out since gridcalc version was not loading, it was taking up room, and I can just mask locally
#
# # 4b. mask the fields based on bathymetry file (this version does not have extra fields on either side, thats just for plotting)
# bathf=interpolate.interp1d(bathdist,bathbath)
# bathonmygrid=bathf(newgrid.distance)
#
# figure()
# plot(bathdist,bathbath,'-')
# plot(newgrid.distance,bathonmygrid,'o')
#
# maskbathy=newgrid.copy()
# for vv in dat:
#     if vv[0]!='d':
#         for dd,adist in enumerate(maskbathy.distance):
#             maskbathy[vv][dd,:,:]=maskbathy[vv][dd,:,:].where(maskbathy[vv][dd,:,:].depth<=bathonmygrid[dd])
#
# ### Save this version of the xarray (for calculating transports)
# pickle.dump(maskbathy,open('../pickles/CF_xarray_gridcalc_notid.pickle','wb'))
