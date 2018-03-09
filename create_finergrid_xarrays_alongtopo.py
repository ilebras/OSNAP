from aux_funcs import *


dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1801.pickle','rb'))

################################################################################
################################################################################
######## Regrid to finer horizotal grid and make nice
################################################################################
################################################################################
dat

#0. fill in one day that is missing from M1
tprob=where(isnan(dat['salinity'].isel(distance=-1).mean(dim='depth'))==True)[0][0]
tprob
for vv in dat:
    if vv[0]!='d':
        dat[vv][-1,:,tprob]=(dat[vv].isel(date=tprob+1,distance=-1)+dat[vv].isel(date=tprob-1,distance=-1))/2

# 1. re-index to finer distance grid
newdist=sort(hstack((dat.distance[:-1]+diff(dat.distance)/4,
        dat.distance[:-1]+diff(dat.distance)/2,
        dat.distance[:-1]+3*diff(dat.distance)/4,
        dat.distance)))

bathf=interpolate.interp1d(bathdist,bathbath)
bathonmygrid=bathf(newdist)
bathatmoor=zeros(len(dat.distance))
#get deepest measurement point at each mooring
for dd in range(len(dat.distance)):
    bathatmoor[dd]=dat.depth[where(isnan(dat['salinity'].mean(dim='date').isel(distance=dd))==False)[0][-1]]

newdpth=sort(hstack((bathatmoor[:-1]+diff(bathatmoor)/4,
        bathatmoor[:-1]+diff(bathatmoor)/2,
        bathatmoor[:-1]+3*diff(bathatmoor)/4,
        bathatmoor)))

pickle.dump([bathatmoor,hstack((newdpth[0],newdpth,newdpth[-1]))],open('../pickles/moordpths.pickle','wb'))

figure()
plot(dat.distance,dat.distance,'o')
plot(newdist,newdist,'x')

figure()
plot(bathdist,bathbath)
plot(newdist,newdpth,'*')
gca().invert_yaxis()

newgrid=dat.reindex(distance=newdist)

#2. insert along-bathymetry interpolation to new grid and then extend downwards
xtenddown=newgrid.copy()
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            distf=interpolate.interp1d(dat.distance,
                  [dat[vv].isel(date=tt,distance=dd,
                  depth=int(argmin(abs(dat.depth-bathatmoor[dd])))) for dd in range(len(dat.distance))])
            botfld=distf(newgrid.distance)
            for dd in range(len(xtenddown.distance)):
                xtenddown[vv][dd,int(argmin(abs(dat.depth-newdpth[dd]))):,tt]=botfld[dd]


# 3. regrid horizontally (linearly) at finer mesh for each time step
dold,zold=meshgrid(xtenddown.distance,xtenddown.depth)
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            for zz in range(len(dat.depth)):
                nonanind=where(isnan(newgrid[vv][:,zz,tt])==False)[0]
                dinterp=interpolate.interp1d(newdist[nonanind],xtenddown[vv][:,zz,tt][nonanind])
                newgrid[vv][:,zz,tt]=dinterp(newdist)

newgrid['salinity'].mean(dim='date').plot()

# 4a. extend distance coord further to include 5 km on either side and repeat outermost measurements on either end of array (for plotting)
distplussides=hstack((-5,newdist,newdist[-1]+5))
plussides=newgrid.reindex(distance=distplussides)
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            plussides[vv][0,:,tt]=plussides[vv][1,:,tt]
            plussides[vv][-1,:,tt]=plussides[vv][-2,:,tt]
dat['salinity'].mean(dim='date').plot()

plussides['salinity'].mean(dim='date').plot()

### Save this version of the xarray (for plotting)
pickle.dump(plussides,open('../pickles/xarray/CF_xarray_gridplot_notid_1803bathy.pickle','wb'))



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
