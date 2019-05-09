from aux_funcs import *

dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1803extrap.pickle','rb'))

################################################################################
################################################################################
######## Regrid to finer horizotal grid and make nice
################################################################################
################################################################################

# 1. extend all bottommost measurements bottomwards

xtenddown=dat.copy()
for vv in dat:
    if vv[0]!='d':
        for dd in dat.distance:
            for tt in dat.date:
                fieldtmp=dat[vv].sel(distance=dd,date=tt)
                if sum(~isnan(fieldtmp))>0:
                    xtenddown[vv].sel(distance=dd,date=tt)[isnan(fieldtmp)]=fieldtmp[~isnan(fieldtmp)][-1]


# 2. re-index to finer distance grid
newdist=sort(hstack((dat.distance[:-1]+diff(dat.distance)/4,dat.distance[:-1]+diff(dat.distance)/2,dat.distance[:-1]+3*diff(dat.distance)/4,dat.distance)))

figure()
plot(dat.distance,dat.distance,'o')
plot(newdist,newdist,'x')

newgrid=xtenddown.reindex(distance=newdist)


# 3. regrid (linearly) at finer mesh for each time step
dold,zold=meshgrid(xtenddown.distance,xtenddown.depth)
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            newgrid[vv][:,:,tt]=griddata((dold.flatten(),zold.flatten()),xtenddown[vv][:,:,tt].T.data.flatten(),(newgrid.distance.data[None,:],newgrid.depth.data[:,None]),method='linear').T


# 4a. extend distance coord further to include 5 km on either side and repeat outermost measurements on either end of array (for plotting)
distplussides=hstack((-5,newdist,newdist[-1]+5))
plussides=newgrid.reindex(distance=distplussides)
for vv in dat:
    if vv[0]!='d':
        for tt,adate in enumerate(dat.date):
            plussides[vv][0,:,tt]=plussides[vv][1,:,tt]
            plussides[vv][-1,:,tt]=plussides[vv][-2,:,tt]


### Save this version of the xarray (for plotting)
pickle.dump(plussides,open('../pickles/xarray/CF_xarray_gridplot_notid_1803extrap_nobathy.pickle','wb'))

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
