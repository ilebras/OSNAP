## Save properties in density space as well as gridded density fields for use in other plotting scripts (so I don't have to re-run all this...)
from aux_funcs import *


#mooring resolution gridded CF data
dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1808lpfilt.pickle','rb'))

#load version with no extrapolation and create a mask that can be applied to "dat"
noextrap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap.pickle','rb'))
dat['extrap_mask']=~isnan(noextrap.salinity)

# just work with my "merged" OOI data that goes all the way to the surface and ends up having more comparable resolution to CF/M1 gridded data
oom=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_denmerged_xray.pickle','rb'))
oom['extrap_mask']=~isnan(oom.sal)


##############################################################################################################################
##############################################################################################################################
#################################### First, get CF and merged OOI data on the same grid
##############################################################################################################################
##############################################################################################################################

dmin=10
dmax=2600
gridded={}
grid_dpth=range(dmin,dmax,5)
grid_date=oom.date.values
dat_valvec=['potential density','salinity','temperature','extrap_mask']
oom_valvec=['pden','sal','ptmp','extrap_mask']
for vv in range(len(oom_valvec)):
    print(vv)
    gridded[vv]=zeros((4,len(grid_dpth),len(oom.date)))
    for ii,dd in enumerate(grid_date):
        ## cf+m1
        jj=0
        for mm in range(4,8):
            if mm!=6:
                mden=dat[dat_valvec[vv]][mm,:,:].sel(date=dd)
                if sum(~isnan(mden))>1:
                    nini=~isnan(mden).values
                    nanind=argwhere(nini)[0][0]
                    fm=interp1d(dat['depth'][nini],mden[nini],bounds_error=False)
                    gridded[vv][jj,:,ii]=fm(grid_dpth)
                else:
                    gridded[vv][jj,:,ii]=NaN*ones(len(grid_dpth))
                jj+=1


        ## oom
        oomden=oom[oom_valvec[vv]].sel(date=dd)
        if sum(~isnan(oomden))>1:
            nini=~isnan(oomden).values
            foom=interp1d(oom.depth[nini].values,oomden[nini],bounds_error=False)
            gridded[vv][-1,:,ii]=foom(grid_dpth)
        else:
            gridded[vv][-1,:,ii]=NaN*ones(len(grid_dpth))

datelen=len(grid_date)
for kk in gridded:
    ## now go through every depth level and interpolate horizontally
    for jj,zz in enumerate(grid_dpth):
        for mm in range(4):
            ##ooi
            nani=~isnan(gridded[kk][mm,jj,:])
            if sum(nani)>1:
                f_time=interp1d(arange(datelen)[nani],gridded[kk][mm,jj,nani],bounds_error=False)
                gridded[kk][mm,jj,:]=f_time(arange(datelen))


grdlon=hstack((CFlon[4:6],CFlon[-1],ooi_lon['fla']))
grdlat=hstack((CFlat[4:6],CFlat[-1],ooi_lat['fla']))

grid_dist=hstack(([distvec[4:6],distvec[-1],oom_dist]))

#make it into an xarray!
grdat=xr.Dataset({'pden': (['distance','depth','date'],gridded[0]),'psal': (['distance','depth','date'],gridded[1]),'ptmp': (['distance','depth','date'],gridded[2]),
                  'mask': (['distance','depth','date'],gridded[3]),'moor': (['distance'],['cf5','cf6','m1','ooi_merged']),
                  'lon': (['distance'],grdlon),'lat': (['distance'],grdlat)},
                  coords={'depth': grid_dpth,'date': grid_date,'distance':grid_dist})


grid_mid_depth=grid_dpth[:-1]+diff(grid_dpth)/2
PVmat=zeros((len(grid_dist),len(grid_mid_depth),len(grid_date)))
for ii in range(len(grdat.distance)):
    PVmat[ii,:,:]=sw.f(grdlat[ii])*grdat.pden[ii,:,:].diff(dim='depth')/mean(diff(grdat.depth))/(grdat.pden.values[ii,:-1,:]+grdat.pden[ii,:,:].diff(dim='depth').values/2+1e3)

PVmat[PVmat<=1e-15]=1e-15

#make it into an xarray with PV!!
grdat_PV=xr.Dataset({'pden': (['distance','depth','date'],gridded[0]),'psal': (['distance','depth','date'],gridded[1]),'ptmp': (['distance','depth','date'],gridded[2]),
                  'mask': (['distance','depth','date'],gridded[3]),'PV': (['distance','mid_depth','date'],PVmat),'moor': (['distance'],['cf5','cf6','m1','ooi_merged']),
                  'lon': (['distance'],grdlon),'lat': (['distance'],grdlat)},coords={'depth': grid_dpth,'date': grid_date,'distance':grid_dist,'mid_depth':grid_mid_depth})


##############################################################################################################################
########### Save gridded data xarray as a netcdf file:
grdat_PV.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_5m.nc','w',format='netCDF4')
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
#################################### Next, bin the gridded data in density space
##############################################################################################################################
##############################################################################################################################

denvec=arange(27.5,27.85,0.005)
middenvec=denvec[:-1]+diff(denvec)/2
varvec=['psal','ptmp','dpth','thick','mask']
denmats={}
for var in varvec:
    denmats[var]=zeros((len(grdat.distance),len(middenvec),len(grdat.date)))
    for moornum in range(len(grdat.distance)):
        for dd in range(len(middenvec)):
            if var=='dpth':
                denmats[var][moornum,dd,:]=grdat.depth.where(grdat['pden'][moornum,:,:]>denvec[dd]).min(dim='depth');
                minnan=grdat['pden'][moornum,:,:].min(dim='depth')>denvec[dd] #
                denmats[var][moornum,dd,minnan]=0
            elif var=='thick':
                alldpths=grdat.depth.where(grdat['pden'][moornum,:,:]>denvec[dd]).where(grdat['pden'][moornum,:,:]<=denvec[dd+1])
                denmats[var][moornum,dd,:]=alldpths.max(dim='depth')-alldpths.min(dim='depth');
            elif var=='mask':
                denmats[var][moornum,dd,:]=grdat['mask'][moornum,:,:].where(grdat['pden'][moornum,:,:]>denvec[dd]).where(grdat['pden'][moornum,:,:]<=denvec[dd+1]).min(dim='depth')
            else:
                denmats[var][moornum,dd,:]=grdat[var][moornum,:,:].where(grdat['pden'][moornum,:,:]>denvec[dd]).where(grdat['pden'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');


def makedendat():
    dendat=xr.Dataset({'depth': (['distance','den_bnd','date'],  denmats['dpth']),
                       'thickness': (['distance','den','date'],  denmats['thick']),
                       'psal': (['distance','den','date'],  denmats['psal']),
                       'ptmp': (['distance','den','date'],  denmats['ptmp']),
                       'moor': (['distance'],['cf5','cf6','m1','ooi_merged'])},
                       coords={'den': middenvec,
                               'den_bnd': denvec[:-1],
                               'distance': grdat.distance.values,
                               'date': grdat.date.values})

    return dendat

dendat=makedendat()


##############################################################################################################################
########### Save gridded data xarray binned into density space as a netcdf file:
dendat.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/density_gridded_props_cf5-oom_from5m.nc','w',format='netCDF4')
##############################################################################################################################
