from aux_funcs import *

grdat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_10m.nc')

#make sure to work on masked density
grdat['mden']=grdat.pden.where(grdat.mask==True)

for ii in range(4):
    figure()
    grdat.mden[ii,:].plot()
    (grdat['ML'][ii,:]/10).plot()
    axhline(30,color='k')
    axhline(60,color='k')
    axhline(150,color='k')
    ylim(200,0)

grdat.depth
MLmat={}
denthreshvec=arange(0.001,0.01,0.001)
for dd in denthreshvec:
    MLmat[dd]=NaN*zeros((len(grdat.distance),len(grdat.date)))
    for ii in range(len(grdat.distance)):
        for jj in range(len(grdat.date)):
            mind=150
            denvec=grdat.mden[ii,grdat.depth>=mind,jj]
            if ~isnan(denvec[0]):
                diff_from_top=denvec-denvec[0]
                MLmat[dd][ii,jj]=grdat.depth[grdat.depth>=mind][diff_from_top>dd][0]
            # nind=~isnan(denvec)
            # if sum(nind)>=10:
            #     diff_from_top=denvec-denvec[nind][0]
            #     MLmat[dd][ii,jj]=grdat.depth[grdat.depth>=mind][diff_from_top>dd][0]

for ii in range(len(grdat.distance)):
    figure(figsize=(25,4))
    for dd in denthreshvec:
        plot(grdat.date,MLmat[dd][ii,:],label=dd)
        ylim(1500,0)
    plot(grdat.date,MLmat[0.007][ii,:],'k',linewidth=4)
    legend()


##############################################################################################################################
########### Make a decision about which mixed layer to use
grdat['ML']=(['distance','date'],MLmat[0.007])

grdat['ML'][-2,:].plot()

grdat['ML'][-2,:]

##############################################################################################################################
########### Save 25m gridded data xarray binned into density space as a netcdf file:
grdat.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_10m_wML.nc','w',format='netCDF4')
##############################################################################################################################
