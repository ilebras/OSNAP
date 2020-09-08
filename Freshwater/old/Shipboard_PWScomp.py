from firstfuncs_1618 import *

CTD_14=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_CTD_14.nc')
CTD_16=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_CTD_16.nc')

ADCP_14=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_SADCP_14_10m.nc')
ADCP_16=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_SADCP_16_10m.nc')

#double check all are in the same spot:
plot(CTD_16.lon,CTD_16.lat,'o')
plot(CTD_14.lon,CTD_14.lat,'o')
plot(ADCP_14.lon,ADCP_14.lat,'o')
plot(ADCP_16.lon,ADCP_16.lat,'o')
#I think differences above are because I'm using the longer section for CTD, whereas ADCP also uses inland part when available. Only use overlapping lon ranges!

#calculate PWS salinity in 2014 and 2016
# where is sigma=27.54?
sigmax=27.54
def denplot(CTD_14,ADCP_14):
    CTD_14.sal.plot()
    contour(CTD_14.dist,CTD_14.depth,CTD_14.pden,levels=[sigmax])
    ylim(1000,0)
    dmk=CTD_14.dist.where(CTD_14.lon>-40)
    axvline(dmk[~isnan(dmk)][0])

    figure()
    ADCP_14.uacross.plot()
    contour(CTD_14.dist,CTD_14.depth,CTD_14.pden,levels=[sigmax])
    ylim(1000,0)
    dmk=CTD_14.dist.where(CTD_14.lon>-40)
    distlim=dmk[~isnan(dmk)][0]
    axvline(distlim)
    return distlim

dlim14=denplot(CTD_14,ADCP_14)
dlim16=denplot(CTD_16,ADCP_16)

#get CTD onto ADCP depth grid:
CTD_14_int=CTD_14.interp(depth=ADCP_14.depth)
CTD_16_int=CTD_16.interp(depth=ADCP_16.depth)
#and ADCP onto CTD distance grid:
ADCP_14_int=ADCP_14.interp(dist=CTD_14.dist)
ADCP_16_int=ADCP_16.interp(dist=CTD_16.dist)
All_14=xr.merge([CTD_14_int,ADCP_14_int],compat='override')
All_16=xr.merge([CTD_16_int,ADCP_16_int],compat='override')

All_14.uacross[0,:]=All_14.uacross[1,:]

All_14.uacross.plot()
All_14.sal.plot()

u_14=All_14.uacross.where(All_14.pden<sigmax).where(All_14.dist<dlim14)
s_14=All_14.sal.where(All_14.pden<sigmax).where(All_14.dist<dlim14)

All_16.uacross[:2,:]=All_16.uacross[2,:]
All_16.uacross.plot()
All_16.sal.plot()

u_16=All_16.uacross.where(All_16.pden<sigmax).where(All_16.dist<dlim16)
sal_16=All_16.uacross.where(All_16.pden<sigmax).where(All_16.dist<dlim16)

PWS_S_14=(u_14*All_14.sal).sum(dim='depth').sum(dim='dist')/(u_14).sum(dim='depth').sum(dim='dist')
PWS_S_16=(u_16*All_16.sal).sum(dim='depth').sum(dim='dist')/(u_16).sum(dim='depth').sum(dim='dist')

PWS_S_14
PWS_S_16
