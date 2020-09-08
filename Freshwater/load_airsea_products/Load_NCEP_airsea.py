from firstfuncs_1618 import *

## CFSv1 (2000-2011)


## CFSv2 (2011-)
#land mask
lmask=xr.open_dataset(datadir+'aux_data/Air-Sea_FW/NCEP_CFSv2/LandCover/pgbh.gdas.201101.grb2.nc')
lmask['mask']=(('lat','lon'),lmask.LAND_L1_Avg.values[0,:,:])
lmask.mask.plot()
#eminusp
eplist=sort(glob.glob(datadir+'aux_data/Air-Sea_FW/NCEP_CFSv2/Average/EminusP/*.nc'))
ep=xr.open_dataset(eplist[0])
for ii,xx in enumerate(eplist[1:]):
    eptmp=xr.open_dataset(xx)
    ep=xr.concat([ep,eptmp],dim='time')
ep=ep.rename({'EMNP_L1_FcstAvg6hr':'ep'})
#surfacefluxes
hflist=sort(glob.glob(datadir+'aux_data/Air-Sea_FW/NCEP_CFSv2/Average/SurfaceFluxes/*.nc'))
hf=xr.open_dataset(hflist[0])
for ii,xx in enumerate(hflist[1:]):
    hftmp=xr.open_dataset(xx)
    hf=xr.concat([hf,hftmp],dim='time')


hf=hf.rename({'SHTFL_L1_FcstAvg6hr':'shf','LHTFL_L1_FcstAvg6hr':'lhf','DSWRF_L1_FcstAvg6hr':'dsw','USWRF_L1_FcstAvg6hr':'usw','DLWRF_L1_FcstAvg6hr':'dlw','ULWRF_L1_FcstAvg6hr':'ulw'})


hf=hf.where(lmask.mask==0)

hf.usw.mean(dim='time').plot()
hf.dsw.mean(dim='time').plot()
hf.lhf.mean(dim='time').plot()
hf.shf.mean(dim='time').plot()
ep.ep.mean(dim='time').plot()
#load NORESM boundaries and extract data within box
osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')

bnd_lon=hstack((fs.LONGITUDE.values,bso.LONGITUDE.values,ns.LONGITUDE.values[::-1],osnap.LONGITUDE.values[::-1],-45,-40,fs.LONGITUDE.values[0]))
bnd_lat=hstack((fs.LATITUDE.values,bso.LATITUDE.values,ns.LATITUDE.values[::-1],osnap.LATITUDE.values[::-1],60,75,fs.LATITUDE.values[0]))

plot(bnd_lon,bnd_lat,linewidth=3,color='k')

import regionmask

bnd_tpl=array(tuple((zip(bnd_lon,bnd_lat))))

mybox = regionmask.Regions([bnd_tpl], name='mybox')
ax = mybox.plot()

#create masks
epmask = mybox.mask(ep.lon.values,ep.lat.values)
hfmask = mybox.mask(hf.lon.values,hf.lat.values)

cut_ep=ep.where(epmask==0)
cut_hf=hf.where(hfmask==0)

cut_ep

cut_ep['time']=cut_hf['time']

cut_ep.ep.mean(dim='time').plot()
cut_hf.lhf.mean(dim='time').plot()
cut_hf.shf.mean(dim='time').plot()


def get_area(dat):
    lonmat,latmat=meshgrid(dat.lon,dat.lat)
    loni=vstack((lonmat.T-0.25,lonmat[:,-1].T+0.25)).T
    lathold=vstack((latmat.T,latmat[:,-1].T)).T
    lati=vstack((latmat+0.25,latmat[-1,:]+0.25))
    lonhold=vstack((latmat,latmat[-1,:]))
    # plot(lonmat[:10,:10],latmat[:10,:10],'.');
    # plot(loni[:10,:10],lati[:10,:10],'o');
    latdist=gsw.distance(lonhold,lati,axis=0)
    londist=gsw.distance(loni,lathold)
    dat['area']=(('lat','lon'),londist*latdist)
    return dat


cut_ep=get_area(cut_ep)
cut_hf=get_area(cut_hf)



cut_ep['ep_int']=(cut_ep.ep*cut_ep.area).sum(dim='lon').sum(dim='lat')/100/24/60/60/1e6
#convert from cm/day to Sv (10^6 m3/s)


cut_hf['hf_int']=((cut_hf.lhf+cut_hf.shf+cut_hf.ulw+cut_hf.usw-cut_hf.dsw-cut_hf.dlw)*cut_hf.area).sum(dim='lon').sum(dim='lat')/1e12
#convert to TW


cut_hf.area.plot()

int_all=xr.merge([cut_hf.hf_int,cut_ep.ep_int])
int_all.ep_int.plot()


int_all.to_netcdf(datadir+'aux_data/Air-Sea_FW/NCEP_CFSv2/NCEP_ephf_cutint_ilebras.nc','w')
