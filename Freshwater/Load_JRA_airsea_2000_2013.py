from firstfuncs_1618 import *

def load_and_rename(varshort,jraname,newname,lnum='4'):
    evaplist=sort(glob.glob(datadir+'aux_data/Air-Sea_FW/JRA55/*_'+varshort+'*.nc'))
    evap=xr.open_dataset(evaplist[0])
    for ii,xx in enumerate(evaplist[1:]):
        eptmp=xr.open_dataset(xx)
        evap=xr.concat([evap,eptmp],dim='initial_time0_hours')
    evap=evap.rename({jraname:newname,'initial_time0_hours':'time','g'+lnum+'_lon_2':'lon','g'+lnum+'_lat_1':'lat'})
    evap['time']=array([tt+np.timedelta64(14,'D') for tt in evap['time'].values])
    lind=where(evap.lon>180)[0][-1]
    evap['lon']=hstack((evap['lon'][:lind+1]-360,evap['lon'][lind+1:]))
    evap[newname].mean(dim='time').plot()
    return evap


evap=load_and_rename('evp','EVP_GDS4_SFC_S130','evap')
prec=load_and_rename('tprat','TPRAT_GDS4_SFC_S130','prec')

sens=load_and_rename('shtfl','SHTFL_GDS4_SFC_S130','sens')
late=load_and_rename('lhtfl','LHTFL_GDS4_SFC_S130','late')

downlong=load_and_rename('dlwrf','DLWRF_GDS4_SFC_S130','downlong')
uplong=load_and_rename('ulwrf','ULWRF_GDS4_SFC_S130','uplong')

upshort=load_and_rename('uswrf','USWRF_GDS4_SFC_S130','upshort')

downshort=load_and_rename('dswrf','DSWRF_GDS4_SFC_S130','downshort')

ep=evap.evap-prec.prec
hf=sens.sens+late.late+upshort.upshort+uplong.uplong-downshort.downshort-downlong.downlong

lmask=load_and_rename('land','LAND_GDS4_SFC','lm')


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

cut_ep.mean(dim='time').plot()
cut_hf.mean(dim='time').plot()

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

cut_ep.area.plot()

cut_ep
noresm=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_1912.nc')
noresm['ep_int']=noresm.evap-noresm.liqprec-noresm.solprec
cut_ep['ep_int']=(cut_ep*cut_ep.area*(1-lmask.lm[0,:,:].values)).sum(dim='lon').sum(dim='lat')/1000/24/60/60/1e6
#convert from mm/day to Sv (10^6 m3/s)

(1-lmask.lm[0,:,:]).plot()

figure(figsize=(14,3))
cut_ep.ep_int.plot()
noresm.ep_int.plot()

cut_hf['hf_int']=(cut_hf*cut_hf.area*(1-lmask.lm[0,:,:])).sum(dim='lon').sum(dim='lat')/1e12
#convert to TW
cut_hf.hf_int.plot()

int_all=xr.Dataset({'hf_int':(['time'],cut_hf.hf_int.values),'ep_int':(['time'],cut_ep.ep_int.values)},coords={'time':cut_hf.time.values},)

int_all.to_netcdf(datadir+'aux_data/Air-Sea_FW/JRA55/JRA55_2000_2013_ephf_cutint_ilebras.nc','w')
