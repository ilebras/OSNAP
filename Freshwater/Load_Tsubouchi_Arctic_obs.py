from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/???/'

################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################
datlist=glob.glob(datadir+'aux_data/Tsubouchi-etal-2018/section*')

def load_convert():
    dat=xr.open_dataset(sort(datlist)[0])
    timvec=dat.data_time
    for ii,dd in enumerate(sort(datlist)[1:]):
        dat_tmp=xr.open_dataset(dd)
        # print(ii,dat_tmp.data_time)
        dat=xr.concat([dat,dat_tmp],dim='data_time')
        timvec=hstack((timvec,dat_tmp.data_time))

    timvec=[datetime.datetime(2005,9,15)+datetime.timedelta(days=30.5*ii) for ii in range(12)]

    dat=dat.rename({'data_time':'TIME'})
    dat['TIME']=timvec

    dat=dat.rename({'npres':'DEPTH'})
    dat=dat.assign_coords(DEPTH=(['DEPTH'],-gsw.z_from_p(dat.pres[0,:],75)))
    dat=dat.rename({'ngrid':'LONGITUDE'})
    dat=dat.assign_coords(LONGITUDE=(['LONGITUDE'],dat.lon[0,:]))
    dat['LATITUDE']=(['LONGITUDE'],dat.lat[0,:])
    dat=dat.rename({'npair':'MLON'})
    dat=dat.assign_coords(MLON=(['MLON'],dat.mlon[0,:]))
    dat['MLAT']=(['MLON'],dat.mlat[0,:])
    dat=dat.drop('lat').drop('lon').drop('pres').drop('mlon').drop('mlat')

    dat=dat.rename({'ptmp':'PTMP','sal':'PSAL','temp':'TEMP'})
    return dat

dat=load_convert()

def map_split():
    plot(dat.LONGITUDE,dat.LATITUDE,'o');
    axhline(78)
    axhline(68)

map_split()


def get_section(minlat,maxlat):
    latind=(dat.LATITUDE>minlat)& (dat.LATITUDE<=maxlat)
    mlind=(dat.MLAT>minlat)& (dat.MLAT<=maxlat)
    sec=xr.Dataset({'PTMP':(['TIME','LONGITUDE','DEPTH'],dat.PTMP.values[:,latind,:]),
                    'PSAL':(['TIME','LONGITUDE','DEPTH'],dat.PSAL.values[:,latind,:]),
                    'LATITUDE':(['LONGITUDE'],dat.LATITUDE.values[latind]),
                    'MLAT':(['MLON'],dat.MLAT.values[mlind]),
                    'VELO_I':(['TIME','MLON','DEPTH'],dat.finabsvel.values[:,mlind,:]),
                    'frac':(['TIME','MLON','DEPTH'],dat.frac.values[:,mlind,:]),},
                    coords={'TIME':dat.TIME.values,'DEPTH':dat.DEPTH.values,'LONGITUDE':dat.LONGITUDE.values[latind],'MLON':dat.MLON.values[mlind]})
    return sec

fs=get_section(78,82)
fs.PSAL.mean(dim='TIME').plot()
fs['VELO_I']=

bso=get_section(68,78)
bso.PSAL.mean(dim='TIME').plot()

fs=fs.transpose('TIME','DEPTH','LONGITUDE','MLON')
bso=bso.transpose('TIME','DEPTH','LONGITUDE','MLON')

def add_PDEN(xray):
    if 'PRES' in list(xray.data_vars):
            PRES_out=xray['PRES']
    else:
            PRES_out=gsw.p_from_z(-xray['DEPTH'],60)
    SA_out=NaN*xray['PSAL'].copy()

    for ii,pp in enumerate(PRES_out):
        for jj,ll in enumerate(xray.LONGITUDE.values):
            SA_out[:,ii,jj]=gsw.SA_from_SP(xray['PSAL'][:,ii,jj],pp,ll,xray.LATITUDE[jj])
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['PDEN']=(('TIME','DEPTH','LONGITUDE'),PD_out)
    return xray

fs=add_PDEN(fs)
bso=add_PDEN(bso)

fs



fs.to_netcdf(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_fs_xray_1912.nc')
bso.to_netcdf(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_bso_xray_1912.nc')
