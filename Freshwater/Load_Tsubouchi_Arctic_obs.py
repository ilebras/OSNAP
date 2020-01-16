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
    plot(dat.MLON,dat.MLAT,'x');
    axhline(82)
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
                    'VELO':(['TIME','MLON','DEPTH'],dat.finabsvel.values[:,mlind,:]),
                    'frac':(['TIME','MLON','DEPTH'],dat.frac.values[:,mlind,:]),},
                    coords={'TIME':dat.TIME.values,'DEPTH':dat.DEPTH.values,'LONGITUDE':dat.LONGITUDE.values[latind],'MLON':dat.MLON.values[mlind]})

    sec_final=sec.interp(LONGITUDE=sec.MLON.values).drop('LONGITUDE').rename({'MLON':'LONGITUDE'})
    sec_final=sec_final.assign_coords(LON_LONG=sec.LONGITUDE.values)
    sec_final['LAT_LONG']=('LON_LONG',sec.LATITUDE.values)
    return sec_final

fs=get_section(78,82)
bso=get_section(68,78)


bso.LON_LONG.values

fs=fs.transpose('TIME','DEPTH','LONGITUDE','LON_LONG')
bso=bso.transpose('TIME','DEPTH','LONGITUDE','LON_LONG')

ind=10
plot(fs.LONGITUDE[:ind],fs.VELO.mean('DEPTH').mean('TIME')[:ind],'o')
plot(fs.LONGITUDE[:ind],fs.PSAL.mean('DEPTH').mean('TIME')[:ind],'o')
plot(fs.LON_LONG[:ind],fs.LAT_LONG[:ind],'x')

ind=-20
plot(fs.LONGITUDE[ind:],fs.LATITUDE[ind:],'o')
plot(fs.LON_LONG[ind:],fs.LAT_LONG[ind:],'x')

########################################################################################################
##################################  ADD AREA AND TRANSPORT ################################################
########################################################################################################

def get_AREA_TRANS(sec):
    sec['DISTDIFF']=('LONGITUDE',sw.dist(sec.LAT_LONG.values,sec.LON_LONG.values)[0]*1e3)
    sec['DEPTHDIFF']=('DEPTH',diff(hstack((0,sec.DEPTH.values[:-1]+diff(sec.DEPTH.values)/2,sec.DEPTH.values[-1]+diff(sec.DEPTH.values)[-1]/2))))
    sec['AREA']=sec['DISTDIFF']*sec['DEPTHDIFF']*sec['frac']
    sec['TRANS']=sec['AREA']*sec['VELO']/1e6
    return sec


fs=get_AREA_TRANS(fs)
bso=get_AREA_TRANS(bso)

sec=bso

plot(sec.LAT_LONG.values,sec.LON_LONG.values,'o')

sw.dist(sec.LAT_LONG.values,sec.LON_LONG.values)[0]*1e3

bso.DISTDIFF.plot()

fs.TRANS.sum('DEPTH').sum('LONGITUDE').plot()
bso.TRANS.sum('DEPTH').sum('LONGITUDE').plot()

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

fs.to_netcdf(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_fs_xray_2001.nc')
bso.to_netcdf(datadir+'aux_data/Tsubouchi-etal-2018/Tsubouchi2018_bso_xray_2001.nc')
