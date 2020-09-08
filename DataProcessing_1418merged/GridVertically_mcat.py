from firstfuncs_1618 import *

###### Load microcat data:
def load_mcat(moornum):
    if moornum==1:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_recon_2016recovery_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_recon_2018recovery_daily.nc')
    elif moornum==2:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_2016recovery_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_recon_2018recovery_daily.nc')
    elif (moornum==4) | (moornum==7):
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_corr_2016recovery_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_corr_2018recovery_daily.nc')
    elif (moornum==5) | (moornum==6):
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_2016recovery_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_corr_2018recovery_daily.nc')
    elif moornum==8:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/M1_mcat_2016recovery_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/M1_mcat_2018recovery_daily.nc')
    else:
        dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_2016recovery_daily.nc')
        dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_mcat_2018recovery_daily.nc')
    return dat16,dat18

def add_SA_CT_PT(xray):
    SA_out=gsw.SA_from_SP(xray['PSAL'].values,prsgrid,xray.lon.values,xray.lat.values)
    CT_out=gsw.CT_from_pt(SA_out,xray['PTMP'].values)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['ASAL']=(('date','depth'),SA_out)
    xray['CTMP']=(('date','depth'),CT_out)
    xray['PDEN']=(('date','depth'),PD_out)


d16,d18=load_mcat(5)

d16.PSAL.sel(DEPTH=750).plot()
d16.PSAL.sel(DEPTH=1000).plot(label='1000m')
d16.PSAL.sel(DEPTH=1300).plot(label='1300m')
legend()

for moornum in range(1,9):
    d16,d18=load_mcat(moornum)
    if 'LATITUDE' in d18:
        d18=d18.drop('LATITUDE').drop('LONGITUDE')

    dat=xr.merge([d16,d18])
    dat['depth_tvar']=(('TIME','DEPTH'),-gsw.z_from_p(dat['PRES'],dat.LATITUDE.values))

    #create depth grid:
    dmax=float(dat['depth_tvar'].max())
    dpthgrid=arange(0,dmax,2)
    prsgrid=gsw.p_from_z(-dpthgrid,dat.LATITUDE.values)

    dat_grid=xr.Dataset(coords={'date': dat.TIME.values,'depth':dpthgrid,'lon':dat.LONGITUDE.values,'lat':dat.LATITUDE.values,})
    for var in ['PSAL','PTMP']: # just grid the PSAL and PTMP, extrapolate top 20m gradient, then fill in the derived sal,ptmp, den below
        gridded_var=NaN*ones((len(dat.TIME.values),len(dpthgrid)))
        for tt,dateval in enumerate(dat_grid.date):
            dat_tt=dat[var][tt,:]
            if sum(~isnan(dat_tt))!=0:
                dpth_nonan=dat.depth_tvar[tt,~isnan(dat_tt)]
                f1=interpolate.interp1d(dpth_nonan,dat_tt[~isnan(dat_tt)],kind='linear')
                dpth_ind=((dpthgrid<=max(dpth_nonan).values)&(dpthgrid>=min(dpth_nonan).values))
                interdpth=dpthgrid[dpth_ind]
                pan1=f1(interdpth)
                # as before, extrapolate top 20m gradient to surface
                z1=dpthgrid[dpthgrid>=min(dpth_nonan).values][0]
                z2=z1+20
                s1=pan1[interdpth==z1]
                s2=pan1[interdpth==z2]
                s0=s1-(s2-s1)*z1/(z2-z1)
                f2=interpolate.interp1d(hstack((0,interdpth,max(dpthgrid))),hstack((s0,pan1,pan1[-1])),kind='linear')
                gridded_var[tt,:]=f2(dpthgrid)
        dat_grid[var]=(('date','depth'),gridded_var)

    print(moornum)
    figure()
    dat_grid.PSAL.plot()
    figure()
    dat_grid.PTMP.plot()

    add_SA_CT_PT(dat_grid)
    figure()
    dat_grid.PDEN.plot()

    dat_grid.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(moornum)+'_mcat_vertgrid_daily_2m.nc','w',format='netCDF4')

dat_grid
