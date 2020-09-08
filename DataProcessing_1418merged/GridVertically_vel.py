from firstfuncs_1618 import *

###### Load microcat data:
def load_vel(moornum):
    if moornum==8:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/M1_vel_2016recovery_daily.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/M1_vel_2018recovery_daily.nc')
    elif moornum==7:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_vel_2016recovery_daily_dropadcpbin.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_vel_2018recovery_daily_dropadcpbin_drop_u_500.nc')
    else:
            dat16=xr.open_dataset(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(moornum)+'_vel_2016recovery_daily_dropadcpbin.nc')
            dat18=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(moornum)+'_vel_2018recovery_daily_dropadcpbin.nc')

    return dat16,dat18


for moornum in range(1,9):
    d16,d18=load_vel(moornum)
    if 'LATITUDE' in d18:
        d18=d18.drop('LATITUDE').drop('LONGITUDE')

    dat=xr.merge([d16,d18])
    dat
    dat['depth_tvar']=(('TIME','DEPTH'),-gsw.z_from_p(dat['PRES'],mean(dat.LATITUDE.values)))

    ### Plot pressures at which velocity is measured for each mooring for easy reference:
    figure(figsize=(12,3))
    plot(dat.TIME,dat.PRES)
    ylabel('pressure [db]')
    title('CF'+str(moornum))
    savefig(figdir+'pressure_overview_vel/CF'+str(moornum)+'_velpres_dropadcpbin.png')

    #create depth grid:
    dmax=float(dat['depth_tvar'].max())
    dpthgrid=arange(0,dmax,2)
    prsgrid=gsw.p_from_z(-dpthgrid,mean(dat.LATITUDE.values))

    dat_grid=xr.Dataset(coords={'date': dat.TIME.values,'depth':dpthgrid,'lon':dat.LONGITUDE.values,'lat':dat.LATITUDE.values,})
    for var in ['UCUR','VCUR']: # just grid the PSAL and PTMP, extrapolate top 20m gradient, then fill in the derived sal,ptmp, den below
        gridded_var=NaN*ones((len(dat.TIME.values),len(dpthgrid)))
        for tt,dateval in enumerate(dat_grid.date):
            nonan_dat=~isnan(dat[var][tt,:])
            dat_tt_tmp=dat[var][tt,nonan_dat]
            dpth_nonan_tmp=dat.depth_tvar[tt,nonan_dat]
            nanind_dpth=~isnan(dpth_nonan_tmp)
            dat_tt=dat_tt_tmp[nanind_dpth]
            dpth_nonan=dpth_nonan_tmp[nanind_dpth]
            if sum(dat_tt)!=0:
                f1=interpolate.interp1d(dpth_nonan,dat_tt,kind='linear')
                dpth_ind=((dpthgrid<=max(dpth_nonan).values)&(dpthgrid>=min(dpth_nonan).values))
                interdpth=dpthgrid[dpth_ind]
                pan1=f1(interdpth)
                # for velocity, just take the shallowest and deepest recorded value and repeat to fill
                f2=interpolate.interp1d(hstack((0,interdpth,max(dpthgrid))),hstack((pan1[0],pan1,pan1[-1])),kind='linear')
                gridded_var[tt,:]=f2(dpthgrid)

        dat_grid[var]=(('date','depth'),gridded_var)

    print(moornum)
    figure()
    dat_grid.UCUR.T.plot()
    title('CF'+str(moornum)+': w/ dropped bin')
    savefig(figdir+'pressure_overview_vel/CF'+str(moornum)+'_ucur_dropadcpbin.png',bbox_inches='tight')
    figure()
    dat_grid.VCUR.T.plot()
    title('CF'+str(moornum)+': w/ dropped bin')
    savefig(figdir+'pressure_overview_vel/CF'+str(moornum)+'_vcur_dropadcpbin.png',bbox_inches='tight')

    dat_grid.to_netcdf(datadir+'OSNAP_CFgridded_2014-2018/CF'+str(moornum)+'_vel_vertgrid_daily_2m.nc','w',format='netCDF4')
