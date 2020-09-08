from firstfuncs_1618 import *

dat=xr.open_dataset('../data/OSNAP2016recovery/pickles/xarray/CFmoorings_notid_1804.nc')

dat

zeta_test=0.6/1e4
zeta_test
sw.f(60)
zeta_test/sw.f(60)

#doing a rough overview of PV components with gridded daily fields. can go into instrument data later.

## q/f = (1 + zeta/f) b_z - (u_z^2 + v_z^2)

# calculate each of these terms at mid-depth points.
dat=dat.assign_coords(mid_depth=dat.depth.values[:-1]+diff(dat.depth.values)/2)

#b_z=g/rho_0(-drho/dz) Using g/rho0=100
dat['bz']=(('distance','mid_depth','date'),(dat['potential density'].diff(dim='depth')/dat.depth.diff(dim='depth')).values/100)
dat['uz2']=(('distance','mid_depth','date'),((dat['across track velocity'].diff(dim='depth')/dat.depth.diff(dim='depth')).values)**2)
dat['vz2']=(('distance','mid_depth','date'),((dat['along track velocity'].diff(dim='depth')/dat.depth.diff(dim='depth')).values)**2)

dat=dat.assign_coords(mid_distance=dat.distance.values[:-1]+diff(dat.distance.values)/2)

dat['zeta_f']=(('mid_distance','depth','date'),(dat['across track velocity'].diff(dim='distance')/dat.distance.diff(dim='distance')/1e3/sw.f(60)).values)

def plot_zeta_f():
    (dat.zeta_f.mean(dim='date').T).plot()
    contour(dat.distance.values,dat.depth.values,dat['across track velocity'].mean(dim='date').values.T,colors='k',levels=arange(-0.4,0.8,0.05))
    ylim(2000,0)
    ylabel('distance [km]')
    text(45,-20,'CF5')
    text(60,-20,'CF6')
    savefig(figdir+'SI/zeta_f_mean.png',bbox_inches='tight')

plot_zeta_f()

figdir
#### plot time series of PV (components at moorings CF4-CF6)
def plot_pv_comps(moornum):
    f,axx=subplots(3,1,sharex=True,sharey=True,figsize=(10,10))
    dlim=150
    vmaxi=0.5
    (dat.bz[moornum,dlim:,:]*1e5).plot(ax=axx[0],vmin=0,vmax=vmaxi,cbar_kwargs={"label":'$b_z$'})
    ((dat.uz2[moornum,dlim:,:]+dat.vz2[moornum,:,:])*1e5).plot(ax=axx[1],vmin=0,vmax=vmaxi,cbar_kwargs={"label":'$u_z^2+v_z^2$'})
    ((dat.bz[moornum,dlim:,:]-dat.uz2[moornum,dlim:,:]-dat.vz2[moornum,dlim:,:])*1e5).plot(ax=axx[2],vmin=0,vmax=vmaxi,cbar_kwargs={"label":'$b_z-u_z^2-v_z^2$'})
    ylim(1500,0)
    for xx in axx:
        xx.set_xlabel('')
        xx.set_ylabel('depth [m]')
        xx.set_title('')
    savefig(figdir+'SI/PV_hovmueller_CF'+str(moornum+1)+'.png',bbox_inches='tight')

    #plot time series of PV at 625m.
    dind=int(624/2)
    figure(figsize=(12,3))
    plot(dat.date,(dat.bz[moornum,dind,:])*1e5,label='$b_z$',color='black')
    plot(dat.date,(dat.uz2[moornum,dind,:]+dat.vz2[moornum,dind,:])*1e5,label='$u_z^2+v_z^2$',color='C1')
    plot(dat.date,(dat.bz[moornum,dind,:]-dat.uz2[moornum,dind,:]-dat.vz2[moornum,dind,:])*1e5,label='$b_z-(u_z^2+v_z^2)$',color='red')
    legend()
    xlim([datetime.datetime(2014,8,15),datetime.datetime(2016,8,15)])
    savefig(figdir+'SI/PV_tseries_CF'+str(moornum+1)+'_'+str(dind*2)+'m.png',bbox_inches='tight')
    xlim([datetime.datetime(2015,2,1),datetime.datetime(2015,7,1)])
    savefig(figdir+'SI/PV_tseries_CF'+str(moornum+1)+'_'+str(dind*2)+'m_zoom.png',bbox_inches='tight')


plot_pv_comps(4)
plot_pv_comps(5)


era=xr.open_dataset('../data/aux_data/ERA_1804/tau_hflux_180411.nc')

era=era.rename({'inss':'tauy','iews':'taux','longitude':'lon','latitude':'lat','time':'date','sshf':'sensible','slhf':'latent'})
era['lon']=era['lon']-360

era['tau along']=-(era['taux']*cos(theta)+era['tauy']*sin(theta))
era['heat flux']=era['sensible']+era['latent']

era['tau along']
