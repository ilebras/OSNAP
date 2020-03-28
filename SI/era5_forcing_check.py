from firstfuncs_1618 import *
figdir

dat=xr.open_dataset('../data/aux_data/ERA_1804/tau_hflux_180411.nc')

dat.UCUR*cos(theta)+dat.VCUR*sin(theta)

era=dat.rename({'inss':'tauy','iews':'taux','longitude':'lon','latitude':'lat','time':'date','sshf':'sensible','slhf':'latent'})
era['lon']=era['lon']-360

era['tau along']=-(era['taux']*cos(theta)+era['tauy']*sin(theta))
era['heat flux']=era['sensible']+era['latent']

def pltseries(var):
    figure(figsize=(12,3))
    era[var].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(color='k')
    xlim([datetime.datetime(2014,8,15),datetime.datetime(2016,8,15)])
    title('')
    axhline(0,color='k')
    savefig(figdir+'SI/Forcing_tseries_'+var+'.png',bbox_inches='tight')
    # xlim([datetime.datetime(2015,2,1),datetime.datetime(2015,7,1)])
    # savefig(figdir+'SI/Forcing_tseries_'+var+'_zoom.png',bbox_inches='tight')

pltseries('taux')
pltseries('tauy')
pltseries('heat flux')
era
era['heat flux'].mean(dim='date').plot()
