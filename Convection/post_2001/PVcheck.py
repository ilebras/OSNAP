from firstfuncs_1618 import *

moornum=5

vel=xr.open_dataset(datadir+'OSNAP2016recovery/Hourly_netcdf/CF'+str(moornum)+'_vel_2016recovery_hourly_dropadcpbin.nc')

mcat=xr.open_dataset(datadir+'OSNAP2016recovery/Hourly_netcdf/CF'+str(moornum)+'_mcat_2016recovery_hourly.nc')
vel['ug']=vel.UCUR*cos(theta)+vel.VCUR*sin(theta)

vel['z']=(['TIME','DEPTH'],gsw.z_from_p(vel['PRES'],60))

mcat['z']=(['TIME','DEPTH'],gsw.z_from_p(mcat['PRES'],60))

figdir='/home/isabela/Documents/projects/OSNAP/figures/MixedLayer/Slantwise_Convection/'



def plot_depths():
    figure(figsize=(15,3))
    vel['z'].sel(DEPTH=250).plot()
    vel['z'].sel(DEPTH=500).plot()
    vel['z'].sel(DEPTH=750).plot()
    vel['z'].sel(DEPTH=1000).plot()
    mcat['z'].sel(DEPTH=250).plot()
    mcat['z'].sel(DEPTH=500).plot()
    mcat['z'].sel(DEPTH=750).plot()
    mcat['z'].sel(DEPTH=1000).plot()
    title('CF5: instrument positions')
    savefig(figdir+'Inst_pos.png',bbox_inches='tight',dpi=300)



plot_depths()


def plot_PVcomps_2day(d1,d2):
    f,axx=subplots(2,1,figsize=(15,8),sharex=True)
    bz=((mcat['PDEN'].sel(DEPTH=d2)-mcat['PDEN'].sel(DEPTH=d1))/(-mcat['z'].sel(DEPTH=d2)+mcat['z'].sel(DEPTH=d1))/100)*10**5
    vz2=((vel['ug'].sel(DEPTH=d1)-vel['ug'].sel(DEPTH=d2))**2/(-vel['z'].sel(DEPTH=d2)+vel['z'].sel(DEPTH=d1))**2)*10**5
    # bz.plot(ax=axx[0],alpha=0.75,color='limegreen',label='$\Delta b/\Delta z$',linewidth=0,marker='.',)
    # vz2.plot(ax=axx[0],alpha=0.75,color='C1',label='$(\Delta v/\Delta z)^2$',linewidth=0,marker='.',)
    bz.resample(TIME='2D').mean().plot(ax=axx[0],color='limegreen',label='$\Delta b/\Delta z$',linewidth=2)
    vz2.resample(TIME='2D').mean().plot(ax=axx[0],color='C1',label='$(\Delta v/\Delta z)^2$',linewidth=2)
    # (bz-vz2).plot(linewidth=0,marker='.',ax=axx[1],color='yellow',label='')
    # (bz-vz2).where(bz>0).plot(linewidth=0,marker='.',ax=axx[1],color='purple',label='$\Delta b/\Delta z-(\Delta v/\Delta z)^2$',)
    (bz-vz2).resample(TIME='2D').mean().plot(ax=axx[1],color='purple',label='$\Delta b/\Delta z-(\Delta v/\Delta z)^2$',linewidth=2)
    axhline(0,color='r',zorder=10)
    for axi in axx:
        axi.set_xlabel('')
        axi.set_title('')
        axi.set_ylabel('x 10$^{-5}$ [1/s$^2$]')
        axi.legend()
        axi.set_ylim(-0.25,1)
    axx[0].set_title('CF5: '+str(d2)+'m - '+str(d1)+'m')
    xlim(datetime.datetime(2014,9,1),datetime.datetime(2016,8,1))
    savefig(figdir+'PV_'+str(d2)+'-'+str(d1)+'_twodaymeans.png',bbox_inches='tight',dpi=300)


def plot_PVcomps(d1,d2):
    f,axx=subplots(2,1,figsize=(15,8),sharex=True)
    bz=((mcat['PDEN'].sel(DEPTH=d2)-mcat['PDEN'].sel(DEPTH=d1))/(-mcat['z'].sel(DEPTH=d2)+mcat['z'].sel(DEPTH=d1))/100)*10**5
    vz2=((vel['ug'].sel(DEPTH=d1)-vel['ug'].sel(DEPTH=d2))**2/(-vel['z'].sel(DEPTH=d2)+vel['z'].sel(DEPTH=d1))**2)*10**5
    bz.plot(ax=axx[0],alpha=0.75,color='limegreen',label='$\Delta b/\Delta z$',linewidth=0,marker='.',)
    vz2.plot(ax=axx[0],alpha=0.75,color='C1',label='$(\Delta v/\Delta z)^2$',linewidth=0,marker='.',)
    bz.resample(TIME='2D').mean().plot(ax=axx[0],color='grey',label='$\Delta b/\Delta z$',linewidth=2)
    vz2.resample(TIME='2D').mean().plot(ax=axx[0],color='r',label='$(\Delta v/\Delta z)^2$',linewidth=2)
    (bz-vz2).plot(linewidth=0,marker='.',ax=axx[1],color='yellow',label='')
    (bz-vz2).where(bz>0).plot(linewidth=0,marker='.',ax=axx[1],color='purple',label='$\Delta b/\Delta z-(\Delta v/\Delta z)^2$',)
    (bz-vz2).resample(TIME='2D').mean().plot(ax=axx[1],color='purple',label='$\Delta b/\Delta z-(\Delta v/\Delta z)^2$',linewidth=2)
    axhline(0,color='r',zorder=10)
    for axi in axx:
        axi.set_xlabel('')
        axi.set_title('')
        axi.set_ylabel('x 10$^{-5}$ [1/s$^2$]')
        axi.legend()
        axi.set_ylim(-0.25,1)
    axx[0].set_title('CF5: '+str(d2)+'m - '+str(d1)+'m')
    xlim(datetime.datetime(2014,9,1),datetime.datetime(2016,8,1))
    savefig(figdir+'PV_'+str(d2)+'-'+str(d1)+'.png',bbox_inches='tight',dpi=300)
    xlim(datetime.datetime(2014,10,1),datetime.datetime(2015,6,1))
    savefig(figdir+'PV_'+str(d2)+'-'+str(d1)+'_zoom.png',bbox_inches='tight',dpi=300)
    xlim(datetime.datetime(2015,3,1),datetime.datetime(2015,6,1))
    savefig(figdir+'PV_'+str(d2)+'-'+str(d1)+'_zoom_more.png',bbox_inches='tight',dpi=300)

plot_PVcomps(250,500)
plot_PVcomps(500,750)
plot_PVcomps(750,1000)

def PV_in_space():
    figure(figsize=(20,5))
    for d1 in [250,500,750]:
        d2=d1+250
        bz=((mcat['PDEN'].sel(DEPTH=d2)-mcat['PDEN'].sel(DEPTH=d1))/(-mcat['z'].sel(DEPTH=d2)+mcat['z'].sel(DEPTH=d1))/100)*10**5
        vz2=((vel['ug'].sel(DEPTH=d1)-vel['ug'].sel(DEPTH=d2))**2/(-vel['z'].sel(DEPTH=d2)+vel['z'].sel(DEPTH=d1))**2)*10**5
        scatter(mcat.TIME.values,(mcat.z.sel(DEPTH=d2).values+mcat.z.sel(DEPTH=d1).values)/2, c=(bz-vz2).where(bz>0).values,vmin=-0.2,vmax=0.2,cmap=cm.RdBu_r)
    colorbar(label='$\Delta b/\Delta z-(\Delta v/\Delta z)^2$ x 10$^{-5}$ [1/s$^2$]')
    ylabel('depth [m]')
    xlim(datetime.datetime(2014,9,1),datetime.datetime(2016,8,1))
    savefig(figdir+'PV_in_space.png',bbox_inches='tight',dpi=300)
    xlim(datetime.datetime(2014,10,1),datetime.datetime(2015,6,1))
    savefig(figdir+'PV_in_space_zoom.png',bbox_inches='tight',dpi=300)
    xlim(datetime.datetime(2015,3,1),datetime.datetime(2015,6,1))
    savefig(figdir+'PV_in_space_zoom_more.png',bbox_inches='tight',dpi=300)

PV_in_space()

def PV_in_space():
    figure(figsize=(20,5))
    for d1 in [250,500,750]:
        d2=d1+250
        bz=((mcat['PDEN'].sel(DEPTH=d2)-mcat['PDEN'].sel(DEPTH=d1))/(-mcat['z'].sel(DEPTH=d2)+mcat['z'].sel(DEPTH=d1))/100)*10**5
        vz2=((vel['ug'].sel(DEPTH=d1)-vel['ug'].sel(DEPTH=d2))**2/(-vel['z'].sel(DEPTH=d2)+vel['z'].sel(DEPTH=d1))**2)*10**5
        mid_depth=(mcat.z.sel(DEPTH=d2).values+mcat.z.sel(DEPTH=d1).values)/2
        scatter(mcat.TIME.values,mid_depth, c=(bz-vz2).where(bz>0).values,vmin=-0.2,vmax=0.2,cmap=cm.RdBu_r)
    colorbar(label='$\Delta b/\Delta z-(\Delta v/\Delta z)^2$ [x 10$^{-5}$ 1/s$^2$]')
    ylabel('depth [m]')
    xlim(datetime.datetime(2014,9,1),datetime.datetime(2016,8,1))
    savefig(figdir+'PV_in_space.png',bbox_inches='tight',dpi=300)
    xlim(datetime.datetime(2014,10,1),datetime.datetime(2015,6,1))
    savefig(figdir+'PV_in_space_zoom.png',bbox_inches='tight',dpi=300)
    xlim(datetime.datetime(2015,3,1),datetime.datetime(2015,6,1))
    savefig(figdir+'PV_in_space_zoom_more.png',bbox_inches='tight',dpi=300)

    figure()
    plot(mid_depth-max(mid_depth),bz-vz2,'.')

PV_in_space()


for d1 in [250,500,750]:
    d2=d1+250
    bz=((mcat['PDEN'].sel(DEPTH=d2)-mcat['PDEN'].sel(DEPTH=d1))/(-mcat['z'].sel(DEPTH=d2)+mcat['z'].sel(DEPTH=d1))/100)*10**5
    vz2=((vel['ug'].sel(DEPTH=d1)-vel['ug'].sel(DEPTH=d2))**2/(-vel['z'].sel(DEPTH=d2)+vel['z'].sel(DEPTH=d1))**2)*10**5
    mcat['PV']=bz-vz2

    figure(figsize=(12,3))
    (mcat['PV']<0).resample(TIME='1D').sum().plot(label='all negative PV')
    (mcat['PV'].where(bz>0)<0).resample(TIME='1D').sum().plot(label='negative PV w/out inversions')
    ylabel('#')
    legend()
    title('Number of hours per day with negative PV\n CF5: '+str(d2)+'m-'+str(d1)+'m')
    savefig(figdir+'Count_negPV_'+str(d2)+'-'+str(d1)+'.png',bbox_inches='tight',dpi=300)


dat=xr.open_dataset('../data/aux_data/ERA_1804/tau_hflux_180411.nc')
era=dat.rename({'inss':'tauy','iews':'taux','longitude':'lon','latitude':'lat','time':'date'})
era['tau along']=(era['taux']*cos(theta)+era['tauy']*sin(theta))
era['lon']=era['lon']-360


def pltseries():
    figure(figsize=(10,4))
    era['tau along'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(label='At CF moorings',color='k')
    xlim([datetime.datetime(2014,8,15),datetime.datetime(2016,8,15)])
    legend()
    axhline(0,color='k')
    ylim(0,1)

pltseries()

def PV_blowdown():
    for d1 in [250,500,750]:
        d2=d1+250
        bz=((mcat['PDEN'].sel(DEPTH=d2)-mcat['PDEN'].sel(DEPTH=d1))/(-mcat['z'].sel(DEPTH=d2)+mcat['z'].sel(DEPTH=d1))/100)*10**5
        vz2=((vel['ug'].sel(DEPTH=d1)-vel['ug'].sel(DEPTH=d2))**2/(-vel['z'].sel(DEPTH=d2)+vel['z'].sel(DEPTH=d1))**2)*10**5
        mid_depth=(mcat.z.sel(DEPTH=d2).values+mcat.z.sel(DEPTH=d1).values)/2
        figure(figsize=(8,8))
        plot(mid_depth-max(mid_depth),bz-vz2,'.')

PV_blowdown()
