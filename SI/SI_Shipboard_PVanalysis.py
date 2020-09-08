from aux_funcs_2020 import *


yrvec=[14,16,18]
dat={}
for ii in yrvec:
    dat[ii]=xr.open_dataset(datadir+'Shipboard/netcdf/CTD_SADCP_10m_20'+str(ii)+'.nc')

#######################################################################################
# Calculate all the components necessary for Ertel PV and interpolate to the same grid:
#######################################################################################
### Coriolis parameter:
for ii in dat:
    dat[ii]['f']=(['dist'],sw.f(dat[ii].lat.values))

#relative vorticity (horizontal gradient of velocity)
for ii in dat:
    dat[ii]=dat[ii].assign_coords(middist=dat[ii].dist.values[:-1]+diff(dat[ii].dist.values)/2)
    dat[ii]['zeta_I']=(['depth','middist'],diff(dat[ii].uacross.values,axis=1)/diff(dat[ii].dist.values*1e3))
    dat[ii]['zeta']=(['depth','dist'],dat[ii]['zeta_I'].interp(middist=dat[ii].dist.values))

# horizontal density gradient
for ii in dat:
    dat[ii]['drhodx_I']=(['depth','middist'],diff(dat[ii].pden.values,axis=1)/diff(dat[ii].dist.values*1e3))
    dat[ii]['drhodx']=(['depth','dist'],dat[ii]['drhodx_I'].interp(middist=dat[ii].dist.values))


#vertical density gradient and buoyancy frequency (N2)
g=9.8
for ii in dat:
    dat[ii]=dat[ii].assign_coords(middepth=dat[ii].depth.values[:-1]+diff(dat[ii].depth.values)/2)
    dat[ii]['drhodz_I']=(['middepth','dist'],(diff(dat[ii].pden.values,axis=0).T/diff(dat[ii].depth.values).T).T)
    #remove inversions and interpolate over them
    dat[ii]['drhodz_I']=dat[ii]['drhodz_I'].where(dat[ii]['drhodz_I']>0).interpolate_na(dim='middepth')
    dat[ii]['drhodz']=(['depth','dist'],dat[ii]['drhodz_I'].interp(middepth=dat[ii].depth.values))
    dat[ii]['N2']=(g*dat[ii]['drhodz']/(dat[ii]['pden']+1e3)).where(dat[ii].depth>50)

# vertical velocity gradient
for ii in dat:
    dat[ii]['dvdz_I']=(['middepth','dist'],(diff(dat[ii].uacross.values,axis=0).T/diff(dat[ii].depth.values).T).T)
    dat[ii]['dvdz']=(['depth','dist'],dat[ii]['dvdz_I'].interp(middepth=dat[ii].depth.values))

#geostrophic velocity
for ii in dat:
    dat[ii]['geoshear']=(-dat[ii]['drhodx']*g/((dat[ii]['f'].values)*(dat[ii]['pden']+1e3))).where(dat[ii].depth>50)
    dat[ii]['geovel_0bot']=(dat[ii]['geoshear'][::-1,:].T*diff(hstack((0,dat[ii].depth.values[:-1]+diff(dat[ii].depth.values)/2,2495)))[::-1]).cumsum(dim='depth').T
    dat[ii]['geovel']=dat[ii]['geovel_0bot']-dat[ii]['geovel_0bot'].sel(depth=100)+dat[ii]['uacross'].sel(depth=slice(50,150)).mean(dim='depth')
    dat[ii]['geovel']=dat[ii]['geovel'].where(~isnan(dat[ii]['pden']))

#geostrophic zeta
for ii in dat:
    dat[ii]['geozeta_I']=(['depth','middist'],diff(dat[ii].geovel.values,axis=1)/diff(dat[ii].dist.values*1e3))
    dat[ii]['geozeta']=(['depth','dist'],dat[ii]['geozeta_I'].interp(middist=dat[ii].dist.values))


#######################################################################################
# Subsample at approximate mooring sites
#######################################################################################

moordistvec=range(37,110,15)
moordepthvec=range(50,2000,150)

#mooring site based geostrophic zeta
for ii in dat:
    dat[ii]['moorvel']=dat[ii]['geovel'].interp(dist=moordistvec).interp(depth=moordepthvec).interp(depth=dat[ii].depth).interp(dist=dat[ii].dist)
    dat[ii]['moorzeta']=(['depth','middist'],diff(dat[ii].moorvel.values,axis=1)/diff(dat[ii].dist.values*1e3))
    # dat[ii]['moorzeta']=(['depth','dist'],dat[ii]['moorzeta_I'].interp(middist=dat[ii].dist.values))

for ii in dat:
    dat[ii]['moorN2']=(['depth','middist'],dat[ii]['N2'].interp(dist=moordistvec).interp(depth=moordepthvec).interp(depth=dat[ii].depth).interp(dist=dat[ii].middist))
    dat[ii]['moorshear']=(['depth','middist'],dat[ii]['geoshear'].interp(dist=moordistvec).interp(depth=moordepthvec).interp(depth=dat[ii].depth).interp(dist=dat[ii].middist))

#######################################################################################
# Calculate PV for different approximations:
#######################################################################################

#Ertel PV (note, leaving out the horizontal Coriolis component - not important at Gulf Stream and should be even less important up here)
for ii in dat:
    dat[ii]['EPV'] = (dat[ii]['f']+dat[ii]['zeta'])*dat[ii]['N2'] + dat[ii]['dvdz']*g*dat[ii]['drhodx']/(dat[ii]['pden']+1e3)
    dat[ii]['log(EPV)'] = log10((dat[ii]['f']+dat[ii]['zeta'])*dat[ii]['N2'] + dat[ii]['dvdz']*g*dat[ii]['drhodx']/(dat[ii]['pden']+1e3))
    dat[ii]['fxN2'] = dat[ii]['f']*dat[ii]['N2']
    dat[ii]['EPV_zetaxN2'] = dat[ii]['zeta']*dat[ii]['N2']
    dat[ii]['EPV_shearterm'] = dat[ii]['dvdz']*g*dat[ii]['drhodx']/(dat[ii]['pden']+1e3)

# Geostrophic Ertel PV
for ii in dat:
    dat[ii]['geoPV'] = (dat[ii]['f']+dat[ii]['geozeta'])*dat[ii]['N2'] - dat[ii]['f']*dat[ii]['geoshear']**2
    dat[ii]['log(geoPV)'] = log10((dat[ii]['f']+dat[ii]['geozeta'])*dat[ii]['N2'] - dat[ii]['f']*dat[ii]['geoshear']**2)
    dat[ii]['geoPV_zetaxN2'] = (dat[ii]['geozeta'])*dat[ii]['N2']
    dat[ii]['geoPV_shearterm'] =  - dat[ii]['f']*dat[ii]['geoshear']**2

#Mooring site based geostrophic Ertel PV calculation
for ii in dat:
    dat[ii]['moorPV'] = (['depth','middist'],((dat[ii]['f'].interp(dist=dat[ii].middist)+dat[ii]['moorzeta'])*dat[ii]['moorN2'] - dat[ii]['f'].interp(dist=dat[ii].middist)*dat[ii]['moorshear']**2).T)
    dat[ii]['log(moorPV)'] = (['depth','middist'],(log10((dat[ii]['f'].interp(dist=dat[ii].middist)+dat[ii]['moorzeta'])*dat[ii]['moorN2'] - dat[ii]['f'].interp(dist=dat[ii].middist)*dat[ii]['moorshear']**2)).T)
    dat[ii]['moorPV_fxN2'] =(['depth','middist'],(dat[ii]['f'].interp(dist=dat[ii].middist)*dat[ii]['moorN2']).T)
    dat[ii]['moorPV_zetaxN2'] =(['depth','middist'],( (dat[ii]['moorzeta'])*dat[ii]['moorN2']))
    dat[ii]['moorPV_shearterm'] = (['depth','middist'],(- dat[ii]['f'].interp(dist=dat[ii].middist)*dat[ii]['moorshear']**2).T)

## Differences between total PVs
for ii in dat:
    dat[ii]['EPV-geoPV']=dat[ii]['EPV']-dat[ii]['geoPV']
    dat[ii]['geoPV-moorPV']=dat[ii]['geoPV'].interp(dist=dat[ii].middist)-dat[ii]['moorPV']
    dat[ii]['(EPV-geoPV)_zetaxN2']=dat[ii]['EPV_zetaxN2']-dat[ii]['geoPV_zetaxN2']
    dat[ii]['(geoPV-moorPV)_zetaxN2']=dat[ii]['geoPV_zetaxN2'].interp(dist=dat[ii].middist)-dat[ii]['moorPV_zetaxN2']
    dat[ii]['zetadiff']=dat[ii]['geozeta'].interp(dist=dat[ii].middist)-dat[ii]['moorzeta']


#######################################################################################
################################### PLOT ##############################################
#######################################################################################

def make_uni(varit,cmapo,vmini,vmaxi):
    uni[varit]={}
    uni[varit]['cmap']=cmapo
    uni[varit]['vmin']=vmini
    uni[varit]['vmax']=vmaxi

    return uni

uni=make_uni('uacross',cm.RdBu_r,-0.6,0.6)
uni=make_uni('log',cm.rainbow,-11.5,-8.5)
uni=make_uni('PV',cm.YlGnBu,0,1e-9)
uni=make_uni('symmterm',cm.RdBu_r,-1e-9,1e-9)
uni=make_uni('zeta',cm.RdBu_r,-5e-5,5e-5)
uni=make_uni('zetadiff',cm.RdBu_r,-5e-5,5e-5)

def plot_3_panel(var):
    if 'log' in var:
        varit='log'
    elif (var=='moorPV') | (var=='EPV') | (var=='geoPV') | ('fxN2' in var):
        varit='PV'
    elif ('zetaxN2' in var) | ('shearterm' in var)| ('-' in var):
        varit='symmterm'
    elif ('moorzeta' in var) | ('geozeta' in var):
        varit='zeta'
    else:
        varit=var
    f,axx=subplots(3,1,sharex=True,sharey=True,figsize=(12,15))
    if ('moor' in var) | ('zetadiff'==var):
        xi='middist'
    else:
        xi='dist'
    for ii,yy in enumerate(yrvec):
        dat[yy][var].plot(ax=axx[ii],x=xi,y='depth',cmap=uni[varit]['cmap'],vmin=uni[varit]['vmin'],vmax=uni[varit]['vmax'])
        axx[ii].set_ylabel('depth [m]')
        axx[ii].set_xlabel('')
        axx[ii].set_title('20'+str(yy))
        for xx in moordistvec:
            axx[ii].axvline(xx,color='grey')
        for yy in moordepthvec:
            axx[ii].axhline(yy,color='grey')

    axx[2].set_xlabel('distance [km]')
    xlim(36,107)
    if ('EPV' in var) | (var=='uacross'):
        ylim(450,0)
    else:
        ylim(2000,0)
    savefig(figdir+'SI/Shipboard/Shipboard_3panel_'+var+'.png',bbox_inches='tight')

yrvec

def lineplotcomp(var,dpsel):
    if var=='zeta':
        extravec=['','geo','moor']
    else:
        extravec=['EPV_','geoPV_','moorPV_']
    f,axx=subplots(3,1,sharex=True,sharey=True,figsize=(8,10))
    for ii,yy in enumerate(yrvec):
        for extra in extravec:
            dat[yy][extra+var].sel(depth=dpsel).plot(ax=axx[ii],label=extra+var,linewidth=2)
        multi=3
        (dat[yy][extravec[-1]+var].sel(depth=dpsel)*multi).plot(ax=axx[ii],label=extra+var+' x'+str(multi),linewidth=2)
        axx[ii].set_ylabel(var)
        axx[ii].set_xlabel('')
        axx[ii].set_title('20'+str(yy)+' ('+str(dpsel)+'m)')
        axx[ii].set_xlim(40,95)
    axx[2].set_xlabel('distance [km]')
    legend()
    savefig(figdir+'SI/Shipboard/Shipboard_linecomp_'+var+'_'+str(dpsel)+'m.png',bbox_inches='tight')



lineplotcomp('zeta',50)

lineplotcomp('zeta',750)
lineplotcomp('shearterm',250)

lineplotcomp('zetaxN2',750)


for var in ['moorzeta','geozeta','zetadiff']:
    plot_3_panel(var)

for var in ['(EPV-geoPV)_zetaxN2','(geoPV-moorPV)_zetaxN2']:
    plot_3_panel(var)

for var in ['EPV-geoPV','geoPV-moorPV']:
    plot_3_panel(var)

for var in ['moorPV','log(moorPV)','moorPV_fxN2','moorPV_zetaxN2','moorPV_shearterm']:
    plot_3_panel(var)

for var in ['geoPV','log(geoPV)','geoPV_zetaxN2','geoPV_shearterm']:
    plot_3_panel(var)

for var in ['EPV','log(EPV)','fxN2','EPV_zetaxN2','EPV_shearterm']:
    plot_3_panel(var)
