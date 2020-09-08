from firstfuncs_1618 import *

dat=xr.open_dataset(datadir+'aux_data/FramStrait/Fram_Strait_gridded_monthly_mean_V_S_2002-2015-v1.0.nc')

dat.PRES.sum(dim='Z').plot()

[lonlen,zlen,tlen]=shape(dat.PSAL)

zlen

lonmat=tile(dat.LONGITUDE.values,(zlen,1)).T

def quicksec(var):
    contourf(lonmat,dat.PRES.mean(dim='TIME'),dat[var].mean(dim='TIME'));
    colorbar()
    gca().invert_yaxis()
    xlabel('longitude')
    ylabel('presure [db]')
    title(var)

dat

quicksec('PSAL')

quicksec('V_VEL')


xport_WSC=xr.open_dataset(datadir+'aux_data/FramStrait/OS_1_FRAM_WSC_transport_1997_2012_D.NC')

xport_WSC.summary

xport_EGC_0=xr.open_dataset(datadir+'aux_data/FramStrait/Fram_Strait_freshwater_transport_1997-2002-v1.0.nc')
xport_EGC_1=xr.open_dataset(datadir+'aux_data/FramStrait/Fram_Strait_freshwater_transport_2002-2015-v1.0.nc')
xport_EGC_2=xr.open_dataset(datadir+'aux_data/FramStrait/Fram_Strait_freshwater_transport_2003-2015-v1.0.nc')

xport_EGC_0.FWT.plot()
xport_EGC_1.FWT.plot()
xport_EGC_2.FWT.plot()

xr.open_dataset(datadir+'aux_data/FramStrait/Fram_Strait_freshwater_transport_2002-2015-v1.0.nc')


##############################################################################
##############################################################################
#### Get my own volume transport and salinity from the monthly gridded data
#### Also compare FWT to above to be able to contextualize
##############################################################################
#############################################################################

# make an area field:
dat.LONGITUDE.values[:-1]+
dat.LATITUDE.mean().values*ones(len(dat.midlon))
dat=dat.assign_coords(midlon=(hstack((dat.LONGITUDE.values[0]-0.125,hstack((dat.LONGITUDE.values[:-1]+diff(dat.LONGITUDE.values)/2,dat.LONGITUDE.values[-1]+0.125))))))
dat['distance']=(['LON'],sw.dist(dat.LATITUDE.mean().values*ones(len(dat.midlon)),dat.midlon)[0])
dat['depth']=(['LON','Z','TIME'],-gsw.z_from_p(dat.PRES,dat.LATITUDE.mean().values))

0*dat.depth[:,0,:]

dat.depth[:,:-1,:]+dat.depth.diff(dim='Z')/2

dat.depth[:,-1,:]+dat.depth.diff(dim='Z')[:,-1,:]/2

dat['depthdiff']=xr.concat([0*dat.depth[:,0,:],dat.depth[:,:-1,:]+dat.depth.diff(dim='Z')/2,dat.depth[:,-1,:]+dat.depth.diff(dim='Z')[:,-1,:]/2],dim='Z').diff(dim='Z')



dat['area']=dat.depthdiff*dat.distance/1e3


dat['TRANS']=(dat.area*dat.V_VEL/100).sum(dim='Z').sum(dim='LON')

dat['S'] = (dat.area*dat.PSAL*dat.V_VEL/100).sum(dim='Z').sum(dim='LON')/(dat.area*dat.V_VEL/100).sum(dim='Z').sum(dim='LON')

sref=34.9
dat['FWT']=(dat.area*(sref-dat.PSAL)*dat.V_VEL/100).sum(dim='Z').sum(dim='LON')/sref*1e3

xport_EGC_0.summary
xport_EGC_1.summary
xport_EGC_2.summary
dat.summary
def fwcomp():
    figure(figsize=(12,3))
    xport_EGC_0.FWT.plot(label='1997-2002 (79N, 6.5W-1W)')
    xport_EGC_1.FWT.plot(label='2002-2015 (78.83N, 6.5W-1W)')
    xport_EGC_2.FWT.plot(label='2003-2015(78.83N, 8W-1W)')
    dat.FWT.plot(label='calculated from monthly fields')
    legend()

fwcomp()



figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/FramStrait/'


def plot_FSprops():
    figure(figsize=(12,7))
    subplot(211)
    dat.TRANS.plot(label='mean = -9.2 Sv')
    xlabel('')
    ylabel('Transport [Sv]')
    title('Fram Strait Polar Water (8W-1W)')
    legend()
    subplot(212)
    dat.S.plot(label='mean = 34.63')
    ylabel('Transport weighted salinity, S')
    xlabel('')
    legend()
    savefig(figdir+'FramStrait_SU_tseries.png',bbox_inches='tight')



plot_FSprops()



dat.S.plot(figsize=(12,3))
