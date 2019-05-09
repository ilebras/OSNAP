from aux_funcs import *


dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1803extrap.pickle','rb'))

#0. fill in one day that is missing from M1
tprob=where(isnan(dat['salinity'].isel(distance=-1).mean(dim='depth'))==True)[0][0]
tprob
for vv in dat:
    if vv[0]!='d':
        dat[vv][-1,:,tprob]=(dat[vv].isel(date=tprob+1,distance=-1)+dat[vv].isel(date=tprob-1,distance=-1))/2

from eofs.xarray import Eof



for ii in range(3,8):
    if ii<4:
        skipnum=10
    else:
        skipnum=50
    solver_across = Eof(dat['across track velocity'][:,::skipnum,:].T[:,:,ii])
    solver_along = Eof(dat['along track velocity'][:,::skipnum,:].T[:,:,ii])

    figure(figsize=(10,5))
    subplot(121)
    plot(squeeze(solver_across.eofs(neofs=1)),-dat.depth[::skipnum],label='mode 1: '+str(int(solver_across.varianceFraction()[0].values*100))+'% of variance')
    plot(squeeze(solver_across.eofs()[1,:]),-dat.depth[::skipnum],label='mode 2: '+str(int(solver_across.varianceFraction()[1].values*100))+'% of variance')
    ylabel('depth [m]')
    xlabel('along-stream velocity [m/s]')
    xlim([-1,1])
    legend(loc=(0.2,1.01))
    subplot(122)
    plot(squeeze(solver_along.eofs(neofs=1)),-dat.depth[::skipnum],label='mode 1: '+str(int(solver_along.varianceFraction()[0].values*100))+'% of variance')
    plot(squeeze(solver_along.eofs()[1,:]),-dat.depth[::skipnum],label='mode 2: '+str(int(solver_along.varianceFraction()[1].values*100))+'% of variance')
    gca().set_yticklabels('')
    xlabel('across-stream velocity [m/s]')
    xlim([-1,1])
    legend(loc=(0.2,1.01))
    suptitle('CF'+str(ii+1))
    savefig('../figures/VertModes/eof_cf'+str(ii+1)+'.png')

    figure(figsize=(12,4))
    solver_across.pcs().sel(mode=0).plot(label='along-stream')
    solver_along.pcs().sel(mode=0).plot(label='across-stream')
    legend()
    title('CF'+str(ii+1)+' first mode time series')
    savefig('../figures/VertModes/eof_1804/eof_cf'+str(ii+1)+'_tseries.png')
    #
    # figure(figsize=(12,4))
    # solver.pcs().sel(mode=1).plot()
    # title('CF'+str(ii+1)+' second mode time series')

############################################################################
## Fit exponential to mean density and see what the solution looks like!
############################################################################

def alldenfig():
    figure(figsize=(4,6))
    plot(dat['potential density'].mean(dim='date').where(dat.distance>30).T,dat.depth)
    ylim([2000,0])
    xlabel('potential density [kg/m$^3$]')
    ylabel('depth [m]]')
    title('Mean density profiles at CF4-M1')
    savefig('../figures/VertModes/expden_1805/Mean_denprof.png')

alldenfig()

from scipy.optimize import curve_fit

def expfunc(x, a, b, c):
    return a * np.exp(-b * x) + c

dat.depth[25:]

def fitexp(denprof):
    fromd=25
    nonan_zero=~isnan(denprof)
    p0test,pcoc_test=curve_fit(expfunc,dat.depth[nonan_zero],denprof[nonan_zero])
    nonan=~isnan(denprof[fromd:])
    popt,pcoc=curve_fit(expfunc,dat.depth[fromd:][nonan],denprof[fromd:][nonan],p0=p0test)
    figure(figsize=(4,6))
    plot(denprof.T,dat.depth,label='original')
    plot(expfunc(dat.depth[nonan_zero],*popt).T,dat.depth[nonan_zero],label='fitted')
    legend()
    gca().invert_yaxis()
    xlabel('potential density [kg/m$^3$]')
    ylabel('depth [m]]')

    return popt


popt={}
for ii in range(3,8):
    popt[ii+1]=fitexp(dat['potential density'].mean(dim='date')[ii,:].values)
    title('CF'+str(ii+1)+', decay scale = '+str(int(1/popt[ii+1][1]))+'m')
    savefig('../figures/VertModes/expden_1805/Mdens_expfit_cf'+str(ii+1)+'_from50m.png')

from scipy import special

dat.depth[25]

for ii in range(3,8):
    if ii<4:
        skipnum=10
    else:
        skipnum=50
    solver_across = Eof(dat['across track velocity'][:,::skipnum,:].T[:,:,ii]-dat['across track velocity'][:,::skipnum,:].T[:,:,ii].mean(dim='date'))
    solver_along = Eof(dat['along track velocity'][:,::skipnum,:].T[:,:,ii]-dat['along track velocity'][:,::skipnum,:].T[:,:,ii].mean(dim='date'))
    denprof=dat['potential density'].mean(dim='date')[ii,:].values
    nonan=~isnan(denprof)
    popt,pcoc=curve_fit(expfunc,range(len(denprof[nonan])),denprof[nonan])


    figure(figsize=(10,5))
    subplot(121)
    plot(abs(squeeze(solver_across.eofs(neofs=1))),dat.depth[::skipnum],label='mode 1: '+str(int(solver_across.varianceFraction()[0].values*100))+'% of variance')
    plot(abs(dat['across track velocity'][:,::skipnum,:].mean(dim='date')[ii,:].T),dat.depth[::skipnum],label='Mean along-stream velocity')
    plot(abs(popt[0]*exp(-popt[1]/2*range(len(denprof[nonan])))*special.j1(2.4048*exp(-popt[1]/2*range(len(denprof[nonan]))))),dat.depth[nonan],label='Exponential density fit solution')
    ylabel('depth [m]')
    xlabel('along-stream velocity [m/s]')
    xlim([-0.1,1])
    legend(loc=(0.2,1.01))
    gca().invert_yaxis()
    subplot(122)
    plot(abs(squeeze(solver_along.eofs(neofs=1))),dat.depth[::skipnum],label='mode 1: '+str(int(solver_along.varianceFraction()[0].values*100))+'% of variance')
    plot(abs(dat['along track velocity'][:,::skipnum,:].mean(dim='date')[ii,:].T),dat.depth[::skipnum],label='Mean across-stream velocity')
    plot(abs(popt[0]*exp(-popt[1]/2*range(len(denprof[nonan])))*special.j1(2.4048*exp(-popt[1]/2*range(len(denprof[nonan]))))),dat.depth[nonan],label='Exponential density fit solution')
    gca().set_yticklabels('')
    xlabel('across-stream velocity [m/s]')
    xlim([-0.1,1])
    legend(loc=(0.2,1.01))
    suptitle('CF'+str(ii+1))
    gca().invert_yaxis()
    savefig('../figures/VertModes/expden_1805/eof_wfit_cf'+str(ii+1)+'.png',bbox_inches='tight')

dat

for ii in range(3,8):
    if ii<4:
        skipnum=10
    else:
        skipnum=50
    fromd=25
    solver_across = Eof(dat['across track velocity'][:,fromd::skipnum,:].T[:,:,ii]-dat['across track velocity'][:,fromd::skipnum,:].T[:,:,ii].mean(dim='date'))
    solver_along = Eof(dat['along track velocity'][:,fromd::skipnum,:].T[:,:,ii]-dat['along track velocity'][:,fromd::skipnum,:].T[:,:,ii].mean(dim='date'))
    denprof=dat['potential density'].mean(dim='date')[ii,:].values
    fromd=25
    nonan_zero=~isnan(denprof)
    p0test,pcoc_test=curve_fit(expfunc,dat.depth[nonan_zero],denprof[nonan_zero])
    nonan=~isnan(denprof[fromd:])
    popt,pcoc=curve_fit(expfunc,dat.depth[fromd:][nonan],denprof[fromd:][nonan],p0=p0test)


    figure(figsize=(10,5))
    subplot(121)
    plot(abs(squeeze(solver_across.eofs(neofs=1))),dat.depth[fromd::skipnum],label='mode 1: '+str(int(solver_across.varianceFraction()[0].values*100))+'% of variance')
    plot(abs(dat['across track velocity'][:,fromd::skipnum,:].mean(dim='date')[ii,:].T),dat.depth[fromd::skipnum],label='Mean along-stream velocity')
    plot(abs(popt[0]*exp(-popt[1]/2*dat.depth[nonan_zero])*special.j1(2.4048*exp(-popt[1]/2*dat.depth[nonan_zero]))),dat.depth[nonan_zero],label='Exponential density fit solution')
    ylabel('depth [m]')
    xlabel('along-stream velocity [m/s]')
    xlim([-0.1,1])
    legend(loc=(0.2,1.01))
    gca().invert_yaxis()
    subplot(122)
    plot(abs(squeeze(solver_along.eofs(neofs=1))),dat.depth[fromd::skipnum],label='mode 1: '+str(int(solver_along.varianceFraction()[0].values*100))+'% of variance')
    plot(abs(dat['along track velocity'][:,fromd::skipnum,:].mean(dim='date')[ii,:].T),dat.depth[fromd::skipnum],label='Mean across-stream velocity')
    plot(abs(popt[0]*exp(-popt[1]/2*dat.depth[nonan_zero])*special.j1(2.4048*exp(-popt[1]/2*dat.depth[nonan_zero]))),dat.depth[nonan_zero],label='Exponential density fit solution')
    gca().set_yticklabels('')
    xlabel('across-stream velocity [m/s]')
    xlim([-0.1,1])
    legend(loc=(0.2,1.01))
    suptitle('CF'+str(ii+1))
    gca().invert_yaxis()
    savefig('../figures/VertModes/expden_1805/eof_wfit_cf'+str(ii+1)+'from50.png',bbox_inches='tight')

### Try using shipboard data stratification!

CTD=pickle.load(open('../pickles/Shipboard/CTD_xarray_1805_finer.pickle','rb'))
deepdist=[int(dd/5)*5 for dd in distvec[3:]]
deepdist

CTD.occupation
popt
def fitexp_ship(denprof,ii):
    nonan=~isnan(denprof)
    popti,pcoc=curve_fit(expfunc,CTD.pressure[nonan],denprof[nonan],p0=popt[ii+4])
    return popti,nonan


popt_ship={}
popt_ship[2014]={}
popt_ship[2016]={}
for ii,dd in enumerate(deepdist):
    figure()
    plot(CTD.density.sel(distance=dd).sel(occupation='2014 (KN221)'),CTD.pressure,label='2014')
    plot(CTD.density.sel(distance=dd).sel(occupation='2016'),CTD.pressure,label='2016')
    plot(dat['potential density'].mean(dim='date').isel(distance=ii+3).T,dat.depth,label='mooring mean')
    popt_ship[2014][ii+4],nonan=fitexp_ship(CTD.density.sel(distance=dd).sel(occupation='2014 (KN221)').values,ii)
    popt_ship[2016][ii+4],nonan=fitexp_ship(CTD.density.sel(distance=dd).sel(occupation='2016').values,ii)
    plot(expfunc(CTD.pressure[nonan],*popt_ship[2014][ii+4]).T,CTD.pressure[nonan],label='fitted to 2014 profile')
    plot(expfunc(CTD.pressure[nonan],*popt_ship[2016][ii+4]).T,CTD.pressure[nonan],label='fitted to 2016 profile')
    xlabel('potential density [kg/m$^3$]')
    ylabel('depth [m]]')
    title('CF'+str(ii+4)+', mean decay scale = '+str(int(1/(popt_ship[2014][ii+4][1]+popt_ship[2016][ii+4][1])))+'m')
    gca().invert_yaxis()
    legend()
    savefig('../figures/VertModes/shipden_1805/Denprof_cf'+str(ii+3)+'.png',bbox_inches='tight')


popt_ship['mean']={}
for ii in popt_ship[2014]:
    popt_ship['mean'][ii]=(popt_ship[2014][ii]+popt_ship[2016][ii])/2
popt

for ii in range(3,8):
    if ii<4:
        skipnum=10
    else:
        skipnum=50
    solver_across = Eof(dat['across track velocity'][:,::skipnum,:].T[:,:,ii]-dat['across track velocity'][:,::skipnum,:].T[:,:,ii].mean(dim='date'))
    solver_along = Eof(dat['along track velocity'][:,::skipnum,:].T[:,:,ii]-dat['along track velocity'][:,::skipnum,:].T[:,:,ii].mean(dim='date'))
    nonan=~isnan(CTD.density.sel(distance=deepdist[ii-3]).sel(occupation='2016'))
    nonan_moor=~isnan(dat['potential density'].mean(dim='date').isel(distance=ii))

    figure(figsize=(10,5))
    subplot(121)
    plot(abs(squeeze(solver_across.eofs(neofs=1))),dat.depth[::skipnum],label='mode 1: '+str(int(solver_across.varianceFraction()[0].values*100))+'% of variance')
    plot(abs(dat['across track velocity'][:,::skipnum,:].mean(dim='date')[ii,:].T),dat.depth[::skipnum],label='Mean along-stream velocity')
    plot(abs(popt_ship['mean'][ii+1][0]*exp(-popt_ship['mean'][ii+1][1]/2*CTD.pressure[nonan])*special.j1(2.4048*exp(-popt_ship['mean'][ii+1][1]/2*CTD.pressure[nonan]))),CTD.pressure[nonan],label='Exponential density from ship')
    plot(abs(popt[ii+1][0]*exp(-popt[ii+1][1]/2*dat.depth[nonan_moor])*special.j1(2.4048*exp(-popt[ii+1][1]/2*dat.depth[nonan_moor]))),dat.depth[nonan_moor],label='Exponential density from mooring')
    ylabel('depth [m]')
    xlabel('along-stream velocity [m/s]')
    xlim([-0.1,1])
    legend(loc=(0.2,1.01))
    gca().invert_yaxis()
    subplot(122)
    plot(abs(squeeze(solver_along.eofs(neofs=1))),dat.depth[::skipnum],label='mode 1: '+str(int(solver_along.varianceFraction()[0].values*100))+'% of variance')
    plot(abs(dat['along track velocity'][:,::skipnum,:].mean(dim='date')[ii,:].T),dat.depth[::skipnum],label='Mean across-stream velocity')
    plot(abs(popt_ship['mean'][ii+1][0]*exp(-popt_ship['mean'][ii+1][1]/2*CTD.pressure[nonan])*special.j1(2.4048*exp(-popt_ship['mean'][ii+1][1]/2*CTD.pressure[nonan]))),CTD.pressure[nonan],label='Exponential density from ship')
    plot(abs(popt[ii+1][0]*exp(-popt[ii+1][1]/2*dat.depth[nonan_moor])*special.j1(2.4048*exp(-popt[ii+1][1]/2*dat.depth[nonan_moor]))),dat.depth[nonan_moor],label='Exponential density from mooring')
    gca().set_yticklabels('')
    xlabel('across-stream velocity [m/s]')
    xlim([-0.1,1])
    legend(loc=(0.2,1.01))
    suptitle('CF'+str(ii+1))
    gca().invert_yaxis()
    savefig('../figures/VertModes/shipden_1805/eof_shipden_wfit_cf'+str(ii+1)+'.png',bbox_inches='tight')

dat['potential density'].sel(date=slice('2015-1-1','2015-3-1')).mean(dim='date')

def stratfig():
    plot(dat['potential density'].sel(date=slice('2015-1-1','2015-3-1')).mean(dim='date').T,dat.depth,'b');
    plot(dat['potential density'].sel(date=slice('2015-8-1','2015-9-1')).mean(dim='date').T,dat.depth,'r');
    gca().invert_yaxis()

stratfig()
