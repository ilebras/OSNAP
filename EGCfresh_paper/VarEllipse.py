from aux_funcs import *


dat=pickle.load(open('../../pickles/xarray/CF_xarray_notid_1803extrap.pickle','rb'))

#0. fill in one day that is missing from M1
tprob=where(isnan(dat['salinity'].isel(distance=-1).mean(dim='depth'))==True)[0][0]
tprob
for vv in dat:
    if vv[0]!='d':
        dat[vv][-1,:,tprob]=(dat[vv].isel(date=tprob+1,distance=-1)+dat[vv].isel(date=tprob-1,distance=-1))/2

## Get some variance ellipses for k and l relationship and along-stream justification

#transform back to u and v
u=dat['across track velocity']*cos(theta)-dat['along track velocity']*sin(theta)
v=dat['across track velocity']*sin(theta)+dat['along track velocity']*cos(theta)

#start with depth mean vels
umean=u.where(u.depth<4e3).mean(dim='depth')
vmean=v.where(u.depth<4e3).mean(dim='depth')

umeana=umean-umean.mean(dim='date')
vmeana=vmean-vmean.mean(dim='date')

savename='full'
klrat={}
for ii in range(8):
    covu=cov(umeana[ii,:]**2)
    covv=cov(vmeana[ii,:]**2)
    covuv=cov(umeana[ii,:]*vmeana[ii,:])

    varmax=0.5*(covu+covv+ sqrt((covu-covv)**2+4*covuv**2))
    varmax
    varmin=covu**2+covv**2-varmax
    varmin

    import xarray as xray
    uvmat=xray.concat((umeana,vmeana),dim='uv')

    from eofs.xarray import Eof
    solver= Eof(uvmat.isel(distance=ii).T)

    evec1x=solver.eofs(neofs=1)[0][0].values
    evec1y=solver.eofs(neofs=1)[0][1].values

    maxtheta_rad=arctan(evec1y/evec1x)
    maxtheta=arctan(evec1y/evec1x)*180/pi
    maxtheta

    eig1=solver.eigenvalues()[0].values
    eig2=solver.eigenvalues()[1].values
    maxvar=2*sqrt(5.991*eig1)
    minvar=2*sqrt(5.991*eig2)


    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    hist2d(umean[ii,:],vmean[ii,:],30);
    colorbar()
    e=Ellipse(xy=hstack((umean.mean(dim='date')[ii],vmean.mean(dim='date')[ii])),
    width=maxvar,height=minvar,angle=maxtheta,
    edgecolor='w',facecolor='none',lw=3)
    ax.add_artist(e)
    mabsrad=abs(maxtheta_rad)
    klrat[ii]=(eig1*sin(mabsrad)+eig2*cos(mabsrad))/(eig1*cos(mabsrad)+eig2*sin(mabsrad))
    title('CF'+str(ii+1)+', '+savename+'\n k='+str('{:5.2f}'.format(klrat[ii]))+'x l, angle = '+str('{:5.2f}'.format(maxtheta)))
    savefig('../../figures/VertModes/varellipse/varellipse_'+savename+'_cf'+str(ii+1)+'.png')

#Use a constant klrat for simplicity, but have all available if want to go that route
mean([klrat[kk] for kk in range(3,8)])
klrat_con=1.6

help(sw.bfrq)
