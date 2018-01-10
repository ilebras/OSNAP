#### Find mean flow direction

#from aux_funcs import *

from map_funcs import *

dat2=pickle.load(open('../pickles/CF_xarray_notid_ptcorr.pickle','rb'))
trackv=dat2['across track velocity'].mean(dim='date').mean(dim='depth')
tracku=dat2['along track velocity'].mean(dim='date').mean(dim='depth')

dat=pickle.load(open('../pickles/CF_xarray_notid_newtheta.pickle','rb'))
flowv=dat['across track velocity'].mean(dim='date').mean(dim='depth')
flowu=dat['along track velocity'].mean(dim='date').mean(dim='depth')

#just checking that new flow rotation is really more aligned with flow. looks good!
quiver(flowu,flowv,color='r')
quiver(tracku,trackv)
grid('on')

u=dat['across track velocity']*cos(theta)-dat['along track velocity']*sin(theta)

v=dat['across track velocity']*sin(theta)+dat['along track velocity']*cos(theta)

umean=u.mean(dim='date').mean(dim='depth')
vmean=v.mean(dim='date').mean(dim='depth')

def quivermap():
    map=makeMap()
    map.quiver(CFlon,CFlat,umean,vmean,latlon=True,linewidth=0.5,scale=2,zorder=100)
    title('Time and depth averaged flow at Cape Farewell moorings')
    savefig('../figures/flowdirection/Map_wvelarrows.pdf',bbox_inches='tight')
    savefig('../figures/flowdirection/Map_wvelarrows.png',bbox_inches='tight')
quivermap()

xdist=sw.dist([CFlat[-2],CFlat[-2]],[CFlon[0],CFlon[-2]])[0][0]
ydist=sw.dist([CFlat[0],CFlat[-2]],[CFlon[-2],CFlon[-2]])[0][0]
prevtheta=abs(arctan((xdist)/(ydist)))

figure()
quiver(CFlon,CFlat,umean,vmean,scale=1,units='xy')
plot(CFlon,CFlat,'g-')
ii=0
plot([CFlon[ii],CFlon[ii]-0.5*cos(theta)],[CFlat[ii],CFlat[ii]-0.5*sin(theta)],'r-',label='along flow direction (new projection)')
plot([CFlon[ii],CFlon[ii]-0.5*cos(prevtheta)],[CFlat[ii],CFlat[ii]-0.5*sin(prevtheta)],'g-',label='across track direction')
[plot([CFlon[ii],CFlon[ii]-0.5*cos(theta)],[CFlat[ii],CFlat[ii]-0.5*sin(theta)],'r-') for ii in range(8)];
[plot([CFlon[ii],CFlon[ii]-0.5*cos(prevtheta)],[CFlat[ii],CFlat[ii]-0.5*sin(prevtheta)],'g-') for ii in range(8)];
axis('equal')
legend()
savefig('../figures/flowdirection/newprojdir.pdf',bbox_inches='tight')
savefig('../figures/flowdirection/newprojdir.png',bbox_inches='tight')
