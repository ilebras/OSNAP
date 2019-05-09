#### Find mean flow direction

#from aux_funcs import *

from map_funcs import *

dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1801.pickle','rb'))

u=dat['across track velocity']*cos(theta)-dat['along track velocity']*sin(theta)

v=dat['across track velocity']*sin(theta)+dat['along track velocity']*cos(theta)

#plot mean over top 150m
umean=u[:,:76,:].mean(dim='date').mean(dim='depth')
vmean=v[:,:76,:].mean(dim='date').mean(dim='depth')

sqrt(vmean**2+umean**2)

def quivermap():
    map=makeMap('CF')
    map.plot(CFlon,CFlat,'ko', markersize=10,zorder=98,latlon=True)
    qq=map.quiver(CFlon,CFlat,umean,vmean,latlon=True,linewidth=0.5,scale=2,zorder=100)#,scale_units='xy')
    title('Mean velocity top 150m',fontsize=25)
    # savefig('../figures/flowdirection/Map_wvelarrows_150_nogrid.pdf',bbox_inches='tight')
    # savefig('../figures/flowdirection/Map_wvelarrows_150_nogrid.png',bbox_inches='tight')

quivermap()

####################################################################################
##### Make quiver map of shipboard ADCP top 150
####################################################################################

adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')
adcp_16=io.loadmat(datadir+'Shipboard/ar07_2016/ar07_2016_vm_adcp_ilebras.mat')

def dpthmn(vmdic,choosed,var):
    meanout=zeros(len(vmdic['lon'][0]))
    for ii in range(len(meanout)):
        dlim=where(vmdic['depth'][:,ii]>choosed)[0][0]
        meanout[ii]=sum(vmdic[var][:dlim,ii]*diff(hstack((-median(diff(vmdic['depth'][:dlim,ii]))/2,vmdic['depth'][:dlim,ii]))))/vmdic['depth'][dlim,ii]
    return meanout

def quivermap_vm(vmdic,yr):
    map=makeMap('CF')
    # map.plot(CFlon,CFlat,'ko', markersize=10,zorder=98,latlon=True)
    qq=map.quiver(vmdic['lon'],vmdic['lat'],dpthmn(vmdic,160,'u'),dpthmn(vmdic,160,'v'),latlon=True,linewidth=0.5,scale=2,zorder=100)#,scale_units='xy')
    title('Mean velocity top 150m \n 20'+yr+' vessel mounted ADCP',fontsize=25)
    savefig('../figures/flowdirection/Map_wvelarrows_150_nogrid_VMADCP'+yr+'.pdf',bbox_inches='tight')
    savefig('../figures/flowdirection/Map_wvelarrows_150_nogrid_VMADCP'+yr+'.png',bbox_inches='tight')

quivermap_vm(adcp_14,'14')

quivermap_vm(adcp_16,'16')


# This is a vestige of when I was trying to figure out the best projection direction
# xdist=sw.dist([CFlat[-2],CFlat[-2]],[CFlon[0],CFlon[-2]])[0][0]
# ydist=sw.dist([CFlat[0],CFlat[-2]],[CFlon[-2],CFlon[-2]])[0][0]
# prevtheta=abs(arctan((xdist)/(ydist)))
#
# figure()
# quiver(CFlon,CFlat,umean,vmean,scale=1,units='xy')
# plot(CFlon,CFlat,'g-')
# ii=0
# plot([CFlon[ii],CFlon[ii]-0.5*cos(theta)],[CFlat[ii],CFlat[ii]-0.5*sin(theta)],'r-',label='along flow direction (new projection)')
# plot([CFlon[ii],CFlon[ii]-0.5*cos(prevtheta)],[CFlat[ii],CFlat[ii]-0.5*sin(prevtheta)],'g-',label='across track direction')
# [plot([CFlon[ii],CFlon[ii]-0.5*cos(theta)],[CFlat[ii],CFlat[ii]-0.5*sin(theta)],'r-') for ii in range(8)];
# [plot([CFlon[ii],CFlon[ii]-0.5*cos(prevtheta)],[CFlat[ii],CFlat[ii]-0.5*sin(prevtheta)],'g-') for ii in range(8)];
# axis('equal')
# legend()
# savefig('../figures/flowdirection/newprojdir.pdf',bbox_inches='tight')
# savefig('../figures/flowdirection/newprojdir.png',bbox_inches='tight')
