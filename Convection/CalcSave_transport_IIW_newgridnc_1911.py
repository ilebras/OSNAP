#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

# versname='1810JHIL'
# daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_'+versname+'.pickle','rb'))
#
# daily['across track velocity']=-1*daily['across track velocity']
#
# from firstfuncs_1618 import *

daily=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid.nc')
daily['across track velocity']=-1*daily['across track velocity']
daily=daily.where(daily['across track velocity']!=0)

daily=daily.sel(date=slice('2014-8-15','2016-11-1')) #just to avoid some zeros at the beginning which mess up filtering.
#   PLUS HERE IM RESTRICTING DATES TO FIRST DEPLOYMENT(ISH)

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

onesxr=daily.salinity/daily.salinity


srefb=34.9
sep=9


daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3



egic={}
egic['trans']=daily['xport over 27.8'][sep:,:,:].sum('depth').sum('distance')
egic['area']=(onesxr.where(daily['potential density']<27.8)[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport over 27.8'][sep:,:-1,:]*daily['salinity'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport over 27.8'][sep:,:-1,:]*daily['temperature'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport over 27.8'][sep:,:-1,:]*daily['potential density'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']

for ii in range(0,200,20):
    figure()
    (daily['potential density'][sep:,:-1,:].where((daily['potential density']<d2)&(daily['potential density']>=d1)).isel(date=ii).T).plot()
    (daily['potential density'][sep:,:-1,:].where((daily['potential density']<d3)&(daily['potential density']>=d2)).isel(date=ii).T).plot()

uIIW={}
uIIW['trans']=daily.xport.where((daily['potential density']<d2)&(daily['potential density']>=d1)).sum('distance').sum('depth')
uIIW['area']=(onesxr.where((daily['potential density']<d2)&(daily['potential density']>=d1))[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
uIIW['trans cf5+']=daily.xport.where(daily.distance>=45).where((daily['potential density']<d2)&(daily['potential density']>=d1)).sum('distance').sum('depth')
uIIW['meanvel']=uIIW['trans']/uIIW['area']

dIIW={}
dIIW['trans']=daily.xport.where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')
dIIW['area']=(onesxr.where((daily['potential density']<d3)&(daily['potential density']>=d2))[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
dIIW['meanvel']=dIIW['trans']/dIIW['area']
dIIW['trans cf5+']=daily.xport.where(daily.distance>=45).where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')

IIW={}
IIW['trans']=daily.xport.where((daily['potential density']<d3)&(daily['potential density']>=d1)).sum('distance').sum('depth')

lessthan={}
lessthan['trans']=daily.xport.where(daily['potential density']<d1).sum('distance').sum('depth')

morethan={}
morethan['trans']=daily.xport.where((daily['potential density']>=d3)&(daily['potential density']<27.8)).sum('distance').sum('depth')


pickle.dump([uIIW,dIIW,IIW,egic,lessthan,morethan],open(datadir+'OSNAP2016recovery/pickles/convection_xport/IIWtrans_direct.pickle','wb'),protocol=2)




dat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_rotnewfieldnames.nc')
dat=dat.rename({'potential density':'pden'})
###############################################################################################################################
############################################## Calc and add PV ####################################################
###############################################################################################################################


PVmat=NaN*ones((len(dat.distance),len(dat.depth)-1,len(dat.date)))
for ii in range(len(dat.distance)):
    PVmat[ii,:,:]=sw.f(dat.lat[ii])*dat['potential density'][ii,:,:].diff(dim='depth')/mean(diff(dat.depth))/(dat['potential density'].values[ii,:-1,:]+dat['potential density'][ii,:,:].diff(dim='depth').values/2+1e3)

dat=dat.assign_coords(mid_depth=(dat.depth.values[:-1]+diff(dat.depth.values)/2))
PVmat[PVmat<=1e-15]=1e-15
dat['PV']=(('distance','mid_depth','date'),PVmat)
dat.PV.sel(distance=distvec[4]).plot()
