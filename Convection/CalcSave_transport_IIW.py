#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

versname='1810JHIL'
daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_'+versname+'.pickle','rb'))

daily['across track velocity']=-1*daily['across track velocity']

#################################################################################
###################### (Freshwater) Transport #################################
#################################################################################

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

onesxr=daily.salinity/daily.salinity


srefb=34.9
sep=9

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
uIIW['trans cf5+']=daily.xport.where(daily.distance>=45).where((daily['potential density']<d2)&(daily['potential density']>=d1)).sum('distance').sum('depth')

dIIW={}
dIIW['trans']=daily.xport.where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')
dIIW['trans cf5+']=daily.xport.where(daily.distance>=45).where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')

IIW={}
IIW['trans']=daily.xport.where((daily['potential density']<d3)&(daily['potential density']>=d1)).sum('distance').sum('depth')


pickle.dump([uIIW,dIIW,IIW,egic],open(datadir+'OSNAP2016recovery/pickles/convection_xport/IIWtrans_direct.pickle','wb'),protocol=2)


"XXXXXXXXXXXXXXX
# def implement_fitsin(field):
#     figure(figsize=(9,9))
#     # fit a sin to the filtered transport
#     [trans_fit,trans_std,trans_period]=fitsin(tvec,field['trans filt'],mean(field['trans filt']),30,std(field['trans filt']),365.25)
#     subplot(211)
#     plot(daily.date,trans_fit)
#     plot(daily.date,field['trans filt'])
#     # print('trans period: ',trans_period)
#     print('trans corr')
#     print(corrcoef(trans_fit,field['trans filt'])[0,1])
#     print('trans min')
#     print(min(trans_fit))
#     print(daily.date[argmin(trans_fit)].values)
#     print('trans max')
#     print(max(trans_fit))
#     print(daily.date[argmax(trans_fit)].values)
#     print('trans percent increase')
#     print((max(trans_fit)-min(trans_fit))/min(trans_fit))
#     # fit a sin to the filtered freshwater
#     [fresh_fit,fresh_std,fresh_period]=fitsin(tvec,field['fresh filt'],mean(field['fresh filt']),30,std(field['fresh filt']),365.25)
#     subplot(212)
#     plot(daily.date,fresh_fit)
#     plot(daily.date,field['fresh filt'])
#     # print('fresh period: ',fresh_period)
#     print('fresh corr')
#     print(corrcoef(fresh_fit,field['fresh filt'])[0,1])
#     print('fresh min')
#     print(min(fresh_fit))
#     print(daily.date[argmin(fresh_fit)].values)
#     print('fresh max')
#     print(max(fresh_fit))
#     print(daily.date[argmax(fresh_fit)].values)
#     print('fresh percent increase')
#     print((max(fresh_fit)-min(fresh_fit))/min(fresh_fit))
#     return trans_fit,fresh_fit
#
# cc['trans seas'],cc['fresh seas']=implement_fitsin(cc)
# egic['trans seas'],egic['fresh seas']=implement_fitsin(egic)
# eg['trans seas'],eg['fresh seas']=implement_fitsin(eg)
# ic['trans seas'],ic['fresh seas']=implement_fitsin(ic)
