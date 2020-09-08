from aux_funcs import *

#2018: use section 2 which is just upstream as it is in order and closer to synoptic:
sec2_2018=hstack((range(27,49),range(51,55),range(56,58)))
#2020: shorter onshore section
#note section 5, downstream is more complete, but has no processed sadcp
# section EGC 1 is also more complete, but CTD only goes to 500m
#going to leave 2020 out of it for now anyway.
# sec4_2020=range(97,113)

sadcp={}
sadcp[14]=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_SADCP_14_10m.nc')
sadcp[16]=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_SADCP_16_10m.nc')
sadcp[18]=xr.open_dataset(datadir+'Shipboard/netcdf/AllSADCP_18.nc').sel(sta=sec2_2018+1)

ctd={}
ctd[14]=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_CTD_14.nc')
ctd[16]=xr.open_dataset(datadir+'Shipboard/netcdf/CFsec_CTD_16.nc')
ctd[18]=xr.open_dataset(datadir+'Shipboard/netcdf/AllCTD18_2m_Leahproc.nc').sel(sta=sec2_2018)


for dd in [14,16]:
    sadcp[dd]=sadcp[dd].where(sadcp[dd].dist<150)
    ctd[dd]=ctd[dd].where(ctd[dd].dist<150)

yrvec=[14,16,18]
def plot_map():
    for ii,dd in enumerate(yrvec):
        plot(ctd[dd].lon,ctd[dd].lat,'o',label='20'+str(dd),color='C'+str(ii))
        plot(sadcp[dd].lon,sadcp[dd].lat,'*',color='C'+str(ii),mec='k')
    legend()

plot_map()

### make the coordinates more regular
sadcp[18]['depth']=(['prs_u'],-gsw.z_from_p(sadcp[18]['prs_u'],mean(sadcp[18]['lat'])))
sadcp[18]=sadcp[18].swap_dims({'prs_u':'depth'})
# sadcp[18]=sadcp[18].interp(depth=range(10,450,10))
sadcp[18]['uacross']=sadcp[18]['v']*sin(theta)+sadcp[18]['u']*cos(theta) #note, using the same direction as at CF line


sadcp[18]['dist']=(['sta'],ctd[18]['dist'])
sadcp[18]=sadcp[18].swap_dims({'sta':'dist'})

ctd[18]['depth']=(['prs'],-gsw.z_from_p(ctd[18]['prs'],mean(ctd[18]['lat'])))
ctd[18]=ctd[18].swap_dims({'prs':'depth'}).swap_dims({'sta':'dist'})

ctd[18]=ctd[18].sortby('dist').sel(dist=slice(40,120))
sadcp[18]=sadcp[18].sortby('dist').sel(dist=slice(40,120))

sadcp[18]=sadcp[18].drop('prs_u').drop('u').drop('v').drop('sta')

sadcp[18]['uacross']=sadcp[18]['uacross'].interpolate_na(dim='depth')

plot_map()
ctd[18]=ctd[18].rename_vars({'den':'pden'})

for ii in yrvec:
    ctd[ii]=ctd[ii].interp(depth=range(10,2500,10))
    sadcp[ii]=sadcp[ii].interp(depth=range(10,2500,10))
    sadcp[ii]['dist']=ctd[ii]['dist']
    sadcp[ii]['lon']=ctd[ii]['lon']
    sadcp[ii]['lat']=ctd[ii]['lat']


dat={}
for ii in yrvec:
    dat[ii]=xr.merge([ctd[ii],sadcp[ii]])

for ii in yrvec:
    for var in ['uacross','sal','tmp','pden']:
        figure()
        dat[ii][var].plot()
        gca().invert_yaxis()
        title(str(ii))


###### write it out!
yrvec
for ii in yrvec:
    dat[ii].to_netcdf(datadir+'Shipboard/netcdf/CTD_SADCP_10m_20'+str(ii)+'.nc','w')
