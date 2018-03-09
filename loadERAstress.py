from aux_funcs import *

nclist=sort(glob.glob('../data/aux_data/ERA_stress/*.nc'))

#create arrays that are time x lat x lon for taux and tauy

dat={}

dat['2014']=Dataset(nclist[0])
dat['2015']=Dataset(nclist[1])
dat['2016']=Dataset(nclist[2])

lon=dat['2014']['longitude'][:]
lat=dat['2014']['latitude'][:]

taux=vstack((dat['2014'].variables['iews'][:],dat['2015'].variables['iews'][:],dat['2016'].variables['iews'][:]))
tauy=vstack((dat['2014'].variables['inss'][:],dat['2015'].variables['inss'][:],dat['2016'].variables['inss'][:]))

from dateutil.relativedelta import relativedelta

datevec=[datetime.datetime(2014,1,15)+relativedelta(months=ii) for ii in range(36)]

era=xr.Dataset({'taux': (['date','lat','lon'],  taux),
                'tauy': (['date','lat','lon'],  tauy),},
                coords={'date': datevec,
                        'lat': lat,
                        'lon': lon-360})


from mpl_toolkits.basemap import Basemap

lat_start=50
lat_end  =70
lon_start=-50

lon_end  =-20
lon_0= - (abs(lon_end)+abs(lon_start))/2.0

map = Basemap(llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='l',projection='cyl',
                lat_1=lat_start)

lonmat,latmat=meshgrid(era.lon,era.lat)
x, y = map(lonmat,latmat)

# This did not work very well...
# dlon=vstack((gsw.distance(lonmat,latmat)[:,0],gsw.distance(lonmat,latmat).T)).T
# dlat=vstack((gsw.distance(lonmat.T,latmat.T)[:,0],gsw.distance(lonmat.T,latmat.T).T))
# curl=gradient(tauy)[2]/dlon-gradient(taux)[1]/dlat
# era['curl'] = era.taux.copy()
# era['curl'][:]=curl*1e6

skiparr=2

def contmap(date1,date2):
    contit=map.contourf(x,y,sqrt(era['taux']**2+era['tauy']**2).sel(date=slice(date1,date2)).mean(dim='date'),51,cmap=cm.inferno,vmin=0,vmax=0.6,extend='both')
    # map.colorbar(ticks=arange(0,0.8,0.1),label='Wind Stress magnitude [N m$^{-2}$]')
    map.plot(CFlon,CFlat,'w-',latlon=True,linewidth=8)
    map.plot(CFlon,CFlat,'k-',latlon=True,linewidth=2)
    map.quiver(x[::skiparr,::skiparr],y[::skiparr,::skiparr],era['taux'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],era['tauy'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],color='white',scale=2)
    map.fillcontinents(color='lightgrey')
    return contit,map
    # savefig('../../confschools/1802_oceansciences/presentation/figures/'+date1+'-'+date2+'_windstress.pdf',bbox_inches='tight')

labvec=[0,0,0,1]

fig=figure(figsize=(4,9))
ax1=subplot(311)
c1,map=contmap('2014-8-1','2014-10-1')
map.drawmeridians(range(-45,-10,10),linewidth=0.01)
map.drawparallels(range(55,80,10),labels=[1,0,0,0],linewidth=0.01)
text(-48,67,'Summer',fontsize=20)
subplot(312)
c2,map=contmap('2014-10-1','2015-1-1')
map.drawmeridians(range(-45,-10,10),linewidth=0.01)
map.drawparallels(range(55,80,10),labels=[1,0,0,0],linewidth=0.01)
text(-48,67,'Fall',fontsize=20)
cbaxes = fig.add_axes([1.0, 0.3, 0.06, 0.45])
colorbar(c2,ticks=arange(0,0.8,0.1),label='Wind Stress magnitude [N m$^{-2}$]',cax=cbaxes)
subplot(313)
c3,map=contmap('2015-1-01','2015-4-1')
map.drawmeridians(range(-45,-10,10),labels=labvec,linewidth=0.01)
map.drawparallels(range(55,80,10),labels=[1,0,0,0],linewidth=0.01)
text(-48,67,'Winter',fontsize=20)
savefig('../../confschools/1802_oceansciences/presentation/figures/ERA_wstress.pdf',bbox_inches='tight')

contmap('2015-01-01','2015-04-01')

contmap('2015-4-01','2015-9-01')

contmap('2015-9-01','2015-12-01')

contmap('2015-12-01','2016-3-01')

contmap('2016-3-01','2016-8-01')
