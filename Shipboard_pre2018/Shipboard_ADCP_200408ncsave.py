from aux_funcs import *

# adcp_dist=array(pd.DataFrame.from_csv(datadir+'Shipboard/adcp_distances.dat').index)

with open(datadir+'Shipboard/adcp_distances.dat', 'r') as f:
    reader = csv.reader(f)
    adcp_dist = list(reader)

adcp_dist=[float(dd[0]) for dd in adcp_dist]

adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')
adcp_16=io.loadmat(datadir+'Shipboard/ar07_2016/ar07_2016_vm_adcp_ilebras.mat')

# plot(adcp_14['lon'],adcp_14['lat'],'ko');
# plot(adcp_16['lon'],adcp_16['lat'],'ro');
# plot(360+CFlon,CFlat,'yo');
# ylim([60,60.2])
# xlim([316.75,318])

def acrosstrackname(adic):
    adic['across track velocity']=adic['v_dt']*sin(theta)+adic['u_dt']*cos(theta)

    return adic

adcp_14=acrosstrackname(adcp_14)
adcp_16=acrosstrackname(adcp_16)


def contadcp(distvec,adic,tit):
    figure()
    contourf(distvec,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.8,0.05),cmap=cm.RdBu_r);
    colorbar(label='m/s')
    contour(distvec,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2],colors='k');
    fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    gca().invert_yaxis()
    title('Vessel mounted ADCP across-track velocity, '+tit)
    ylim([800,0])
    xlabel('distance (km)')
    ylabel('depth (m)')
    # savefig('../figures/shipboard/ADCPsection_full_'+tit+'.pdf')
    # xlim([-20,90])
    # savefig('../figures/shipboard/ADCPsection_zoom_'+tit+'.pdf')
shape(adcp_14['across track velocity'])
shape(adcp_dist)
#make an adcp xarray to save as netcdf (8 Apr 2020)

xr14=xr.Dataset({'uacross':(['depth','dist'],adcp_14['across track velocity']),'lon':(['dist'],adcp_14['lon'][0]-360),'lat':(['dist'],adcp_14['lat'][0])},coords={'dist':adcp_dist,'depth':adcp_14['depth'][:,0]})
xr16=xr.Dataset({'uacross':(['depth','dist'],adcp_16['across track velocity']),'lon':(['dist'],adcp_16['lon'][0]-360),'lat':(['dist'],adcp_16['lat'][0])},coords={'dist':adcp_dist[:-3],'depth':adcp_16['depth'][:,0]})

xr14_int=xr14.interp(depth=range(10,450,10))
xr16_int=xr16.interp(depth=range(10,450,10))
xr14_int.to_netcdf(datadir+'Shipboard/netcdf/CFsec_SADCP_14_10m.nc')
xr16_int.to_netcdf(datadir+'Shipboard/netcdf/CFsec_SADCP_16_10m.nc')
xr_all=xr.concat([xr14_int,xr16_int],dim='occ')

xr_all['occ']=['2014','2016']
xr_all.to_netcdf(datadir+'Shipboard/netcdf/CFsec_SADCP_1416_10m.nc')


contadcp(adcp_dist,adcp_14,'2014')


contadcp(adcp_dist[:-3],adcp_16,'2016')
