############################################################################################################
###################### Load CTD and LADCP data, plot together and get transports ###########################
############################################################################################################

from aux_funcs import *

LADCP=pickle.load(open('../pickles/Shipboard/LADCP_xarray.pickle','rb'))
CTD=pickle.load(open('../pickles/Shipboard/CTD_xarray.pickle','rb'))

LADCP
CTD
grdat=xr.merge([LADCP,CTD])


d1=0
d2=30

figure(figsize=(10,6))
contourf(grdat['distance [km]'][d1:d2],grdat['pressure [db]'],grdat['across track velocity [m/s]'][:,d1:d2,:].mean(dim='occupation'),1001,vmin=-0.3,vmax=0.3,cmap=cm.RdBu_r,extend="both")
colorbar(ticks=arange(-0.4,0.5,0.1),label='[m/s]')
contour(grdat['distance [km]'][d1:d2],grdat['pressure [db]'],grdat['across track velocity [m/s]'][:,d1:d2,:].mean(dim='occupation'),arange(-0.4,0.5,0.1),colors='k')
denlab=contour(grdat['distance [km]'][d1:d2],grdat['pressure [db]'],grdat['density'][:,d1:d2,:].mean(dim='occupation'),levels=[27.8],colors='k',zorder=100,linewidths=3)
clabel(denlab,fmt='%1.1f')
fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
title('Mean across track velocity from LADCP',fontsize=20)
xlabel('distance [km]')
ylabel('pressure [db]')
ylim([3e3,0])
xlim([0,250])
savefig('../figures/Shipboard_forPenny/wden/LADCP_uacross_mean_wden.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/wden/LADCP_uacross_mean_wden.pdf',bbox_inches='tight')


mid_dist=hstack((0,(diff(grdat['distance [km]'][d1:d2][:-1])+diff(grdat['distance [km]'][d1:d2])[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(grdat['pressure [db]'])-1,len(grdat['occupation']),1])),(0,2,1))
depthdiffmat=transpose((tile(diff(grdat['pressure [db]']),[len(grdat['distance [km]'][d1:d2]),len(grdat['occupation']),1])),(2,0,1))

trans=(grdat.where(grdat['density']>=27.8)['across track velocity [m/s]'][:-1,d1:d2,:]*depthdiffmat*middistmat/1e3).sum('pressure [db]').cumsum('distance [km]').min('distance [km]')

plot(trans,'o')
trans_tit=[str(round(tt*100)/100) for tt in trans.values]

g=grdat['across track velocity [m/s]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-0.4,vmax=0.4,cmap=cm.RdBu_r,figsize=(15,12))
ylim([3e3,0]);
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii]+': '+trans_tit[ii]+' Sv',fontsize=20)
    ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    ax.contour(grdat['distance [km]'][d1:d2],grdat['pressure [db]'],grdat['across track velocity [m/s]'][:,d1:d2,ii],arange(-0.4,0.5,0.1),colors='k')
    denlab=ax.contour(grdat['distance [km]'][d1:d2],grdat['pressure [db]'],grdat['density'][:,d1:d2,ii],levels=[27.8],colors='k',linewidths=3,zorder=100)
    # clabel(denlab,fmt='%1.1f',manual=[()])
savefig('../figures/Shipboard_forPenny/wden/LADCP_panels_uacross_wden.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/wden/LADCP_panels_uacross_wden.pdf',bbox_inches='tight')
