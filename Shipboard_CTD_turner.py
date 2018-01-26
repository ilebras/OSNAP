from aux_funcs import *

datin=pickle.load(open('../pickles/Shipboard/CTD_xarray.pickle','rb'))

prsvec=range(10,3125,10)

SA=datin['salinity'].copy()
for ii in range(6):
    SA[:,:,ii]=gsw.SA_from_SP(datin['salinity'][:,:,ii].T,prsvec,CFlon[3],CFlat[3]).T

CT=gsw.CT_from_pt(SA,datin['temperature [$^\circ$ C]'])

SAint=SA[:,:-1,:]+diff(SA,axis=1)
CTint=CT[:,:-1,:]+diff(CT,axis=1)
alphaint=CTint.copy()
betaint=SAint.copy()
for ii in range(6):
    alphaint[:,:,ii]=gsw.alpha(SAint[:,:,ii].T,CTint[:,:,ii].T,prsvec).T
    betaint[:,:,ii]=gsw.beta(SAint[:,:,ii].T,CTint[:,:,ii].T,prsvec).T

alphaT=-alphaint*diff(datin['temperature [$^\circ$ C]'],axis=1)
betaS=betaint*diff(datin['salinity'],axis=1)

R=alphaT/betaS

dendiff=diff(datin['density'],axis=1)

turner=arctan2((alphaT-betaS)*dendiff/abs(dendiff),(alphaT+betaS)*dendiff/abs(dendiff))*180/pi

turner_noden=arctan2((alphaT-betaS),(alphaT+betaS))*180/pi

datin

grdat=xr.Dataset({'turner angle': (['pressure [db]', 'distance [km]', 'occupation'],  turner),
                'turner angle noden': (['pressure [db]', 'distance [km]', 'occupation'],  turner_noden),
                'density gradient': (['pressure [db]', 'distance [km]', 'occupation'],  dendiff)},
                coords={'distance [km]': datin['distance [km]'][:-1]+diff(datin['distance [km]'])/2,
                        'pressure [db]': datin['pressure [db]'],
                        'occupation': datin.occupation})

grdat

g=grdat['turner angle'].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-90,vmax=90, extend='both',cmap=cm.RdBu_r,figsize=(15,12))
xlim([0,100])
ylim([2000,0])
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    # contourf(grdat['distance'],grdat['depth'],grdat['turner angle'][:,:,ii],arange(-90,100,5),cmap=cm.RdBu_r,zorder=20)
    ax.contour(grdat['distance [km]'],grdat['pressure [db]'],grdat['turner angle'][:,:,ii],levels=[-90,-45,0,45,90],colors='k',linewidths=3,zorder=100)
    # xlabel('distance [km]')
    # ylabel('pressure [db]')
savefig('../figures/shipboard_turb/Turner_horiz_ctd.png',bbox_inches='tight')

g=grdat['turner angle noden'].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-180,vmax=180, extend='both',cmap=cm.RdBu_r,figsize=(15,12))
xlim([0,100])
ylim([2000,0])
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    # contourf(grdat['distance'],grdat['depth'],grdat['turner angle'][:,:,ii],arange(-90,100,5),cmap=cm.RdBu_r,zorder=20)
    ax.contour(grdat['distance [km]'],grdat['pressure [db]'],grdat['turner angle'][:,:,ii],levels=[-180,-90,0,90,180],colors='k',linewidths=3,zorder=100)
    # xlabel('distance [km]')
    # ylabel('pressure [db]')
savefig('../figures/shipboard_turb/Turner_horiz_noden_ctd.png',bbox_inches='tight')
