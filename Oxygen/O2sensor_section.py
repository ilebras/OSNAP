from aux_funcs import *

osnap=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.nc')
# dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))

bathdistnew=hstack((-15,-11,-10,bathdist))
bathbathnew=hstack((0,0,100,bathbath))

oeind=osnap.LONGITUDE.values>CFlon[0]
bathind=osnap_bathy['lon'][0]>CFlon[0]


osnap_dist=0*osnap.LATITUDE[oeind]
for ii in range(len(osnap_dist)):
    latii=osnap.LATITUDE.values[oeind][ii]
    lonii=osnap.LONGITUDE.values[oeind][ii]
    osnap_dist[ii]=sw.dist([CFlat[0],latii],[CFlon[0],lonii])[0][0]-5

osnap_dist_bathy=0*osnap_bathy['lon'][0][bathind]
for ii in range(len(osnap_dist_bathy)):
    osnap_dist_bathy[ii]=sw.dist([CFlat[0],osnap_bathy['lat'][0][bathind][ii]],[CFlon[0],osnap_bathy['lon'][0][bathind][ii]])[0][0]

Lon_tot=hstack((CFlon,MMlon[1:3]))#,ooi_lon['fla'],ooi_lon['flb'],MMlon[3:]))
Lat_tot=hstack((CFlat,MMlat[1:3]))#,ooi_lat['fla'],ooi_lat['flb'],MMlat[3:]))


plot(Lon_tot,Lat_tot,'o')


MMdist=0*array(MMlon)
for ii in range(len(MMdist)):
    MMdist[ii]=sw.dist([CFlat[0],MMlat[ii]],[CFlon[0],MMlon[ii]])[0][0]

# ooifldist=hstack((sw.dist([CFlat[0],ooi_lat['fla']],[CFlon[0],ooi_lon['fla']])[0][0],sw.dist([CFlat[0],ooi_lat['flb']],[CFlon[0],ooi_lon['flb']])[0][0]))

distvec_new=sort(hstack((distvec,MMdist[1:-1])))

adcpdp_new=hstack((adcpdp,NaN,NaN,NaN))

osnapvel=osnap.VELO.mean(dim='TIME')[:,oeind].values
osnapvel[isnan(osnapvel)]=-0.2

sal_cols=[(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_pos=[0,0.9,0.96,0.99,1]

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Oxygen/'


depths['M1']=array([50,50,500,500,1000,1000,1500,
                    1500,1750,1750,2025,2025])
inst['M1']=array(['MCO','AQ','AQ','MCO','AQ','MCO','AQ',
                  'MC','AQ','MCO','MCO','AQ'])

depths['M2']=[1000,1000,1500,1500,2000,2000,2200,2200,2423,2423]
inst['M2']=['AQ','MC','AQ','MCO','AQ','MC','AQ','MCO','AQ','MCO']

depths['M3']=[1000,1000,1500,1500,2000,2000,2250,2250,2557,2557]
inst['M3']=['AQ','MCO','AQ','MC','AQ','MCO','AQ','MC','AQ','MCO']

inst['CF2'][4]='MCO'
inst['CF4'][2]='MCO'
inst['CF4'][4]='MCO'

inst['CF5']=['MC', 'WH300', 'MCO', 'tidbit', 'tidbit', 'MCO', 'AQ', 'tidbit','RBR solo', 'MCO', 'AQ', 'MC', 'AQ', 'MCO', 'AQ', 'MC', 'AQ']
inst['CF6']=['CTDp', 'WH300', 'MCO', 'tidbit', 'tidbit', 'MCO', 'AQ', 'tidbit', 'RBR solo', 'MC', 'AQ', 'MCO', 'AQ', 'MC', 'AQ', 'MCO', 'AQ', 'MC','AQ']
inst['CF7']=['RBR XR-420', 'WH300', 'MCO', 'tidbit', 'tidbit', 'MC', 'AQ', 'tidbit', 'RBR solo', 'MCO', 'AQ', 'MC', 'AQ', 'MCO', 'AQ', 'MC','AQ', 'MC', 'AQ']
def plotinstpos(axchoice):
        for rr in range(len(adcpdp_new)):
            if rr==0:
                axchoice.plot(distvec_new[rr],adcpdp_new[rr],'^',color='lightgreen',markersize=8,zorder=40,label='ADCP',mec='k')
            else:
                axchoice.plot(distvec_new[rr],adcpdp_new[rr],'^',color='lightgreen',markersize=8,zorder=40,mec='k')
        mm=0
        for key in depths:
            for dd in range(len(depths[key])):
                    if ('O' in inst[key][dd]):
                        if (mm==0):
                            axchoice.plot(distvec_new[mm],depths[key][dd],'ro',zorder=35,markersize=10,label='Oxygen sensor')
                        else:
                            axchoice.plot(distvec_new[mm],depths[key][dd],'ro',zorder=35,markersize=10,label='')

                    if ('AQ' in inst[key][dd]):
                        if (mm==0):
                            axchoice.plot(distvec_new[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='Current meter')
                        else:
                            axchoice.plot(distvec_new[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='')

                    if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
                        if (mm==0) & (dd==0):
                            axchoice.plot(distvec_new[mm],depths[key][dd],'o',color='#fee090',zorder=38,markersize=4,mec='k',label='T,S recorder')
                        else:
                            axchoice.plot(distvec_new[mm],depths[key][dd],'o',color='#fee090',zorder=38,markersize=4,mec='k',label='')

            axchoice.plot([distvec_new[mm],distvec_new[mm]],[depths[key][0],depths[key][-1]],'k')

            mm+=1



def plotsec_only():
    f,axx=subplots(1,1,figsize=(8,4.5))
    vel=axx.contourf(osnap_dist,osnap.DEPTH,osnapvel,cmap=cm.RdBu_r,levels=arange(-0.5,0.5,0.025),extend='both')
    salii=axx.contour(osnap_dist,osnap.DEPTH,osnap.PSAL.mean(dim='TIME')[:,oeind],levels=[33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35],cmap=sal_cmap)
    axx.set_ylim(3000,0)
    axx.set_ylabel('depth [m]')
    colorbar(vel,ax=axx,label='Velocity [m/s]',ticks=arange(-0.4,0.5,0.2))
    axx.fill_between(osnap_dist_bathy,-osnap_bathy['bathy'].flatten()[bathind],[4000]*len(osnap_bathy['lon'].flatten()[bathind]),color='k',zorder=22)
    axx.fill_between(bathdistnew,bathbathnew,3500*ones(len(bathbathnew)),color='k',zorder=22)
    plotinstpos(axx)
    axx.set_xlim(-15,165)
    axx.set_xlabel('distance [km]')
    axx.legend(loc=3).set_zorder(102)
    tf=14
    axx.text(-2,425,'CF1,2,3',color='w',fontsize=tf,zorder=101)
    axx.text(distvec_new[1]-2,-100,'1',color='k',fontsize=tf,zorder=101)
    axx.text(20,800,'CF4',color='w',fontsize=tf,zorder=101)
    axx.text(distvec_new[3]-2,-100,'2',color='k',fontsize=tf,zorder=101)
    axx.text(33,1500,'CF5',color='w',fontsize=tf,zorder=101)
    axx.text(distvec_new[4]-2,-100,'4',color='k',fontsize=tf,zorder=101)
    axx.text(distvec_new[5]-2,-100,'4',color='k',fontsize=tf,zorder=101)
    axx.text(distvec_new[6]-2,-100,'3',color='k',fontsize=tf,zorder=101)
    axx.text(distvec_new[7]-2,-100,'5',color='k',fontsize=tf,zorder=101)
    axx.text(distvec_new[8]-2,-100,'3',color='k',fontsize=tf,zorder=101)
    axx.text(distvec_new[9]-2,-100,'3',color='k',fontsize=tf,zorder=101)
    axx.text(50,2050,'CF6',color='w',fontsize=tf,zorder=101)
    axx.text(70,2160,'CF7',color='w',fontsize=tf,zorder=101)
    axx.text(88,2300,'M1',color='w',fontsize=tf,zorder=101)
    axx.text(113,2700,'M2',color='w',fontsize=tf,zorder=101)
    axx.text(140,2850,'M3',color='w',fontsize=tf,zorder=101)
    savefig(figdir+'O2sensor_section.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'O2sensor_section.pdf',bbox_inches='tight')


plotsec_only()

for kk in inst:
    print(kk)
    for ll in inst[kk]:
        if 'O' in ll:
