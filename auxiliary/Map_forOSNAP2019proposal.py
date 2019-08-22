from aux_funcs import *

osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))

bathdistnew=hstack((-15,-11,-10,bathdist))
bathbathnew=hstack((0,0,100,bathbath))

oeind=osnap.LONGITUDE.values>CFlon[0]
bathind=osnap_bathy['lon'][0]>CFlon[0]


osnap_dist=0*osnap.LATITUDE[oeind]
for ii in range(len(osnap_dist)):
    latii=osnap.LATITUDE.values[oeind][ii]
    lonii=osnap.LONGITUDE.values[oeind][ii]
    osnap_dist[ii]=sw.dist([CFlat[0],latii],[CFlon[0],lonii])[0][0]

osnap_dist_bathy=0*osnap_bathy['lon'][0][bathind]
for ii in range(len(osnap_dist_bathy)):
    osnap_dist_bathy[ii]=sw.dist([CFlat[0],osnap_bathy['lat'][0][bathind][ii]],[CFlon[0],osnap_bathy['lon'][0][bathind][ii]])[0][0]

Lon_tot=hstack((CFlon,MMlon[1:3],ooi_lon['fla'],ooi_lon['flb'],MMlon[3:]))
Lat_tot=hstack((CFlat,MMlat[1:3],ooi_lat['fla'],ooi_lat['flb'],MMlat[3:]))
Lat_tot

plot(Lon_tot,Lat_tot,'o')


MMdist=0*array(MMlon)
for ii in range(len(MMdist)):
    MMdist[ii]=sw.dist([CFlat[0],MMlat[ii]],[CFlon[0],MMlon[ii]])[0][0]

ooifldist=hstack((sw.dist([CFlat[0],ooi_lat['fla']],[CFlon[0],ooi_lon['fla']])[0][0],sw.dist([CFlat[0],ooi_lat['flb']],[CFlon[0],ooi_lon['flb']])[0][0]))

ooifldist

distvec_new=sort(hstack((distvec,ooifldist,MMdist[1:])))
distvec_new

adcpdp_new=hstack((adcpdp,NaN,NaN,NaN,500,500,NaN))


depths['M1'][-2:]=2025
inst['M1'][-1]='AQ'


depths['M2']=[1000,1000,1250,1500,1500,1750,2200,2200,2423,2423]
inst['M2']=['AQ','MC','MC','AQ','MC','MC','AQ','MC','AQ','MC']

depths['M3']=[1000,1000,1250,1500,1500,1750,2250,2250,2557,2557]
inst['M3']=['AQ','MC','MC','AQ','MC','MC','AQ','MC','AQ','MC']

depths['FLA']=[30,40,60,90,130,180,250,350,500,750,1000,1500,1700,1700,2000,2300,2300,2700,2700]
inst['FLA']=hstack((['MC']*12,['AQ','MC','MC','AQ','MC','AQ','MC']))

depths['FLB']=[40,60,90,130,180,250,350,500,750,1000,1500,1800,1800,2100,2100,2400,2400,2850,2850]
inst['FLB']=hstack((['MC']*11,['AQ','MC']*4))

depths['M4']=[1500,1500,1750,2000,2250,2250,2500,2750,2985,2985]
inst['M4']=['AQ','MC','MC','MC','AQ','MC','MC','MC','AQ','MC']


def plotinstpos(axchoice):
        for rr in range(len(adcpdp_new)):
            if rr==0:
                axchoice.plot(distvec_new[rr],adcpdp_new[rr],'^',color='lightgreen',markersize=8,zorder=40,label='ADCP',mec='k')
            else:
                axchoice.plot(distvec_new[rr],adcpdp_new[rr],'^',color='lightgreen',markersize=8,zorder=40,mec='k')
        mm=0
        for key in depths:
            for dd in range(len(depths[key])):
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


osnapvel=osnap.VELO.mean(dim='TIME')[:,oeind].values
osnapvel[isnan(osnapvel)]=-0.2

sal_cols=[(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_pos=[0,0.9,0.96,0.99,1]

contour(osnap_dist,osnap.DEPTH,osnap.PSAL.mean(dim='TIME')[:,oeind],levels=[33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35],cmap=sal_cmap)
ylim(3500,0)
xlim(-15,275)
colorbar()


def plotsec_only():
    f,axx=subplots(1,1,figsize=(10,4.5))
    vel=axx.contourf(osnap_dist,osnap.DEPTH,osnapvel,cmap=cm.RdBu_r,levels=arange(-0.5,0.5,0.025),extend='both')
    # axx.contour(osnap_dist,osnap.DEPTH,osnap.VELO.mean(dim='TIME')[:,oeind],levels=[0],colors='k',linewidth=3)
    axx.contourf(dat.distance,dat.depth,dat['across track velocity'].mean(dim='date').T,arange(-0.5,0.5,0.025),cmap=cm.RdBu_r,extend='both')
    salii=axx.contour(osnap_dist,osnap.DEPTH,osnap.PSAL.mean(dim='TIME')[:,oeind],levels=[33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35],cmap=sal_cmap)
    axx.set_ylim(3500,0)
    axx.set_ylabel('depth [m]')
    colorbar(vel,ax=axx,label='Velocity [m/s]',ticks=arange(-0.4,0.5,0.2))
    axx.fill_between(osnap_dist_bathy,-osnap_bathy['bathy'].flatten()[bathind],[4000]*len(osnap_bathy['lon'].flatten()[bathind]),color='k',zorder=22)
    axx.fill_between(bathdistnew,bathbathnew,3500*ones(len(bathbathnew)),color='k',zorder=22)
    plotinstpos(axx)
    axx.set_xlim(-15,275)
    axx.set_xlabel('distance [km]')
    axx.legend(loc=3).set_zorder(102)
    tf=14
    axx.text(0,-100,'EGCC',fontsize=14)
    axx.text(50,-100,'EGC/IC',fontsize=14)
    axx.text(130,2100,'DWBC',fontsize=14,backgroundcolor='white')
    pp=2
    axx.text(28,100,'34',fontsize=8,color='k',backgroundcolor='white', bbox={'facecolor':'white','alpha':0.75, 'pad':pp})#backgroundcolor='#0C2C84',
    axx.text(250,450,'34.9',fontsize=8,color='k',backgroundcolor='white', bbox={'facecolor':'white','alpha':0.75, 'pad':pp})#'#FFEDA0')
    axx.text(90,300,'34.98',fontsize=8,color='k',backgroundcolor='white',zorder=105, bbox={'facecolor':'white','alpha':0.75, 'pad':pp})#'#D95F0E')
    axx.text(-8,425,'CF1,2,3',color='w',fontsize=tf,zorder=101)
    axx.text(20,800,'CF4',color='w',fontsize=tf,zorder=101)
    axx.text(33,1500,'CF5',color='w',fontsize=tf,zorder=101)
    axx.text(50,2050,'CF6',color='w',fontsize=tf,zorder=101)
    axx.text(70,2160,'CF7',color='w',fontsize=tf,zorder=101)
    axx.text(88,2300,'M1',color='w',fontsize=tf,zorder=101)
    axx.text(113,2700,'M2',color='w',fontsize=tf,zorder=101)
    axx.text(140,2850,'M3',color='w',fontsize=tf,zorder=101)
    axx.text(162,3200,'OOI\nFLA',color='w',fontsize=tf,zorder=101)
    axx.text(198,3325,'OOI\nFLB',color='w',fontsize=tf,zorder=101)
    axx.text(238,3225,'M4',color='w',fontsize=tf,zorder=101)#,backgroundcolor='w')
    savefig(figdir+'sections/OSNAP_EAST_BC_SECTION_forproposal2019.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'sections/OSNAP_EAST_BC_SECTION_forproposal2019.pdf',bbox_inches='tight')


plotsec_only()
