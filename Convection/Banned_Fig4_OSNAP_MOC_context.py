from aux_funcs import *

osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

[uIIW,dIIW,IIW,egic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_xport/IIWtrans_direct.pickle','rb'))
uIIW['osnap']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.PDEN>=d1).where(osnap.PDEN<d2).sum(dim='DEPTH')/1e6
dIIW['osnap']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.PDEN>=d2).where(osnap.PDEN<d3).sum(dim='DEPTH')/1e6

uIIW['osnap'].sum(dim='LONGITUDE').mean(dim='TIME')

dIIW['osnap'].sum(dim='LONGITUDE').mean(dim='TIME')



psi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Streamfunction_201408_201604_2018.nc')
SIGmax_east=ones(21)
MOC_east=ones(21)
for ii in range(21):
    MOC_east[ii]=psi['T_EAST'][ii,:].max()
    SIGmax_east[ii]=psi.LEVEL[psi.T_EAST[ii,:].argmax()]

midlevel=(psi.LEVEL[:-1]+diff(psi.LEVEL)/2).values
nettrans_EAST=psi.T_EAST.diff(dim='LEVEL')
moctrans_EAST=nettrans_EAST.copy()
for ii,dd in enumerate(SIGmax_east):
    moctrans_EAST[ii,midlevel>dd]=-nettrans_EAST[ii,midlevel>dd]
uIIW_eastcheck_MOC=moctrans_EAST.where(psi.LEVEL>=d1).where(psi.LEVEL<d2).sum(dim='LEVEL')
dIIW_eastcheck_MOC=moctrans_EAST.where(psi.LEVEL>=d2).where(psi.LEVEL<d3).sum(dim='LEVEL')

def overturning_check():
    figure(figsize=(12,3))
    psi.T_EAST.diff(dim='LEVEL').T.plot(vmin=-2.5,vmax=2.5,cmap=cm.RdBu_r)
    plot(psi.TIME,SIGmax_east,'k')
    ylim(28,27)
    for dd in [d1,d2,d3]:
        axhline(dd,color='g')

    figure(figsize=(12,3))
    moctrans_EAST.T.plot(vmin=-2.5,vmax=2.5,cmap=cm.RdBu_r)
    plot(psi.TIME,SIGmax_east,'k')
    ylim(28,27)
    for dd in [d1,d2,d3]:
        axhline(dd,color='g')


overturning_check()


from map_funcs import *

def mapmeat(lon_start,lat_start,lon_end,lat_end,lon,lat,bathySmoothed):

    lon_0= (lon_end + lon_start)/2.0
    lat_0= - (abs(lat_end)+abs(lat_start))/2.0


    map = Basemap(lat_0=lat_0,lon_0=lon_0,llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='l',projection='stere')

    x, y = map(lon,lat)
    CS0 = map.contour(x,y,bathySmoothed,[-3000,-2000,-1000],colors='grey')

    map.drawcoastlines()
    map.fillcontinents()

    eastind=osnap.LONGITUDE>-43
    map.plot(osnap.LONGITUDE.values[eastind],osnap.LATITUDE.values[eastind],color='k',latlon=True,linewidth=4)

    return map,x,y,CS0


def pltveldensec(axx,draw_moors='yes'):
        vel=axx.contourf(-osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
        axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[27.53],colors='darkgrey',linewidths=3)
        axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=dbnds,colors='k',linewidths=3)
        axx.fill_between(-osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)

        axx.set_ylabel('depth [m]')
        if draw_moors=='yes':
            axx.axvline(-CFlon[4],color='w',linewidth=4)
            axx.axvline(-CFlon[4],color=cf5col,linewidth=2)
            axx.axvline(-CFlon[5],color='w',linewidth=4)
            axx.axvline(-CFlon[5],color=cf6col,linewidth=2)
            axx.axvline(-CFlon[-1],color='w',linewidth=4)
            axx.axvline(-CFlon[-1],color=m1col,linewidth=2)
            axx.axvline(-ooi_lon['fla'],color='w',linewidth=4)
            axx.axvline(-ooi_lon['fla'],color=oomcol,linewidth=2)
        return vel


def plot_Fig4():
    f=figure(figsize=(11,13))
    outer = gridspec.GridSpec(2, 1, height_ratios = [1.25,1],hspace=0.4)
    gs=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[0], height_ratios=[1.5,1],hspace=0.2)
    gs1=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[1], height_ratios=[1,1],hspace=0.2)
    axx=plt.subplot(gs[0])
    bxx=plt.subplot(gs[1])
    dxx=plt.subplot(gs1[0])
    cxx=plt.subplot(gs1[1])

    tf=16
    cax1=f.add_axes([0.92,0.7,0.015,0.175])
    vel=pltveldensec(axx,'no')
    colorbar(vel,ax=axx,cax=cax1,label='Velocity [m/s]')
    axx.set_xticklabels('')
    # axx.set_title('OSNAP East')
    axx.set_xlim(43.5,9)
    axx.set_ylim(3750,0)
    axx.set_yticks(range(0,3500,1000))

    axx.text(40,3500,'Irminger Basin',color='white',fontsize=tf-2,zorder=101)
    axx.text(27,3500,'Iceland Basin',color='white',fontsize=tf-2,zorder=101)
    axx.text(18.8,2000,'Hatton-Rockall\n    Plateau',color='white',fontsize=tf-2,zorder=101)
    axx.text(12,3000,'Rockall\nTrough',color='white',fontsize=tf-2,zorder=101)
    axx.text(28.5,950,'upper IIW',fontsize=tf-2,zorder=101)
    axx.text(28.5,1350,'deep IIW',fontsize=tf-2,zorder=101)
    axx.text(30,275,'$\sigma_{max}$=27.53',color='k',zorder=101,fontsize=tf-3)
    axx.set_title('a) Mean velocity and density structure across OSNAP East',fontsize=14)

    uIIW_cum=uIIW['osnap'].cumsum(dim='LONGITUDE')
    dIIW_cum=dIIW['osnap'].cumsum(dim='LONGITUDE')


    bxx.plot(-osnap.LONGITUDE,-uIIW_cum.mean(dim='TIME'),color=uppercol,linewidth=3,label='upper IIW')
    bxx.fill_between(-osnap.LONGITUDE,-uIIW_cum.mean(dim='TIME')+0.5*uIIW_cum.std(dim='TIME'),-uIIW_cum.mean(dim='TIME')-0.5*uIIW_cum.std(dim='TIME'),color=uppercol,alpha=0.2)
    bxx.plot(-osnap.LONGITUDE,-dIIW_cum.mean(dim='TIME'),color=deepcol,linewidth=3,label='deep IIW')
    bxx.fill_between(-osnap.LONGITUDE,-dIIW_cum.mean(dim='TIME')+0.5*dIIW_cum.std(dim='TIME'),-dIIW_cum.mean(dim='TIME')-0.5*dIIW_cum.std(dim='TIME'),color=deepcol,alpha=0.2)
    # bxx.plot(osnap.LONGITUDE,lnadw['osnap'].cumsum(dim='LONGITUDE').mean(dim='TIME'),color='grey',label='densest waters')

    bxx.set_xlim(43.5,9)
    bxx.set_ylim(2.5,-11,-2.5)
    bxx.set_yticks(arange(2.5,-11,-2.5))
    bxx.axhline(0,color='k')
    bxx.set_ylabel('[Sv]')
    bxx.set_xlabel('Longitude [$^\circ$W]',fontsize=14)
    bxx.set_title('b) Cumulative transport within IIW layers',fontsize=14)

    dxx.plot(osnap.TIME,uIIW_eastcheck_MOC,'-o',color=uppercol,linewidth=2,label='upper IIW')
    dxx.plot(osnap.TIME,dIIW_eastcheck_MOC,'-o',color=deepcol,linewidth=2,label='deep IIW')
    dxx.plot(osnap.TIME,MOC_east,'-o',color='k',linewidth=2,label='total overturning')
    dxx.legend(fontsize=14,loc=(0.7,1.1))
    dxx.set_ylabel('[Sv]')
    dxx.set_xticklabels('')
    dxx.set_title('c) Overturning transport',fontsize=14)

    cxx.plot(osnap.TIME,uIIW_eastcheck_MOC/MOC_east*100,'-o',color=uppercol,linewidth=2)
    cxx.plot(osnap.TIME,dIIW_eastcheck_MOC/MOC_east*100,'-o',color=deepcol,linewidth=2)
    cxx.set_ylabel('Percentage of overturning')
    cxx.set_ylim(-40,65)
    cxx.set_yticks(arange(-40,70,20))
    cxx.set_yticklabels([str(aa)+'%' for aa in range(-40,70,20)])
    cxx.set_title('d) IIW contribution to overturning',fontsize=14)

    for axi in [cxx,dxx]:
        axi.axhline(0,color='grey')
        axi.set_xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,4,15)])
        axi.xaxis.set_major_locator(years)
        axi.xaxis.set_minor_locator(threemonth)
        axi.yaxis.grid(True)

    cxx.xaxis.set_minor_formatter(monthFMT)
    cxx.xaxis.set_major_formatter(yearFMT)

    savefig(figdir+'MixedLayer/paperfigs/Fig4.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig4.pdf',bbox_inches='tight')

plot_Fig4()


# def plot_Fig4():
#     f=figure(figsize=(11,9))
#     outer = gridspec.GridSpec(2, 1, height_ratios = [2.5,1],hspace=0.4)
#     gs=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[0], height_ratios=[1.5,1],hspace=0.2)
#     axx=plt.subplot(gs[0])
#     bxx=plt.subplot(gs[1])
#     cxx=plt.subplot(outer[1])
#
#     tf=16
#     cax1=f.add_axes([0.92,0.65,0.015,0.2])
#     vel=pltveldensec(axx,'no')
#     colorbar(vel,ax=axx,cax=cax1,label='Velocity [m/s]')
#     axx.set_xticklabels('')
#     # axx.set_title('OSNAP East')
#     axx.set_xlim(43.5,9)
#     axx.set_ylim(3750,0)
#     axx.set_yticks(range(0,3500,1000))
#
#     axx.text(40,3500,'Irminger Basin',color='white',fontsize=tf-2,zorder=101)
#     axx.text(27,3500,'Iceland Basin',color='white',fontsize=tf-2,zorder=101)
#     axx.text(18.8,2000,'Hatton-Rockall\n    Plateau',color='white',fontsize=tf-2,zorder=101)
#     axx.text(12,3000,'Rockall\nTrough',color='white',fontsize=tf-2,zorder=101)
#     axx.text(28.5,950,'upper IIW',fontsize=tf-2,zorder=101)
#     axx.text(28.5,1350,'deep IIW',fontsize=tf-2,zorder=101)
#     axx.text(30,275,'$\sigma_{max}$=27.53',color='k',zorder=101,fontsize=tf-3)
#     axx.set_title('a) Mean velocity and density structure across OSNAP East',fontsize=14)
#
#     uIIW_cum=uIIW['osnap'].cumsum(dim='LONGITUDE')
#     dIIW_cum=dIIW['osnap'].cumsum(dim='LONGITUDE')
#
#
#     bxx.plot(-osnap.LONGITUDE,-uIIW_cum.mean(dim='TIME'),color=uppercol,linewidth=3,label='upper IIW')
#     bxx.fill_between(-osnap.LONGITUDE,-uIIW_cum.mean(dim='TIME')+0.5*uIIW_cum.std(dim='TIME'),-uIIW_cum.mean(dim='TIME')-0.5*uIIW_cum.std(dim='TIME'),color=uppercol,alpha=0.2)
#     bxx.plot(-osnap.LONGITUDE,-dIIW_cum.mean(dim='TIME'),color=deepcol,linewidth=3,label='deep IIW')
#     bxx.fill_between(-osnap.LONGITUDE,-dIIW_cum.mean(dim='TIME')+0.5*dIIW_cum.std(dim='TIME'),-dIIW_cum.mean(dim='TIME')-0.5*dIIW_cum.std(dim='TIME'),color=deepcol,alpha=0.2)
#     # bxx.plot(osnap.LONGITUDE,lnadw['osnap'].cumsum(dim='LONGITUDE').mean(dim='TIME'),color='grey',label='densest waters')
#     bxx.legend(fontsize=14)
#     bxx.set_xlim(43.5,9)
#     bxx.set_ylim(2.5,-11,-2.5)
#     bxx.set_yticks(arange(2.5,-11,-2.5))
#     bxx.axhline(0,color='k')
#     bxx.set_ylabel('[Sv]')
#     bxx.set_xlabel('Longitude [$^\circ$W]',fontsize=14)
#     bxx.set_title('b) Cumulative transport within IIW layers',fontsize=14)
#
#     cxx.plot(osnap.TIME,uIIW_eastcheck_MOC/MOC_east*100,'-o',color=uppercol,linewidth=2)
#     cxx.plot(osnap.TIME,dIIW_eastcheck_MOC/MOC_east*100,'-o',color=deepcol,linewidth=2)
#     cxx.axhline(0,color='grey')
#     cxx.set_xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,4,15)])
#     cxx.set_ylabel('Percentage of overturning')
#     cxx.xaxis.set_major_locator(years)
#     cxx.xaxis.set_minor_locator(threemonth)
#     cxx.xaxis.set_minor_formatter(monthFMT)
#     cxx.xaxis.set_major_formatter(yearFMT)
#     cxx.yaxis.grid(True)
#     cxx.set_ylim(-40,60)
#     cxx.set_yticks(arange(-40,70,20))
#     cxx.set_yticklabels([str(aa)+'%' for aa in range(-40,70,20)])
#     cxx.set_title('c) Evolution of IIW contribution to overturning',fontsize=14)
#
#     savefig(figdir+'MixedLayer/paperfigs/Fig4.png',bbox_inches='tight',dpi=300)
#     savefig(figdir+'MixedLayer/paperfigs/Fig4.pdf',bbox_inches='tight')
#
# plot_Fig4()
