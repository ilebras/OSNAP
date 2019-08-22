from aux_funcs import *

osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))

[uIIW,dIIW,IIW,egic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_xport/IIWtrans_direct.pickle','rb'))

uIIW['osnap']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.PDEN>=d1).where(osnap.PDEN<d2).sum(dim='DEPTH')/1e6

dIIW['osnap']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.PDEN>=d2).where(osnap.PDEN<d3).sum(dim='DEPTH')/1e6
lnadw={}
lnadw['osnap']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.PDEN>=d3).sum(dim='DEPTH')/1e6


uIIW['osnap fulltrans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=ooi_lon['fla']).where(osnap.PDEN>=d1).where(osnap.PDEN<d2).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6
dIIW['osnap fulltrans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=ooi_lon['fla']).where(osnap.PDEN>=d2).where(osnap.PDEN<d3).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6

uIIW['osnap bctrans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=CFlon[-1]).where(osnap.PDEN>=d1).where(osnap.PDEN<d2).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6
dIIW['osnap bctrans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=CFlon[-1]).where(osnap.PDEN>=d2).where(osnap.PDEN<d3).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6

egic['osnap bctrans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>=CFlon[2]).where(osnap.LONGITUDE<=CFlon[-1]).where(osnap.PDEN<27.8).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6

uIIW['osnap east trans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=-8).where(osnap.PDEN>=d1).where(osnap.PDEN<d2).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6
dIIW['osnap east trans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=-8).where(osnap.PDEN>=d2).where(osnap.PDEN<d3).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6

uIIW_cum=uIIW['osnap'].cumsum(dim='LONGITUDE').mean(dim='TIME')


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
        axx.contour(-osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[27.53],colors='r',linewidths=3)
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
    # f,[cxx,axx,bxx]=subplots(3,1,sharex=True,figsize=(15,8),gridspec_kw = {'height_ratios':[1.5,1]})
    f=figure(figsize=(13,6))
    outer = gridspec.GridSpec(2, 1, height_ratios = [1.5,1],hspace=0.2)
    # gs0=gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = outer[0], width_ratios=[1,2])
    # gs1=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[1], height_ratios=[1.5,1],hspace=0.1)
    # ax=plt.subplot(gs0[0])
    axx=plt.subplot(outer[0])
    bxx=plt.subplot(outer[1])
    #

    tf=16
    # title('Mean velocity with IIW density bounds')
    cxx=0.92
    cy=0.35
    cax1=f.add_axes([cxx,0.5,0.015,cy])
    vel=pltveldensec(axx,'no')
    colorbar(vel,ax=axx,cax=cax1,label='Velocity [m/s]')
    axx.set_xticklabels('')
    # axx.set_title('OSNAP East')
    axx.set_xlim(43.5,9)
    axx.set_ylim(3750,0)
    axx.set_yticks(range(0,3500,1000))
    axx.text(40,3500,'Irminger Basin',color='white',fontsize=16,zorder=101)
    axx.text(27,3500,'Iceland Basin',color='white',fontsize=16,zorder=101)
    axx.text(18.8,2000,'Hatton-Rockall\n    Plateau',color='white',fontsize=16,zorder=101)
    axx.text(12,3000,'Rockall\nTrough',color='white',fontsize=16,zorder=101)
    axx.text(28.5,950,'upper IIW',fontsize=tf-2,zorder=101)
    axx.text(28.5,1400,'deep IIW',fontsize=tf-2,zorder=101)
    axx.text(30,275,'$\sigma_{max}$=27.53',color='r',zorder=101,fontsize=tf-3)
    # axx.text(35,450,'27.67',backgroundcolor='w',zorder=101)
    # axx.text(35,800,'27.73',backgroundcolor='w',zorder=101)
    # axx.text(35,1400,'27.77',backgroundcolor='w',zorder=101)

    uIIW_cum=uIIW['osnap'].cumsum(dim='LONGITUDE')
    dIIW_cum=dIIW['osnap'].cumsum(dim='LONGITUDE')


    bxx.plot(-osnap.LONGITUDE,uIIW_cum.mean(dim='TIME'),color=uppercol,linewidth=3,label='upper IIW')
    bxx.fill_between(-osnap.LONGITUDE,uIIW_cum.mean(dim='TIME')-0.5*uIIW_cum.std(dim='TIME'),uIIW_cum.mean(dim='TIME')+0.5*uIIW_cum.std(dim='TIME'),color=uppercol,alpha=0.2)
    bxx.plot(-osnap.LONGITUDE,dIIW_cum.mean(dim='TIME'),color=deepcol,linewidth=3,label='deep IIW')
    bxx.fill_between(-osnap.LONGITUDE,dIIW_cum.mean(dim='TIME')-0.5*dIIW_cum.std(dim='TIME'),dIIW_cum.mean(dim='TIME')+0.5*dIIW_cum.std(dim='TIME'),color=deepcol,alpha=0.2)
    # bxx.plot(osnap.LONGITUDE,lnadw['osnap'].cumsum(dim='LONGITUDE').mean(dim='TIME'),color='grey',label='densest waters')
    bxx.legend(loc=(1.01,0.45),fontsize=14)
    bxx.set_xlim(43.5,9)
    bxx.set_yticks(arange(-2.5,11,2.5))
    bxx.axhline(0,color='k')
    bxx.set_ylabel('Cumulative transport [Sv]')
    bxx.set_xlabel('Longitude [$^\circ$W]')
    savefig(figdir+'MixedLayer/paperfigs/Fig4.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig4.pdf',bbox_inches='tight')

plot_Fig4()




def plotsec():
    # f,[cxx,axx,bxx]=subplots(3,1,sharex=True,figsize=(15,8),gridspec_kw = {'height_ratios':[1.5,1]})
    f=figure(figsize=(15,12))
    outer = gridspec.GridSpec(2, 1, height_ratios = [1,1.75],hspace=0.2)
    gs0=gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = outer[0], width_ratios=[1,2])
    gs1=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[1], height_ratios=[1.5,1],hspace=0.1)
    ax=plt.subplot(gs0[0])
    axx=plt.subplot(gs1[0])
    bxx=plt.subplot(gs1[1])


    pltveldensec(ax)
    ax.set_ylim(3100,0)
    ax.set_xlim(43.5,39)
    ax.set_ylabel('depth [m]')
    tf=16
    ax.text(42.75,1200,'CF5',color=cf5col,fontsize=tf,zorder=101,backgroundcolor='w')
    ax.text(42.5,1800,'CF6',color=cf6col,fontsize=tf,zorder=101,backgroundcolor='w')
    ax.text(41.7,2300,'M1',color=m1col,fontsize=tf,zorder=101,backgroundcolor='w')
    ax.text(40.5,2950,'OOI',color=oomcol,fontsize=tf,zorder=101,backgroundcolor='w')
    ax.text(41.1,350,'upper IIW',fontsize=tf-2,zorder=101)
    ax.text(41.1,950,'deep IIW',fontsize=tf-2,zorder=101)
    ax.set_xlabel('Longitude [$^\circ$W]')

    bx=plt.subplot(gs0[1])
    map=SimpleMap()
    x1,y1=map(-27,59)
    bx.text(x1,y1,'OSNAP East',backgroundcolor='white',fontsize=14)

    # title('Mean velocity with IIW density bounds')
    cxx=0.92
    cy=0.2
    cax1=f.add_axes([cxx,0.33,0.015,cy])
    vel=pltveldensec(axx,'no')
    colorbar(vel,ax=axx,cax=cax1,label='Velocity [m/s]')
    axx.set_xticklabels('')
    # axx.set_title('OSNAP East')
    axx.set_xlim(43.5,9)
    axx.set_ylim(3750,0)
    axx.text(40,3500,'Irminger Basin',color='white',fontsize=16,zorder=101)
    axx.text(27,3500,'Iceland Basin',color='white',fontsize=16,zorder=101)
    axx.text(18.8,2000,'Hatton-Rockall\n    Plateau',color='white',fontsize=16,zorder=101)
    axx.text(12,3000,'Rockall\nTrough',color='white',fontsize=16,zorder=101)
    axx.text(27.5,950,'upper IIW',fontsize=tf-2,zorder=101)
    axx.text(27.5,1400,'deep IIW',fontsize=tf-2,zorder=101)
    axx.text(30,250,'$\sigma_{max}$=27.53',color='r',zorder=101)
    axx.text(35,450,'27.67',backgroundcolor='w',zorder=101)
    axx.text(35,800,'27.73',backgroundcolor='w',zorder=101)
    axx.text(35,1400,'27.77',backgroundcolor='w',zorder=101)

    uIIW_cum=uIIW['osnap'].cumsum(dim='LONGITUDE')
    dIIW_cum=dIIW['osnap'].cumsum(dim='LONGITUDE')


    bxx.plot(-osnap.LONGITUDE,uIIW_cum.mean(dim='TIME'),color=uppercol,linewidth=2,label='upper IIW')
    bxx.fill_between(-osnap.LONGITUDE,uIIW_cum.mean(dim='TIME')-uIIW_cum.std(dim='TIME'),uIIW_cum.mean(dim='TIME')+uIIW_cum.std(dim='TIME'),color=uppercol,alpha=0.2)
    bxx.plot(-osnap.LONGITUDE,dIIW['osnap'].cumsum(dim='LONGITUDE').mean(dim='TIME'),color=deepcol,linewidth=2,label='deep IIW')
    bxx.fill_between(-osnap.LONGITUDE,dIIW_cum.mean(dim='TIME')-dIIW_cum.std(dim='TIME'),dIIW_cum.mean(dim='TIME')+dIIW_cum.std(dim='TIME'),color=deepcol,alpha=0.2)
    # bxx.plot(osnap.LONGITUDE,lnadw['osnap'].cumsum(dim='LONGITUDE').mean(dim='TIME'),color='grey',label='densest waters')
    bxx.legend(loc=(1.01,0.2),fontsize=14)
    bxx.set_xlim(43.5,9)
    bxx.axhline(0,color='k')
    bxx.set_ylabel('Cumulative transport [Sv]')
    bxx.set_xlabel('Longitude [$^\circ$W]')
    savefig(figdir+'MixedLayer/paperfigs/OSNAPEAST_secpluscum.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/OSNAPEAST_secpluscum.pdf',bbox_inches='tight')

plotsec()



def plotsec_only():
    f,axx=subplots(1,1,figsize=(18,4.5))
    vel=axx.contourf(osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
    axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[d1,d2,d3],colors='k')
    axx.set_ylim(3000,0)
    axx.set_ylabel('depth [m]')
    colorbar(vel,ax=axx,label='Velocity [m/s]')
    axvline(CFlon[4],color=cf5col,linewidth=2)
    axvline(CFlon[5],color=cf6col,linewidth=2)
    axvline(CFlon[-1],color=m1col,linewidth=2)
    axvline(ooi_lon['fla'],color=oomcol,linewidth=2)
    axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
    axx.set_xlim(-43,-9)
    title('Mean velocity with IIW density bounds')
    # savefig(figdir+'MixedLayer/geotrans/OSNAPeast_meansection_wIIW.png',bbox_inches='tight')

plotsec_only()


XXXXXXXXXXXXXXXXXX

#
# def plotsec():
#     f,axx=subplots(1,1,figsize=(20,4))
#     vel=axx.contourf(osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.mean(dim='TIME'),20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[27.5],colors='g')
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=[d1,d2,d3],colors='k')
#     axx.set_ylim(3000,0)
#     axx.set_ylabel('depth [m]')
#     colorbar(vel,ax=axx,label='Velocity [m/s]')
#     # axvline(CFlon[4],color=cf5col,linewidth=3)
#     # axvline(CFlon[5],color=cf6col,linewidth=3)
#     # axvline(CFlon[-1],color=m1col,linewidth=3)
#     # axvline(ooi_lon['fla'],color=oomcol,linewidth=3)
#     axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
#     axx.set_xlim(-53,-9)
#     title('Mean velocity with IIW density bounds')
#     savefig(figdir+'MixedLayer/geotrans/OSNAPall_meansection_wIIW.png',bbox_inches='tight')
#
#
# plotsec()
#
# def plotsec_ALL():
#     for ii in range(21):
#         f,axx=subplots(1,1,figsize=(20,4))
#         vel=axx.contourf(osnap.LONGITUDE,osnap.DEPTH,osnap.VELO[ii,:,:],20,cmap=cm.RdBu_r,vmin=-0.4,vmax=0.4,extend='both')
#         dcont=axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN[ii,:,:],levels=[d1,d2,d3],colors='k',linewidths=3)
#         # clabel(dcont)
#         axx.set_ylim(2000,0)
#         axx.set_ylabel('depth [m]')
#         colorbar(vel,ax=axx,label='Velocity [m/s]')
#         axvline(CFlon[4],color=cf5col,linewidth=3)
#         axvline(CFlon[5],color=cf6col,linewidth=3)
#         axvline(CFlon[-1],color=m1col,linewidth=3)
#         axvline(ooi_lon['fla'],color=oomcol,linewidth=3)
#         axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
#         axx.set_xlim(-43,-9)
#         title(str(osnap.TIME[ii].values)[:10])
#
# plotsec_ALL()


def plot_context():
    labvec=['upper IIW','deep IIW']
    for ii,var in enumerate([uIIW,dIIW]):
        figure(figsize=(12,3))
        var['osnap fulltrans'].plot(marker='o',label='Transport up to OOI mooring')
        var['osnap bctrans'].plot(marker='o',label='Transport within CF array')
        var['osnap east trans'].plot(marker='o',label='Total OSNAP EAST transport')
        title(labvec[ii])
        legend(loc=(1.05,0.3))
        xlabel('date')
        ylabel('Transport [Sv]')
        savefig(figdir+'MixedLayer/geotrans/OSNAPEast_context_IIWtrans_'+labvec[ii][:4]+'.png',bbox_inches='tight')

plot_context()

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1808lpfilt.pickle','rb'))

def savefilt(field):
    klist=list(field.keys())
    for key in klist:
        field[key+' filt']=sig.filtfilt(B,A,field[key])

    return field

egic=savefilt(egic)
uIIW=savefilt(uIIW)
dIIW=savefilt(dIIW)

def comp_ud():
    figure(figsize=(12,3))
    plot(dat.date,uIIW['trans filt'],linestyle='--',color='C0',linewidth=3,label='uIIW direct calc')
    uIIW['osnap bctrans'].plot(marker='o',color='C0',label='uIIW OSNAP gridded')
    plot(dat.date,dIIW['trans filt'],linestyle='--',color='C1',linewidth=3,label='dIIW direct calc')
    dIIW['osnap bctrans'].plot(marker='o',color='C1',label='dIIW OSNAP gridded')
    legend(loc=(1.05,0.1))
    title('IIW transport up to M1')
    ylabel('Transport [Sv]')
    xlabel('date')
    savefig(figdir+'MixedLayer/geotrans/IIW_trans_comp_withinBC_nogeo.png',bbox_inches='tight')

    figure(figsize=(12,3))
    plot(dat.date,egic['trans filt'],linestyle='--',color='k',linewidth=3,label='EGIC direct calc')
    egic['osnap bctrans'].plot(marker='o',color='k',label='EGIC OSNAP gridded')
    legend(loc=(1.05,0.1))
    title('EGIC transport up to M1')
    ylabel('Transport [Sv]')
    xlabel('date')
    savefig(figdir+'MixedLayer/geotrans/EGIC_trans_comp_withinBC_nogeo.png',bbox_inches='tight')


comp_ud()
psi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Streamfunction_201408_201604_2018.nc')

MOC_east=ones(21)
SIGmax_east=ones(21)
MOC_all=ones(21)
SIGmax_all=ones(21)
for ii in range(21):
    MOC_east[ii]=psi.T_EAST[ii,:].max()
    MOC_all[ii]=psi.T_ALL[ii,:].max()
    SIGmax_east[ii]=psi.LEVEL[psi.T_EAST[ii,:].argmax()]
    SIGmax_all[ii]=psi.LEVEL[psi.T_ALL[ii,:].argmax()]

def plotpsi_ALL():
    for ii in range(21):
        figure(figsize=(2,4))
        plot(psi.T_EAST[ii,:],psi.LEVEL);
        plot(MOC_east[ii],SIGmax_east[ii],'o');
        title(str(psi.TIME[ii].values)[:10])
        ylim(28,27)
        [axhline(dd,color='grey') for dd in [d1,d2,d3]]



nettrans_ALL=psi.T_ALL.diff(dim='LEVEL')
nettrans_EAST=psi.T_EAST.diff(dim='LEVEL')

midlevel=(psi.LEVEL[:-1]+diff(psi.LEVEL)/2).values

moctrans_ALL=nettrans_ALL.copy()
for ii,dd in enumerate(SIGmax_all):
    moctrans_ALL[ii,midlevel<dd]=-nettrans_ALL[ii,midlevel<dd]

moctrans_EAST=nettrans_EAST.copy()
for ii,dd in enumerate(SIGmax_east):
    moctrans_EAST[ii,midlevel<dd]=-nettrans_EAST[ii,midlevel<dd]


uIIW_totcheck=-nettrans_ALL.where(psi.LEVEL>=d1).where(psi.LEVEL<d2).sum(dim='LEVEL')
uIIW_eastcheck=-nettrans_EAST.where(psi.LEVEL>=d1).where(psi.LEVEL<d2).sum(dim='LEVEL')

dIIW_totcheck=-nettrans_ALL.where(psi.LEVEL>=d2).where(psi.LEVEL<d3).sum(dim='LEVEL')
dIIW_eastcheck=-nettrans_EAST.where(psi.LEVEL>=d2).where(psi.LEVEL<d3).sum(dim='LEVEL')

uIIW_totcheck_MOC=-moctrans_ALL.where(psi.LEVEL>=d1).where(psi.LEVEL<d2).sum(dim='LEVEL')
uIIW_eastcheck_MOC=-moctrans_EAST.where(psi.LEVEL>=d1).where(psi.LEVEL<d2).sum(dim='LEVEL')

dIIW_totcheck_MOC=-moctrans_ALL.where(psi.LEVEL>=d2).where(psi.LEVEL<d3).sum(dim='LEVEL')
dIIW_eastcheck_MOC=-moctrans_EAST.where(psi.LEVEL>=d2).where(psi.LEVEL<d3).sum(dim='LEVEL')

def ComptoMOC():
    figure(figsize=(12,3))
    plot(psi.TIME,MOC_east,label='OSNAP East overturning')
    uIIW_eastcheck_MOC.plot(label='uIIW contribution')
    dIIW_eastcheck_MOC.plot(label='dIIW contribution')
    xlabel('')
    legend(loc=(1.05,0.3))
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    savefig(figdir+'MixedLayer/geotrans/Comp_MOC_contrib_upperdeep_east.png',bbox_inches='tight')

uIIW_eastcheck_MOC.mean()
uIIW_eastcheck_MOC.std()

uIIW_eastcheck_MOC.mean()/MOC_east.mean()

dIIW_eastcheck_MOC.mean()

dIIW_eastcheck_MOC.std()

dIIW_eastcheck_MOC.mean()/MOC_east.mean()

ComptoMOC()

def ComptoMOC():
    figure(figsize=(12,3))
    plot(psi.TIME,MOC_all,label='OSNAP total overturning')
    uIIW_totcheck_MOC.plot(label='uIIW contribution')
    dIIW_totcheck_MOC.plot(label='dIIW contribution')
    xlabel('')
    legend(loc=(1.05,0.3))
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    savefig(figdir+'MixedLayer/geotrans/Comp_MOC_contrib_upperdeep_all.png',bbox_inches='tight')

ComptoMOC()

figure(figsize=(12,3))
plot(psi.TIME,(uIIW_eastcheck_MOC/MOC_east)*100,'o-')
plot(psi.TIME,(dIIW_eastcheck_MOC/MOC_east)*100,'o-')
ylabel('Percent contribution to total MOC EAST')
grid('on')
savefig(figdir+'MixedLayer/geotrans/CompPercent_MOC_contrib_upperdeep_east.png',bbox_inches='tight')

figure(figsize=(12,3))
plot(psi.TIME,(uIIW_totcheck_MOC/MOC_all)*100,'o-')
plot(psi.TIME,(dIIW_totcheck_MOC/MOC_all)*100,'o-')
ylabel('Percent contribution to total MOC EAST')
grid('on')
savefig(figdir+'MixedLayer/geotrans/CompPercent_MOC_contrib_upperdeep_all.png',bbox_inches='tight')

def compnet_upper():
    uIIW_eastcheck.plot(figsize=(12,3),label='east')
    uIIW_totcheck.plot(label='total')
    xlabel('')
    legend()
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    title('Net upper IIW transport')
    savefig(figdir+'MixedLayer/geotrans/Comp_EastTot_Net_IIW_upper.png',bbox_inches='tight')

compnet_upper()

def compmoc_upper():
    uIIW_eastcheck_moc.plot(figsize=(12,3),label='east')
    uIIW_totcheck_moc.plot(label='total')
    xlabel('')
    legend()
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    title('upper IIW overturning transport')
    savefig(figdir+'MixedLayer/geotrans/Comp_EastTot_MOC_IIW_upper.png',bbox_inches='tight')

compmoc_upper()

def compnet_deep():
    dIIW_eastcheck.plot(figsize=(12,3),label='east')
    dIIW_totcheck.plot(label='total')
    xlabel('')
    legend()
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    title('Net deep IIW transport')
    savefig(figdir+'MixedLayer/geotrans/Comp_EastTot_Net_IIW_deep.png',bbox_inches='tight')

compnet_deep()

def compmoc_deep():
    dIIW_eastcheck_MOC.plot(figsize=(12,3),label='east')
    dIIW_totcheck_MOC.plot(label='total')
    xlabel('')
    legend()
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    title('deep IIW overturning transport')
    savefig(figdir+'MixedLayer/geotrans/Comp_EastTot_MOC_IIW_deep.png',bbox_inches='tight')

compmoc_deep()

def compMOC_upper():
    uIIW_eastcheck.plot(figsize=(12,3),label='net')
    uIIW_eastcheck_MOC.plot(label='overturning')
    legend()
    xlabel('')
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    title('upper IIW overturning transport')
    savefig(figdir+'MixedLayer/geotrans/Comp_Net_Overturning_IIW_upper.png',bbox_inches='tight')

compMOC_upper()

def plot_MOC_sigmaspace():
    figure(figsize=(12,5))
    pcolormesh(psi.TIME,psi.LEVEL[350:-1],psi.T_EAST.diff(dim='LEVEL')[:,350:].T,vmin=-2,vmax=2,cmap=cm.RdBu_r)
    colorbar(label='Transport [Sv]')
    [axhline(dd,color='k') for dd in [d1,d2,d3]]
    plot(psi.TIME,SIGmax_east,'o-',color='grey')
    ylabel('Density')
    ylim(28,27)
    title('OSNAP EAST overturning in density space')
    savefig(figdir+'MixedLayer/geotrans/Overturning_densityspace_East.png',bbox_inches='tight')

plot_MOC_sigmaspace()


def plot_MOC_sigmaspace():
    figure(figsize=(12,5))
    pcolormesh(psi.TIME,psi.LEVEL[350:-1],psi.T_ALL.diff(dim='LEVEL')[:,350:].T,vmin=-2,vmax=2,cmap=cm.RdBu_r)
    colorbar(label='Transport [Sv]')
    [axhline(dd,color='k') for dd in [d1,d2,d3]]
    plot(psi.TIME,SIGmax_all,'o-',color='grey')
    ylabel('Density')
    title('Total OSNAP overturning in density space')
    ylim(28,27)
    savefig(figdir+'MixedLayer/geotrans/Overturning_densityspace_All.png',bbox_inches='tight')

plot_MOC_sigmaspace()

def plotpsi():
    figure(figsize=(4,8))
    plot(psi.T_EAST.T,psi.LEVEL);
    plot(MOC_east,SIGmax_east,'o');
    ylim(28,27)
    [axhline(dd,color='grey') for dd in [d1,d2,d3]]
    figure(figsize=(12,3))
    plot(psi.TIME,SIGmax_east,'o-')
    [axhline(dd,color='grey') for dd in [d1,d2,d3]]
    title('Density of maximum overturning and IIW density bounds')
    ylabel('potential density [kg/$m^3$]')
    ylim(28,27)
    savefig(figdir+'MixedLayer/geotrans/MaxOverturning_time.png',bbox_inches='tight')
    figure(figsize=(12,3))
    plot(psi.TIME,MOC_east,'o-')
    title('Overturning transport')
    ylabel('Transport [Sv]')
    savefig(figdir+'MixedLayer/geotrans/OSNAPEast_Overturning_time.png',bbox_inches='tight')

plotpsi()
