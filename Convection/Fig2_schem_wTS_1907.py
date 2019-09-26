from aux_funcs import *

dendat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/density_gridded_props_cf5-oom_from5m.nc')

dat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_5m.nc')

dendat.thickness.sel(den=slice(d1,d3)).mean(dim='date').plot()

def TSmrange(ii,axx,da1,da2,col):
    buf=0.02
    salsel=dendat.psal[ii,:,:].sel(date=slice(da1,da2)).sel(den=slice(d1-buf,d3+buf))
    tmpsel=dendat.ptmp[ii,:,:].sel(date=slice(da1,da2)).sel(den=slice(d1-buf,d3+buf))
    thicksel=dendat.thickness[ii,:,:].sel(date=slice(da1,da2)).sel(den=slice(d1-buf,d3+buf))

    salvar=(salsel*thicksel).sum(dim='date')/thicksel.sum(dim='date')
    tmpvar=(tmpsel*thicksel).sum(dim='date')/thicksel.sum(dim='date')
    thickvar=thicksel.mean(dim='date')

    tsm=axx.scatter(salvar.values.flatten(),
                 tmpvar.values.flatten(),thickvar.values.flatten(),
                 color=col,linewidth=3,label='',zorder=3,alpha=0.8,edgecolor=col)


def TScomp(axx,dd1,dd2):
    TSmrange(3,axx,dd1,dd2,ooicol)
    TSmrange(2,axx,dd1,dd2,m1col)
    TSmrange(1,axx,dd1,dd2,cf6col)
    TSmrange(0,axx,dd1,dd2,cf5col)
    # axx.set_title(dd1+' to '+dd2)
    den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',linewidths=3,levels=[d1,d2,d3])
    axx.set_ylim(3.2,5)
    axx.set_xlim(34.865,34.98)
    return den


ooicol='#E45D2B'
m1col='#E4942B'
cf5col='#236493'
cf6col='#1E9E64'

colvec=[cf5col,cf6col,m1col,ooicol]

colors = [(12,44,132),(78,179,211),(237,248,177) ,(247,104,161),(129,15,124)]#(237,248,177),,
sal_cmap = make_cmap(colors,position=[0,0.9,0.92,0.99,1],bit=True)#0.9,


def plotdensec(axx,dd1,dd2):
    denmean=dat.pden.sel(date=slice(dd1,dd2)).mean(dim='date').values
    salmat=dat.psal.sel(date=slice(dd1,dd2))
    salmat_sm=salmat.copy()
    for tt,na in enumerate(salmat.date.values):
        for mm,na in enumerate(salmat.distance.values):
                nanind=~isnan(salmat[mm,:,tt])
                Z,X = sig.butter(2,0.02, output='ba')
                salmat_sm[mm,nanind,tt]=sig.filtfilt(Z,X,salmat[mm,nanind,tt].values)
    salmean=salmat_sm.mean(dim='date').values
    sal_levs=[34,34.2,34.4,34.6,34.8,34.85, 34.9,34.91, 34.92,34.93, 34.94,34.95, 34.96,34.97, 34.98,34.99,35]
    ssal=axx.contourf(dat.distance.values,dat.depth.values,salmean.T,sal_levs,cmap=sal_cmap,extend='both')
    dens=axx.contour(dat.distance.values,dat.depth.values,denmean.T,[d1,d2,d3],colors='k',linewidths=3)
    axx.set_ylim(1250,20)
    for ii,dd in enumerate(dat.distance):
        axx.axvline(dd,linewidth=5,color=colvec[ii])
        for zz in [50,250,500,750,1000]:
            axx.plot(dd,zz,'ko')
    axx.set_xlim(44,172)
    return ssal,dens

moorvec=['CF5','CF6','M1','OOI']

def plot_denschem_TS():
    f,axi=subplots(2,2,sharex='col',sharey='col',figsize=(9,6))
    subplots_adjust(wspace=0.4,hspace=0.3)

    f.text(0.5,0.95, 'Before convection (October 2014)', ha='center',fontsize=16)
    f.text(0.5,0.475, 'During convection (February 2015)', ha='center',fontsize=16)

    dd1='2014-10-1'
    dd2='2014-11-1'
    ssal,dens=plotdensec(axi[0,1],dd1,dd2)
    axi[0,1].clabel(dens,fmt='%1.2f',manual=[(125,100),(125,750),(125,1000)])
    mlab=TScomp(axi[0,0],dd1,dd2)
    axi[0,0].clabel(mlab,fmt='%1.2f',manual=[(34.895,4.6),(34.96,4.3),(34.89,3.5)])
    axi[1,1].set_xlabel('distance [km]')
    axi[1,1].set_ylabel('depth [m]')
    axi[0,1].set_ylabel('depth [m]')
    for ii,dd in enumerate(dat.distance):
        axi[0,1].text(dd-5,-5,moorvec[ii],fontsize=13)

    cbax=f.add_axes([0.95,0.25,0.02,0.5])

    dd1='2015-2-1'
    dd2='2015-3-1'
    ssal,dens=plotdensec(axi[1,1],dd1,dd2)
    colorbar(ssal,cax=cbax,label='salinity')
    mlab=TScomp(axi[1,0],dd1,dd2)

    axi[1,0].set_xlabel('salinity')
    axi[1,0].set_ylabel('pot. temperature [$^\circ$C]')
    axi[0,0].set_ylabel('pot. temperature [$^\circ$C]')

    axi[1,0].scatter(0,0,s=20,color='k',label='20m')
    axi[1,0].scatter(0,0,s=200,color='k',label='200m')
    axi[1,0].scatter(0,0,s=50,color=ooicol,edgecolor=ooicol,label='OOI: Irminger gyre interior',alpha=0.8)
    axi[1,0].scatter(0,0,s=50,color=m1col,edgecolor=m1col,label='M1: Offshore of the boundary current',alpha=0.8)
    axi[1,0].scatter(0,0,s=50,color=cf6col,edgecolor=cf6col,label='CF6: Deeper boundary current',alpha=0.8)
    axi[1,0].scatter(0,0,s=50,color=cf5col,edgecolor=cf5col,label='CF5: Boundary current maximum',alpha=0.8)
    axi[1,0].legend(loc=(0,-0.6),ncol=3)

    fs=13

    axi[1,0].text(34.95,4.6,'uISIW',fontsize=fs,weight='bold')
    axi[1,0].text(34.867,3.44,'dISIW',fontsize=fs,weight='bold')

    axi[1,1].text(110,300,'upper ISIW',fontsize=fs+1,weight='bold')
    axi[1,1].text(110,900,'deep ISIW',fontsize=fs+1,weight='bold')
    axi[1,0].set_xticks(arange(34.88,34.98,0.04))

    savefig(figdir+'MixedLayer/paperfigs/Fig2.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig2.pdf',bbox_inches='tight')


plot_denschem_TS()


def plot_denschem_TS_YR2():

    f,axi=subplots(2,2,sharex='col',sharey='col',figsize=(9,6))
    subplots_adjust(wspace=0.4,hspace=0.3)

    f.text(0.5,0.95, 'Before convection (October 2014)', ha='center',fontsize=16)
    f.text(0.5,0.475, 'During convection (February 2015)', ha='center',fontsize=16)

    dd1='2015-10-1'
    dd2='2015-11-1'
    ssal,dens=plotdensec(axi[0,1],dd1,dd2)
    axi[0,1].clabel(dens,fmt='%1.2f',manual=[(125,100),(125,750),(125,1000)])
    mlab=TScomp(axi[0,0],dd1,dd2)
    axi[0,0].clabel(mlab,fmt='%1.2f',manual=[(34.895,4.6),(34.96,4.3),(34.89,3.5)])
    axi[1,1].set_xlabel('distance [km]')
    axi[1,1].set_ylabel('depth [m]')
    axi[0,1].set_ylabel('depth [m]')
    for ii,dd in enumerate(dat.distance):
        axi[0,1].text(dd-5,-5,moorvec[ii],fontsize=13)

    cbax=f.add_axes([0.95,0.25,0.02,0.5])

    dd1='2016-2-1'
    dd2='2016-3-1'
    ssal,dens=plotdensec(axi[1,1],dd1,dd2)
    colorbar(ssal,cax=cbax,label='salinity')
    mlab=TScomp(axi[1,0],dd1,dd2)

    axi[1,0].set_xlabel('salinity')
    axi[1,0].set_ylabel('pot. temperature [$^\circ$C]')
    axi[0,0].set_ylabel('pot. temperature [$^\circ$C]')

    axi[1,0].scatter(0,0,s=20,color='k',label='20m')
    axi[1,0].scatter(0,0,s=200,color='k',label='200m')
    axi[1,0].scatter(0,0,s=50,color=ooicol,edgecolor=ooicol,label='OOI: Irminger gyre interior',alpha=0.6)
    axi[1,0].scatter(0,0,s=50,color=m1col,edgecolor=m1col,label='M1: Offshore of the boundary current',alpha=0.6)
    axi[1,0].scatter(0,0,s=50,color=cf6col,edgecolor=cf6col,label='CF6: Deeper boundary current',alpha=0.3)
    axi[1,0].scatter(0,0,s=50,color=cf5col,edgecolor=cf5col,label='CF5: Boundary current maximum',alpha=0.3)
    axi[1,0].legend(loc=(0,-0.6),ncol=3)

    fs=13

    axi[1,0].text(34.95,4.5,'uISIW',fontsize=fs,weight='bold')
    axi[1,0].text(34.87,3.45,'dISIW',fontsize=fs,weight='bold')

    axi[1,1].text(110,300,'upper ISIW',fontsize=fs+1,weight='bold')
    axi[1,1].text(110,1000,'deep ISIW',fontsize=fs+1,weight='bold')
    axi[1,0].set_xticks(arange(34.88,34.98,0.04))

    savefig(figdir+'MixedLayer/paperfigs/F3_year2.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/F3_year2.pdf',bbox_inches='tight')


plot_denschem_TS_YR2()
