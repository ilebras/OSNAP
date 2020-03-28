from aux_funcs import *
daily=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid.nc')
daily['across track velocity']=-1*daily['across track velocity']
daily=daily.where(daily['across track velocity']!=0)
#
daily=daily.sel(date=slice('2014-8-15','2016-7-1')) #just to avoid some zeros at the beginning which mess up filtering.
# #   PLUS HERE IM RESTRICTING DATES TO FIRST DEPLOYMENT(ISH)
daily

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

onesxr=daily.salinity/daily.salinity

srefb=34.9
sep=9

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport over 27.8']=daily['across track velocity'].where(daily['potential density']<27.8)[:,:-1,:]*depthdiffmat*middistmat/1e3
daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3


egic={}
egic['trans']=daily['xport over 27.8'][sep:,:,:].sum('depth').sum('distance')
egic['area']=(onesxr.where(daily['potential density']<27.8)[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
egic['sal']=(daily['xport over 27.8'][sep:,:-1,:]*daily['salinity'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['tmp']=(daily['xport over 27.8'][sep:,:-1,:]*daily['temperature'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['den']=(daily['xport over 27.8'][sep:,:-1,:]*daily['potential density'][sep:,:-1,:]).sum('distance').sum('depth')/egic['trans']
egic['meanvel']=egic['trans']/egic['area']

for ii in range(0,200,20):
    figure()
    (daily['potential density'][sep:,:-1,:].where((daily['potential density']<d2)&(daily['potential density']>=d1)).isel(date=ii).T).plot()
    (daily['potential density'][sep:,:-1,:].where((daily['potential density']<d3)&(daily['potential density']>=d2)).isel(date=ii).T).plot()

uIIW={}
uIIW['trans']=daily.xport.where((daily['potential density']<d2)&(daily['potential density']>=d1)).sum('distance').sum('depth')
uIIW['area']=(onesxr.where((daily['potential density']<d2)&(daily['potential density']>=d1))[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
uIIW['trans cf5+']=daily.xport.where(daily.distance>=45).where((daily['potential density']<d2)&(daily['potential density']>=d1)).sum('distance').sum('depth')
uIIW['meanvel']=uIIW['trans']/uIIW['area']

dIIW={}
dIIW['trans']=daily.xport.where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')
dIIW['area']=(onesxr.where((daily['potential density']<d3)&(daily['potential density']>=d2))[sep:,:-1,:]*depthdiffmat[sep:,:,:]*middistmat[sep:,:,:]/1e3).sum('depth').sum('distance')
dIIW['meanvel']=dIIW['trans']/dIIW['area']
dIIW['trans cf5+']=daily.xport.where(daily.distance>=45).where((daily['potential density']<d3)&(daily['potential density']>=d2)).sum('distance').sum('depth')

IIW={}
IIW['trans']=daily.xport.where((daily['potential density']<d3)&(daily['potential density']>=d1)).sum('distance').sum('depth')

lt={}
lt['trans']=daily.xport.where(daily['potential density']<d1).sum('distance').sum('depth')

mt={}
mt['trans']=daily.xport.where((daily['potential density']>=d3)&(daily['potential density']<27.8)).sum('distance').sum('depth')


N  = 2    # Filter order
Wn = 1./30
B, A = sig.butter(N, Wn, output='ba')

def plot_trans_only(axi,field,col,labit):
    hfa=0.3
    axi.plot(field.date,field,color=col,alpha=hfa,label='')
    axi.plot(field.date,sig.filtfilt(B,A,field),color=col,linewidth=3,alpha=1,label=labit)


def plot_trans(axi,field,col,labit):
    hfa=0.3
    axi.plot(field.date,field,color=col,alpha=hfa,label='')
    axi.plot(field.date,sig.filtfilt(B,A,field),color=col,linewidth=3,alpha=1,label=labit)

def plot_Fig4():
    f,[ax0,ax1]=subplots(2,1, sharex=True,figsize=(9,5.5),)
    plot_trans_only(ax0,egic['trans'],'k','total')
    ax0.set_ylim(10,27)
    # ax0.set_ylabel('Transport [Sv]')
    ax0.set_title('a) Total boundary current transport',fontsize=14)
    # ax1.set_yticks(range(0,31,10))
    plot_trans(ax1,lt['trans'],'darkorange','lighter waters')
    plot_trans(ax1,uIIW['trans'],uppercol,'upper ISIW')
    plot_trans(ax1,dIIW['trans'],deepcol,'deep ISIW')
    plot_trans(ax1,mt['trans'],'brown','denser waters')
    ax1.set_title('b) Transport within each density layer',fontsize=14)
    # ax1.set_ylabel('Ratio of transports')
    ax1.legend(loc=(1.02,0.1),fontsize=13)
    # plot_transdecomp(ax3,uIIW,uppercol)
    # plot_transdecomp(ax4,dIIW,deepcol)
    ax2=ax1
    ax2.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])

    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_ylim(0,13)
    f.text(0.05, 0.5, 'Transport [Sv]', va='center', rotation='vertical',fontsize=13)
    savefig(figdir+'MixedLayer/paperfigs/Fig4.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/Fig4.png',bbox_inches='tight',dpi=300)

plot_Fig4()


def plot_Fig4_ml():
    f,[ax0,ax1]=subplots(2,1, sharex=True,figsize=(9,5.5),)
    plot_trans_only(ax0,egic['trans'],'k','total')
    ax0.set_ylim(10,27)
    # ax0.set_ylabel('Transport [Sv]')
    ax0.set_title('a) Total boundary current transport',fontsize=14)
    # ax1.set_yticks(range(0,31,10))
    plot_trans(ax1,lt['trans'],'darkorange','lighter waters')
    plot_trans(ax1,uIIW['trans'],uppercol,'upper ISIW')
    plot_trans(ax1,dIIW['trans'],deepcol,'deep ISIW')
    plot_trans(ax1,mt['trans'],'brown','denser waters')
    ax1.set_title('b) Transport within each density layer',fontsize=14)
    # ax1.set_ylabel('Ratio of transports')
    ax1.legend(loc=(-0.05,-0.5),fontsize=13,ncol=4)
    # plot_transdecomp(ax3,uIIW,uppercol)
    # plot_transdecomp(ax4,dIIW,deepcol)
    ax2=ax1
    ax2.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])

    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_ylim(0,13)
    f.text(0.05, 0.5, 'Transport [Sv]', va='center', rotation='vertical',fontsize=13)
    savefig(figdir+'MixedLayer/paperfigs/Moveleg_Fig4.pdf',bbox_inches='tight')

plot_Fig4_ml()

def plot_Fig4_nolet():
    f,[ax0,ax1]=subplots(2,1, sharex=True,figsize=(9,5.5),)
    plot_trans_only(ax0,egic['trans'],'k','total')
    ax0.set_ylim(10,27)
    # ax0.set_ylabel('Transport [Sv]')
    ax0.set_title('Total boundary current transport',fontsize=14)
    # ax1.set_yticks(range(0,31,10))
    plot_trans(ax1,lt['trans'],'darkorange','lighter waters')
    plot_trans(ax1,uIIW['trans'],uppercol,'upper ISIW')
    plot_trans(ax1,dIIW['trans'],deepcol,'deep ISIW')
    plot_trans(ax1,mt['trans'],'brown','denser waters')
    ax1.set_title('Transport within each density layer',fontsize=14)
    # ax1.set_ylabel('Ratio of transports')
    ax1.legend(loc=(1.02,0.1),fontsize=13)
    # plot_transdecomp(ax3,uIIW,uppercol)
    # plot_transdecomp(ax4,dIIW,deepcol)
    ax2=ax1
    ax2.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])

    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_ylim(0,13)
    f.text(0.05, 0.5, 'Transport [Sv]', va='center', rotation='vertical',fontsize=13)
    savefig(figdir+'MixedLayer/paperfigs/Noletters_Fig4.pdf',bbox_inches='tight')

plot_Fig4_nolet()


def plot_Fig4_onepan():
    f,ax1=subplots(1,1, sharex=True,figsize=(9,3),)
    # plot_trans_only(ax0,egic['trans'],'k','total')
    # ax0.set_ylim(10,27)
    # # ax0.set_ylabel('Transport [Sv]')
    # ax0.set_title('Total boundary current transport',fontsize=14)
    # # ax1.set_yticks(range(0,31,10))
    plot_trans(ax1,lt['trans'],'darkorange','lighter waters')
    plot_trans(ax1,uIIW['trans'],uppercol,'upper ISIW')
    plot_trans(ax1,dIIW['trans'],deepcol,'deep ISIW')
    plot_trans(ax1,mt['trans'],'brown','denser waters')
    ax1.set_title('Transport within each density layer',fontsize=14)
    # ax1.set_ylabel('Ratio of transports')
    ax1.legend(loc=(1.02,0.1),fontsize=13)
    # plot_transdecomp(ax3,uIIW,uppercol)
    # plot_transdecomp(ax4,dIIW,deepcol)
    ax2=ax1
    ax2.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])

    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.set_ylim(0,13)
    f.text(0.05, 0.5, 'Transport [Sv]', va='center', rotation='vertical',fontsize=13)
    savefig(figdir+'MixedLayer/paperfigs/Onepan_Fig4.pdf',bbox_inches='tight')

plot_Fig4_onepan()

def plot_transdecomp(axi,field,col):
    ls=2
    hfa=1
    N  = 2    # Filter order
    Wn = 1./30
    B, A = sig.butter(N, Wn, output='ba')
    axi.plot(field['trans'].date,sig.filtfilt(B,A,field['trans']),color=col,linewidth=ls+1,alpha=hfa-0.2)
    vbar=mean(field['meanvel'])
    Abar=mean(field['area'])
    print(vbar,Abar)
    vprime=field['meanvel']-vbar
    Aprime=field['area']-Abar

    # axi.axhline((vbar*Abar).values,color=col,linewidth=ls-1,alpha=hfa)
    # axi.plot(field['trans'].date,sig.filtfilt(B,A,vbar*Abar+vbar*Aprime),'-',color=col,linewidth=ls,alpha=hfa)
    # axi.plot(field['trans'].date,sig.filtfilt(B,A,vbar*Abar+vprime*Abar),'--',color=col,linewidth=ls,alpha=hfa)
    # axi.plot(field['trans'].date,sig.filtfilt(B,A,vbar*Abar+vprime*Aprime),'-.',color=col,linewidth=ls-1,alpha=hfa)


plot(egic['trans'].date,sig.filtfilt(B,A,(uIIW['trans']/dIIW['trans'])))
plot(egic['trans'].date,sig.filtfilt(B,A,(dIIW['trans'])))

def plot_Fig4_withdecomp():
    f,[ax1,ax2]=subplots(2,1, sharex=True,figsize=(10,8))
    plot_transdecomp(ax1,uIIW,uppercol)
    plot_transdecomp(ax1,dIIW,deepcol)

    ax2.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)

plot_Fig4_withdecomp()
    #
    # lw=3
    # ax0.plot(psi.TIME,MOC_east,color='k',linewidth=lw)
    # axii=ax0.twinx()
    # axii.plot(osnap.TIME,uIIW['osnap east trans'],color=uppercol,linewidth=lw)
    # axii.plot(osnap.TIME,dIIW['osnap east trans'],color=deepcol,linewidth=lw)
    # ax0.text(datetime.datetime(2016,4,1),8,'upper IIW',color=uppercol,fontsize=16)
    # ax0.text(datetime.datetime(2016,4,1),-5,'deep IIW',color=deepcol,fontsize=16)
    # ax0.set_xticklabels('')
    # ax0.set_ylim(-10,25)
    # ax0.axhline(0,color='k')
    # axii.set_ylim(-5,12.5)
    # ax0.set_ylabel('OSNAP East\nOverturning transport [Sv]')
    # axii.set_ylabel('IIW contribution\nto overturning [Sv]')
    #
    # plot_trans(ax1)
    # spowcomp(ax2)
    #
    # for axx in [ax0,ax1,ax2]:
    #     axx.set_xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,7,15)])
    #     axx.xaxis.set_major_locator(years)
    #     axx.xaxis.set_minor_locator(threemonth)
    #     if (axx==ax2):
    #         axx.xaxis.set_minor_formatter(monthFMT)
    #         axx.xaxis.set_major_formatter(yearFMT)




# versname='1810JHcal'
# grd=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_'+versname+'.pickle','rb'))
grd=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_rotnewfieldnames.nc').sel(date=slice('2014-8-15','2016-7-1'))
grd
dat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_5m.nc')


dendiff=dat.pden.diff(dim='distance')
geoshear=dendiff*9.8/gsw.f(60)/1028/diff(dat.distance)[0].T/1e3
geoshear

geovel=geoshear[:,::-1,:].cumsum(dim='depth')[:,::-1,:]*unique(diff(dat.depth))

#compare geovel trans between CF5+6, and directly measured trans between them -- back out the barotropic comp
#do the same between CF6 and M1
diff(dat.distance)

geotrans={}
geotrans['56']=(geovel[0,:,:].sum(dim='depth')*diff(dat.distance)[0]*5/1e3)
geotrans['61']=(geovel[1,:,:].sum(dim='depth')*diff(dat.distance)[1]*5/1e3)

grd

dat.distance
dirtrans={}
dirtrans['56']=(grd['across track velocity'][4,:,:]+grd['across track velocity'][5,:,:]).sum(dim='depth')*-5*diff(dat.distance)[0]/1e3
dirtrans['61']=(grd['across track velocity'][5,:,:]+grd['across track velocity'][-1,:,:]).sum(dim='depth')*-5*diff(dat.distance)[1]/1e3

grd.depth

(grd['across track velocity'][4,:,:]+grd['across track velocity'][5,:,:]).plot()

(dirtrans['56']+dirtrans['61']).plot()

barotrans={}
for kk in geotrans:
    barotrans[kk]=dirtrans[kk]-geotrans[kk]

def complot():
    figure(figsize=(12,3))
    resamp='1D'
    geotrans['56'].resample(date=resamp).mean(dim='date').plot(label='geo 56')
    geotrans['61'].resample(date=resamp).mean(dim='date').plot(label='geo 61')
    (geotrans['56']+geotrans['61']).resample(date=resamp).mean(dim='date').plot(label='geo tot')
    legend()

    figure(figsize=(12,3))
    barotrans['56'].resample(date=resamp).mean(dim='date').plot(label='baro 56')
    barotrans['61'].resample(date=resamp).mean(dim='date').plot(label='baro 61')
    (barotrans['56']+barotrans['61']).resample(date=resamp).mean(dim='date').plot(label='baro tot')
    legend()


    figure(figsize=(12,3))
    (dirtrans['56']+dirtrans['61']).resample(date=resamp).mean(dim='date').plot(label='from dir')
    (geotrans['56']+geotrans['61']+13).resample(date=resamp).mean(dim='date').plot(label='geo tot + 13')
    egic['trans'].resample(date=resamp).mean(dim='date').plot(label='tot')
    legend()

complot()

figure(figsize=(12,3))
ii=0
for dd in range(200,1100,200):
    dendiff[0,:,:].sel(depth=dd).plot(label='',color='C'+str(ii),alpha=0.2)
    dendiff[0,:,:].sel(depth=dd).resample(date='1M').mean(dim='date').plot(label=dd,linewidth=3,color='C'+str(ii))
    ii+=1
axhline(0,color='k')
legend(loc=(1.05,0))
ylim(-0.05,0.2)
title('')


figure(figsize=(12,3))
ii=0
for dd in range(200,1800,200):
    dendiff[1,:,:].sel(depth=dd).plot(label='',color='C'+str(ii),alpha=0.2)
    dendiff[1,:,:].sel(depth=dd).resample(date='1M').mean(dim='date').plot(label=dd,linewidth=3,color='C'+str(ii))
    ii+=1
axhline(0,color='k')
legend(loc=(1.05,0))
ylim(-0.05,0.2)
title('')


grd['across track velocity'][4:6,:,:].mean(dim='distance').isel(depth=100).plot(figsize=(12,3),alpha=0.2)
grd['across track velocity'][4:6,:,:].mean(dim='distance').isel(depth=100).resample(date='1M').mean(dim='date').plot(linewidth=3)
ylim(-0.2,-0.55)


geovel[0,:,:].sel(depth=200).plot(figsize=(12,3),alpha=0.2)
geovel[0,:,:].sel(depth=200).resample(date='1M').mean(dim='date').plot(linewidth=3)
ylim(0,0.35)



geovel[1,:,:].sel(depth=200).plot(figsize=(12,3),alpha=0.2)
geovel[1,:,:].sel(depth=200).resample(date='1M').mean(dim='date').plot(linewidth=3)
ylim(0,0.35)

figure(figsize=(12,3))
(grd['across track velocity'][4:6,:,:].mean(dim='distance').isel(depth=100)*-1-0.2).resample(date='1M').mean(dim='date').plot(linewidth=3)
geovel.sel(depth=200).resample(date='1M').mean(dim='date').plot(linewidth=3)

figure(figsize=(12,3))
(grd['across track velocity'][4:6,:,:].mean(dim='distance').isel(depth=300)*-1-0.2).resample(date='1M').mean(dim='date').plot(linewidth=3)
geovel.sel(depth=600).resample(date='1M').mean(dim='date').plot(linewidth=3)
geovel

((geovel[0,:,:].sel(depth=slice(0,2000)).sum(dim='depth')*diff(dat.distance)[0]+geovel[1,:,:].sel(depth=slice(0,3000)).sum(dim='depth')*diff(dat.distance)[1])*5/1e3).resample(date='1M').mean(dim='date').plot(figsize=(12,3),linewidth=3)
(geovel[1,:,:].sel(depth=slice(0,3000)).sum(dim='depth')*5*diff(dat.distance)[1]/1e3).resample(date='1M').mean(dim='date').plot(figsize=(12,3),linewidth=3)
(grd['across track velocity'][4:7,:,:].mean(dim='distance').sel(depth=slice(200,2000)).sum(dim='depth')*-2).resample(date='1M').mean(dim='date').plot(figsize=(12,3),linewidth=3)

plot(geovel.sel(date=slice('2014-9-1','2014-10-1')).mean(dim='date'),geovel.depth,label='Sept 2014')
plot(geovel.sel(date=slice('2014-10-1','2014-11-1')).mean(dim='date'),geovel.depth,label='Oct 2014')
plot(geovel.sel(date=slice('2014-11-1','2014-12-1')).mean(dim='date'),geovel.depth,label='Nov 2014')
plot(geovel.sel(date=slice('2014-12-1','2015-1-1')).mean(dim='date'),geovel.depth,label='Dec 2014')
plot(geovel.sel(date=slice('2015-1-1','2015-2-1')).mean(dim='date'),geovel.depth,label='Jan 2015')
plot(geovel.sel(date=slice('2015-2-1','2015-3-1')).mean(dim='date'),geovel.depth,label='Feb 2015')
ylim(1250,150)
legend()

gv_3M=geovel.resample(date='2M').mean(dim='date')

for ii in gv_3M.date[:6]:
    plot(gv_3M.sel(date=ii),geovel.depth,label=str(ii.values)[:10])
ylim(1250,150)
legend()

for ii in gv_3M.date[6:]:
    plot(gv_3M.sel(date=ii),geovel.depth,label=str(ii.values)[:10])
ylim(1250,150)
legend()

for ii in gv_3M.date[:6]:
    plot(gv_3M.sel(date=ii)[::-1].cumsum(dim='depth')[::-1]*10,geovel.depth,label=str(ii.values)[:10])
ylim(1250,150)
legend()

for ii in gv_3M.date[6:]:
    plot(gv_3M.sel(date=ii)[::-1].cumsum(dim='depth')[::-1]*10,geovel.depth,label=str(ii.values)[:10])
ylim(1250,150)
legend()


#look at horizontal density gradient between CF5 and CF6
dat.pden.diff(dim='distance').isel(distance=0)[::-1,:].cumsum(dim='depth').plot(figsize=(12,3))
ylim(1250,0)

for dd in range(100,1000,100):
    dat.sel(depth=dd).pden.sel(date=slice('2014-10-1','2014-11-1')).mean(dim='date').plot(figsize=(12,3),label='October 2014')
    dat.sel(depth=dd).pden.sel(date=slice('2015-2-1','2015-3-1')).mean(dim='date').plot(label='February 2015')
    legend()


def plot_trans(axi):
    hfa=0.5
    N  = 2    # Filter order
    Wn = 0.03333333 # Cutoff frequency (30 days)
    B, A = sig.butter(N, Wn, output='ba')
    axi.plot(egic['trans'].date,-egic['trans'],color='k',alpha=hfa)
    axi.plot(egic['trans'].date,-sig.filtfilt(B,A,egic['trans']),color='k',linewidth=3,alpha=0.7)
    axi.set_ylim(8,-37)
    axi.set_ylabel('Boundary current\ntransport [Sv]')
    axx=axi.twinx()
    axx.plot(uIIW['trans'].date,-uIIW['trans'],color=uppercol,alpha=hfa)
    # axx.plot(osnap.TIME,uIIW['osnap east trans'],'o-',color=uppercol,)
    axx.plot(uIIW['trans'].date,-sig.filtfilt(B,A,uIIW['trans']),color=uppercol,linewidth=3)
    # axx.plot(osnap.TIME,dIIW['osnap east trans'],'o-',color=deepcol,)
    axx.plot(dIIW['trans'].date,-dIIW['trans'],color=deepcol,alpha=hfa)
    axx.plot(dIIW['trans'].date,-sig.filtfilt(B,A,dIIW['trans']),color=deepcol,linewidth=3)
    axx.set_ylabel('IIW transport [Sv]')
    axx.set_ylim(0,-19)



    axx.set_xticklabels('')


# grden=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_dengrid/grden_cf3-ooi.pickle','rb'),)

def plot_Fig3_OLD():
    f=figure(figsize=(10,8))

    gs0=gridspec.GridSpec(3,1, hspace = 0.1)
    ax0=plt.subplot(gs0[0])
    ax1=plt.subplot(gs0[1])
    ax2=plt.subplot(gs0[2])

    lw=3
    ax0.plot(psi.TIME,MOC_east,color='k',linewidth=lw)
    axii=ax0.twinx()
    axii.plot(osnap.TIME,uIIW['osnap east trans'],color=uppercol,linewidth=lw)
    axii.plot(osnap.TIME,dIIW['osnap east trans'],color=deepcol,linewidth=lw)
    ax0.text(datetime.datetime(2016,4,1),8,'upper IIW',color=uppercol,fontsize=16)
    ax0.text(datetime.datetime(2016,4,1),-5,'deep IIW',color=deepcol,fontsize=16)
    ax0.set_xticklabels('')
    ax0.set_ylim(-10,25)
    ax0.axhline(0,color='k')
    axii.set_ylim(-5,12.5)
    ax0.set_ylabel('OSNAP East\nOverturning transport [Sv]')
    axii.set_ylabel('IIW contribution\nto overturning [Sv]')

    plot_trans(ax1)
    spowcomp(ax2)

    for axx in [ax0,ax1,ax2]:
        axx.set_xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        if (axx==ax2):
            axx.xaxis.set_minor_formatter(monthFMT)
            axx.xaxis.set_major_formatter(yearFMT)

plot_Fig3_OLD()


    # for axx in [ax0,ax1,ax2]:
    #     for yy in [2015,2016]:
    #         axx.axvspan(datetime.datetime(yy-1,11,1),datetime.datetime(yy,5,1),color='C0',alpha=0.2)



plot_Fig3_OLD()



def plot_eddyvar_plustrans_nosec():
    f=figure(figsize=(10,5))

    gs0=gridspec.GridSpec(2,1, hspace = 0.1)
    ax1=plt.subplot(gs0[0])
    ax2=plt.subplot(gs0[1])

    plot_trans(ax1)
    spowcomp(ax2)

    for axx in [ax1,ax2]:
        axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        if (axx==ax2):
            axx.xaxis.set_minor_formatter(monthFMT)
            axx.xaxis.set_major_formatter(yearFMT)

    for axx in [ax1,ax2]:
        for yy in [2015,2016]:
            axx.axvspan(datetime.datetime(yy-1,11,1),datetime.datetime(yy,5,1),color='C0',alpha=0.2)

    savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans_nosec.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans_nosec.png',bbox_inches='tight')


plot_eddyvar_plustrans_nosec()

XXXXXXXXXXXXXXXXX

#
# cfgrid=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))
#
#
# def plotdensec(dl,axx):
#     dmin=datetime.datetime.strptime(dl, '%Y-%m-%d')-datetime.timedelta(days=15)
#     dmax=datetime.datetime.strptime(dl, '%Y-%m-%d')+datetime.timedelta(days=15)
#     velslice=cfgrid['across track velocity'].sel(date=slice(dmin,dmax))
#     vel=axx.contourf(cfgrid.distance,cfgrid.depth,velslice.mean(dim='date').T,13,vmin=-0.6,vmax=0,cmap=cm.Blues_r,extend='both')
#     # axx.contour(cfgrid.distance,cfgrid.depth,velslice.mean(dim='date').T,arange(-0.6,0.1,0.1),colors='w')#arange(-0.6,-0.1,0.1),colors='C0',linewidths=3)
#     axx.contour(cfgrid.distance,cfgrid.depth,cfgrid['potential density'].sel(date=slice(dmin,dmax)).mean(dim='date').T,levels=dbnds,colors='k',linewidths=4)
#     axx.contour(cfgrid.distance,cfgrid.depth,cfgrid['potential density'].sel(date=slice(dmin,dmax)).mean(dim='date').T,levels=arange(27.61,d3,0.015),colors='k')
#     axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
#     axx.axvline(distvec[4],color=cf5col,linewidth=4)
#     axx.axvline(distvec[5],color=cf6col,linewidth=4)
#     axx.axvline(distvec[7],color=m1col,linewidth=4)
#     # axx.axvline(oom_dist,color=ooicol,linewidth=4)
#     axx.set_ylim(1500,0)
#     axx.set_xlim(30,100)
#     # axx.set_xticks(arange(-42.5,-41,0.5))
#     axx.set_yticks(arange(0,1600,500))
#     return velslice.date.values,vel
#
# def plot_eddyvar_plustrans():
#     f=figure(figsize=(12,10))
#     outer = gridspec.GridSpec(2, 1, height_ratios = [2,1],hspace=0.35)
#     gs0=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[0], hspace = 0.1)
#     ax1=plt.subplot(gs0[0])
#     ax2=plt.subplot(gs0[1])
#     #make nested gridspec at bottom for examples of vel/den sections when eddy activity high
#     gs1 = gridspec.GridSpecFromSubplotSpec(1, 11, subplot_spec = outer[1])
#     ax3=plt.subplot(gs1[1:5])
#     ax4=plt.subplot(gs1[6:10])
#
#     plot_trans(ax1)
#     spowcomp(ax2)
#
#     for axx in [ax1,ax2]:
#         axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
#         axx.xaxis.set_major_locator(years)
#         axx.xaxis.set_minor_locator(threemonth)
#         if (axx==ax2):
#             axx.xaxis.set_minor_formatter(monthFMT)
#             axx.xaxis.set_major_formatter(yearFMT)
#
#     date1='2015-8-15'
#     date2='2016-1-15'
#     ov1,v1=plotdensec(date1,ax3)
#     ov2,v2=plotdensec(date2,ax4)
#     cax1=f.add_axes([0.85,0.125,0.02,0.2])
#     colorbar(v2,ax=ax4,cax=cax1,ticks=arange(-0.6,0,0.2),label='velocity [m/s]')
#     ax4.set_yticklabels('')
#     ax3.set_title('August 2015',fontsize=15)
#     ax4.set_title('January 2016',fontsize=15)
#     ax3.set_ylabel('depth [m]')
#     f.text(0.4, 0.05, 'distance [km]', va='center', rotation='horizontal',fontsize=14)
#     for ii in [ax1,ax2]:
#         for vv in [date1,date2]:
#             [ii.axvline(vv,color='b',linewidth=4,alpha=0.4)]
#     savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans.pdf',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans.png',bbox_inches='tight')
#
#
# plot_eddyvar_plustrans()

#
#
# def plotdensec(dl,axx):
#     dmin=datetime.datetime.strptime(dl, '%Y-%m-%d')-datetime.timedelta(days=20)
#     dmax=datetime.datetime.strptime(dl, '%Y-%m-%d')+datetime.timedelta(days=20)
#     velslice=osnap.VELO.sel(TIME=slice(dmin,dmax))
#     vel=axx.contourf(osnap.LONGITUDE,osnap.DEPTH,velslice.mean(dim='TIME'),31,vmin=-0.5,vmax=0.5,cmap=cm.RdBu_r,extend='both')#arange(-0.6,-0.1,0.1),colors='C0',linewidths=3)
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.sel(TIME=slice(dmin,dmax)).mean(dim='TIME'),levels=dbnds,colors='k',linewidths=4)
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.sel(TIME=slice(dmin,dmax)).mean(dim='TIME'),levels=arange(27.61,d3,0.015),colors='k')
#     axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),2500*ones(len(osnap_bathy['bathy'].flatten())),color='k',zorder=22)
#     axx.axvline(CFlon[4],color=cf5col,linewidth=4)
#     axx.axvline(CFlon[5],color=cf6col,linewidth=4)
#     axx.axvline(CFlon[7],color=m1col,linewidth=4)
#     axx.axvline(ooi_lon['fla'],color=ooicol,linewidth=4)
#     axx.set_ylim(1500,0)
#     axx.set_xlim(-42.5,-41)
#     axx.set_xticks(arange(-42.5,-41,0.5))
#     axx.set_yticks(arange(0,1600,500))
#     return velslice.TIME.values,vel
