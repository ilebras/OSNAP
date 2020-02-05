#################################################################################
#################################################################################
#################################################################################
######################## CALCULATE TRANSPORT  ####################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

daily=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_finergrid.nc')
daily['across track velocity']=-1*daily['across track velocity']
daily=daily.where(daily['across track velocity']!=0)

daily=daily.sel(date=slice('2014-8-15','2016-11-1')) #just to avoid some zeros at the beginning which mess up filtering.
#   PLUS HERE IM RESTRICTING DATES TO FIRST DEPLOYMENT(ISH)

50*1e3/(3*30*24*60**2)

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


##############################################################################################################################
##############################################################################################################################
#################################### Add the PV criterion
##############################################################################################################################
##############################################################################################################################


dat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/gridded_props_cf5-oom_5m_wML.nc')

#### smooth out the PV and density fields for plotting purposes
dat['PV_sm']=NaN*dat['PV']
dat['pden_sm']=NaN*dat['pden']


### smooth out in depth first
for tt,na in enumerate(dat.date.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['PV'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.02, output='ba')
            dat['PV_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['PV'][mm,nanind,tt].values)
        else:
            dat['PV_sm'][mm,:,tt]=NaN
        nanind=~isnan(dat['mden'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.02, output='ba')
            dat['pden_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['mden'][mm,nanind,tt].values)
        else:
            dat['pden_sm'][mm,:,tt]=NaN

## then smooth out in time
for dd,na in enumerate(dat.depth.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['pden_sm'][mm,dd,:])
        if sum(nanind)>10:
            B, A = sig.butter(2,0.05, output='ba')
            dat['pden_sm'][mm,dd,nanind]=sig.filtfilt(B,A,dat['pden_sm'][mm,dd,nanind].values)
        else:
            dat['pden_sm'][mm,dd,:]=NaN

for dd,na in enumerate(dat.mid_depth.values):
    for mm,na in enumerate(dat.distance.values):
            nanind=~isnan(dat['PV_sm'][mm,dd,:])
            if sum(nanind)>10:
                B, A = sig.butter(2,0.1, output='ba')
                dat['PV_sm'][mm,dd,nanind]=sig.filtfilt(B,A,dat['PV_sm'][mm,dd,nanind].values)
                # dat['PV_sm'][mm,dd,nanind]=dat['PV_sm'][mm,dd,nanind].values
            else:
                dat['PV_sm'][mm,dd,:]=NaN

dat

#### Want to use PV criterion to identify newly ventilated waters
#### Specifically, I want to use the smoothed PV and the red line I show on Fig 1 as the threshold (10^-11) (log10(PV)<=-11)
#### Challenge is that PV_sm is on a different grid.
#### First pass approach: Calculate the area of waters that satisfy that criterion within BC, and just use meanvel with this new area to get the new transport.
#### Do a first test to make sure the areas make sense; calculate full area from spareser product, and compare

datones=dat.psal/dat.psal

datones

logPVlim=-11

depthdiff=int(unique(diff(dat.depth)))
distdiff=diff(dat.distance)[:2]
dat_thick_upper=(datones.where((dat['pden']<d2)&(dat['pden']>=d1))*depthdiff).sum('depth')
dat_area_upper=(((dat_thick_upper[:1,:].values+dat_thick_upper[1:2,:].values)/2*distdiff[0]+(dat_thick_upper[1:2,:].values+dat_thick_upper[2:3,:].values)/2*distdiff[1])/1e3).flatten()

dat_thick_deep=(datones.where((dat['pden']<d3)&(dat['pden']>=d2))*depthdiff).sum('depth')
dat_area_deep=(((dat_thick_deep[:1,:].values+dat_thick_deep[1:2,:].values)/2*distdiff[0]+(dat_thick_deep[1:2,:].values+dat_thick_deep[2:3,:].values)/2*distdiff[1])/1e3).flatten()

dat['PV_interp']=(['distance','depth','date'],dat.PV_sm.interp(mid_depth=dat.depth.values).values)
dat_thick_upper_PV=(datones.where((dat['pden']<d2)&(dat['pden']>=d1)&(log10(dat['PV_interp'])<logPVlim))*depthdiff).sum('depth')

dat_area_upper_PV=(((dat_thick_upper_PV[:1,:].values+dat_thick_upper_PV[1:2,:].values)/2*distdiff[0]+(dat_thick_upper_PV[1:2,:].values+dat_thick_upper_PV[2:3,:].values)/2*distdiff[1])/1e3).flatten()

dat_thick_deep_PV=(datones.where((dat['pden']<d3)&(dat['pden']>=d2)&(log10(dat['PV_interp'])<logPVlim))*depthdiff).sum('depth')
dat_area_deep_PV=(((dat_thick_deep_PV[:1,:].values+dat_thick_deep_PV[1:2,:].values)/2*distdiff[0]+(dat_thick_deep_PV[1:2,:].values+dat_thick_deep_PV[2:3,:].values)/2*distdiff[1])/1e3).flatten()
#
# def compareas():
#     figure(figsize=(15,5))
#     uIIW['area'].plot()
#     plot(dat.date,dat_area_upper)
#     plot(dat.date,dat_area_upper_PV)
#     ylabel('upper ISIW area (x $10^6 m^3$)')
#
#     figure(figsize=(15,5))
#     dIIW['area'].plot()
#     plot(dat.date,dat_area_deep)
#     plot(dat.date,dat_area_deep_PV)
#     ylabel('deep ISIW area (x $10^6 m^3$)')
#
# compareas()

dat['upper_area_PV']=(['date'],dat_area_upper_PV)
dat['deep_area_PV']=(['date'],dat_area_deep_PV)

uIIW['trans_PV']=uIIW['meanvel']*dat['upper_area_PV']
dIIW['trans_PV']=(dIIW['meanvel']*dat['deep_area_PV'])

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
    f=figure(figsize=(9,8))
    gs0 = gridspec.GridSpec(3, 1, height_ratios = [1.5,1.5,1],hspace=0.3)
    ax0=plt.subplot(gs0[0])
    ax1=plt.subplot(gs0[1])
    ax2=plt.subplot(gs0[2])
    # f,[ax0,ax1,ax2]=subplots(3,1, sharex=True,figsize=(9.5,8),)
    plot_trans_only(ax0,egic['trans'],'k','total')
    ax0.set_ylim(12,26)
    # ax0.set_ylabel('Transport [Sv]')
    ax0.set_title('a) Total boundary current transport',fontsize=14)
    # ax1.set_yticks(range(0,31,10))
    plot_trans(ax1,lt['trans'],'darkorange','lighter waters')
    plot_trans(ax1,uIIW['trans'],uppercol,'upper ISIW')
    plot_trans(ax1,dIIW['trans'],deepcol,'deep ISIW')
    plot_trans(ax1,mt['trans'],'brown','denser waters')
    ax1.set_title('b) Transport within each density layer',fontsize=14)
    ax1.set_ylim(0,14)
    # ax1.legend(loc=(1.02,-0.3),fontsize=13)
    plot_trans(ax2,uIIW['trans_PV'],uppercol,'upper ISIW')
    plot_trans(ax2,dIIW['trans_PV'],deepcol,'deep ISIW')
    ax2.set_title('c) Boundary current transport of low PPV ISIW',fontsize=14)
    ax2.set_ylim(0,7)
    for axx in [ax0,ax1,ax2]:
        axx.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(threemonth)

    ax2.set_xticklabels(['October','January\n2015','April','July','October','January\n2016','April','July',])
    ax2.set_yticks(arange(0,7,2.5))
    ax2.set_yticklabels(['0','2.5','5'])
    ax1.set_xticklabels([''])
    ax0.set_xticklabels([''])
    ax1.legend(loc=(-0.05,-1.5),fontsize=13,ncol=4)
    f.text(0.05, 0.5, 'Transport [Sv]', va='center', rotation='vertical',fontsize=14)
    savefig(figdir+'MixedLayer/paperfigs/Fig4.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/Fig4.png',bbox_inches='tight',dpi=300)

plot_Fig4()

def botpanel():
    fig,axx=subplots(1,1,figsize=(9,3))
    plot_trans(axx,uIIW['trans_PV'],uppercol,'upper ISIW')
    plot_trans(axx,dIIW['trans_PV'],deepcol,'deep ISIW')
    axx.legend()
    axx.set_ylabel('Transport [Sv]')

    axx.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
    axx.xaxis.set_major_locator(threemonth)

    axx.set_xticklabels(['October','January\n2015','April','July','October','January\n2016','April','July',])
    axx.set_yticks(arange(0,7,2.5))
    savefig(figdir+'MixedLayer/paperfigs/Botpanel_Fig4.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/Botpanel_Fig4.png',bbox_inches='tight',dpi=300)

botpanel()
