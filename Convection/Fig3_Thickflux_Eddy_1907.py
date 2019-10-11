from aux_funcs import *

dendat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/density_gridded_props_cf5-oom_from5m.nc')

WM={}
WM['upper']=dendat.thickness[:,p1:p2,:].sum(dim='den')
WM['deep']=dendat.thickness[:,p2:p3,:].sum(dim='den')

WM['outcrop_upper']=isnan(dendat.thickness[:,p1,:])
WM['outcrop_deep']=isnan(dendat.thickness[:,p2,:])

### Smooth the thickness time series
Z,X = sig.butter(2,1./10, output='ba')

WM['u_sm']=NaN*WM['upper']
WM['d_sm']=NaN*WM['deep']
WM['u_sm_out']=NaN*WM['outcrop_upper']
WM['d_sm_out']=NaN*WM['outcrop_deep']
for ii in range(len(dendat.distance)):
    WM['u_sm'][ii,:]=sig.filtfilt(Z,X,WM['upper'][ii,:].values)
    WM['d_sm'][ii,:]=sig.filtfilt(Z,X,WM['deep'][ii,:].values)
    WM['u_sm_out'][ii,:]=sig.filtfilt(Z,X,WM['outcrop_upper'][ii,:].values)
    WM['d_sm_out'][ii,:]=sig.filtfilt(Z,X,WM['outcrop_deep'][ii,:].values)


WM['u_sm_out'].plot()

CF5=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/AQD_cf5_500m.nc')
CF6=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/AQD_cf6_500m.nc')
M1=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/AQD_m1_500m.nc')

for xxx in [CF5,CF6,M1]:
    xxx['across track']=xxx['u']*cos(theta)+xxx['v']*sin(theta)
    xxx['along track']=xxx['v']*cos(theta)-xxx['u']*sin(theta)

cos(theta)
sin(theta)

### colors
# cf5col='#F74E4A'
# cf6col='#F7984A'
# ooicol='#3AC243'
# m1col='#2C9494'

ooicol='#E45D2B'
m1col='#E4942B'
cf5col='#236493'
cf6col='#1E9E64'


colvec=[cf5col,cf6col,m1col]

moorvec=['CF5','CF6','M1']


def plot_eddyvar(axx,tit):
    field='u'
    for ii,xx in enumerate([CF5,CF6,M1]):
        dstep=24
        if ii==2:
            elfieldo=xx[field]
        else:
            elfieldo=xx[field][::2]
        Z,X = sig.butter(2,1./dstep/20, output='ba')
        xx_resamp=elfieldo.rolling(date=dstep).std()
        xx_max=xx_resamp.resample(date='1D').mean()
        axx.plot(xx_max.date,(xx_max.values)**2,color=colvec[ii],alpha=0.4)
        xx_sm=sig.filtfilt(Z,X,xx_resamp[~isnan(xx_resamp)].values)
        axx.plot(xx_resamp.date[~isnan(xx_resamp)],xx_sm**2,linewidth=2,color=colvec[ii],zorder=3,label=moorvec[ii])
        axx.set_title(tit,fontsize=14)
        axx.set_ylabel('[m$^2$/s$^2$]')
        leg=axx.legend(loc=(1.05,0.2),fontsize=14)
        for line in leg.get_lines():
            line.set_linewidth(3)
    axx.set_ylim(0,0.02)
    axx.set_yticks(arange(0,0.025,0.01))


def cont_thick(field1,field2,axx,maxt,hlevs,labit):
    hh=axx.contourf(dendat.date.values,dendat.distance.values,field1.values,levels=linspace(0,maxt,50),cmap=cm.viridis,extend='max')
    dd=axx.contour(dendat.date.values,dendat.distance.values,field1.values,levels=hlevs,colors='k',linewidths=2)
    axx.contourf(dendat.date.values,dendat.distance.values,field2.values, [0.75,1.25], colors='grey',alpha=0.5,hatches=['\\\\'],)
    axx.clabel(dd,fmt='%d')
    axx.set_title(labit,fontsize=14)
    colorbar(hh,ax=axx,ticks=linspace(0,maxt,4),pad=0.02,label='layer thickness [m]')
    axx.set_ylabel('distance [km]')
    ax2=axx.twinx()
    if 'Deep' in labit:
        ax2.set_yticks(dendat.distance)
        ax2.set_yticklabels(['CF5','CF6','M1','OOI'])
        ax2.set_ylim(45,170)
    elif 'Upper' in labit:
        ax2.set_yticks(dendat.distance)
        ax2.set_yticklabels(['CF5','CF6','M1','OOI'])
        ax2.set_ylim(45,170)
        axx.set_ylim(45,170)

def plotprs(axx,tit):
    for ii,xx in enumerate([CF5]):
        axx.plot(xx.date,xx.prs,color=colvec[ii])
    axx.set_ylabel('pressure [db]')
    axx.set_title(tit,fontsize=14)

def plot_thickhov_eddy():
    f,axi=subplots(3,1,figsize=(10,6),constrained_layout=True)
    cont_thick(WM['d_sm'],WM['d_sm_out'],axi[1],1500,[500,1000],'b) Deep ISIW: $\sigma_{\Theta}=$ 27.73 - 27.77')
    cont_thick(WM['u_sm'],WM['u_sm_out'],axi[0],750,[500],'a) Upper ISIW: $\sigma_{\Theta}=$ 27.65 - 27.73')
    plot_eddyvar(axi[2],'c) Variance in cross-stream velocity at 500m instrument')
    # plotprs(axi[3],'d) Pressure at 500m instrument on CF5 mooring')

    for axx in axi:
        axx.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        axx.xaxis.set_minor_formatter(monthFMT)
        axx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter(''))
        for yy in [2014,2015]:
            axx.axvline(datetime.datetime(yy,12,1),color='k',linewidth=3)
        for yy in [2015,2016]:
            axx.axvline(datetime.datetime(yy,4,30),color='k',linewidth=3)

    axi[-1].xaxis.set_minor_formatter(monthFMT)
    axi[-1].xaxis.set_major_formatter(yearFMT)
    # axi[-1].set_ylim(0,3000)

    savefig(figdir+'MixedLayer/paperfigs/Fig3.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Fig3.pdf',bbox_inches='tight')

plot_thickhov_eddy()

def plot_thickhov_eddy_nolet():
    f,axi=subplots(3,1,figsize=(10,6),constrained_layout=True)
    cont_thick(WM['d_sm'],WM['d_sm_out'],axi[1],1500,[500,1000],'Deep ISIW: $\sigma_{\Theta}=$ 27.73 - 27.77')
    cont_thick(WM['u_sm'],WM['u_sm_out'],axi[0],750,[500],'Upper ISIW: $\sigma_{\Theta}=$ 27.65 - 27.73')
    plot_eddyvar(axi[2],'Variance in cross-stream velocity at 500m instrument')
    for axx in axi:
        axx.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        axx.xaxis.set_minor_formatter(monthFMT)
        axx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter(''))
        for yy in [2014,2015]:
            axx.axvline(datetime.datetime(yy,12,1),color='k',linewidth=3)
        for yy in [2015,2016]:
            axx.axvline(datetime.datetime(yy,4,30),color='k',linewidth=3)

    axi[-1].xaxis.set_minor_formatter(monthFMT)
    axi[-1].xaxis.set_major_formatter(yearFMT)
    axi[-1].set_yticks(arange(0,0.025,0.01))
    savefig(figdir+'MixedLayer/paperfigs/Noletters_Fig3.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Noletters_Fig3.pdf',bbox_inches='tight')

plot_thickhov_eddy_nolet()


def plot_thickhov_eddy_twopan():
    f,axi=subplots(2,1,figsize=(10,4),constrained_layout=True)
    cont_thick(WM['d_sm'],WM['d_sm_out'],axi[1],1500,[500,1000],'Deep ISIW: $\sigma_{\Theta}=$ 27.73 - 27.77')
    cont_thick(WM['u_sm'],WM['u_sm_out'],axi[0],750,[500],'Upper ISIW: $\sigma_{\Theta}=$ 27.65 - 27.73')
    # plot_eddyvar(axi[2],'Variance in cross-stream velocity at 500m instrument')
    for axx in axi:
        axx.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        axx.xaxis.set_minor_formatter(monthFMT)
        axx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter(''))
        for yy in [2014,2015]:
            axx.axvline(datetime.datetime(yy,12,1),color='k',linewidth=3)
        for yy in [2015,2016]:
            axx.axvline(datetime.datetime(yy,4,30),color='k',linewidth=3)

    axi[-1].xaxis.set_minor_formatter(monthFMT)
    axi[-1].xaxis.set_major_formatter(yearFMT)
    # axi[-1].set_yticks(arange(0,0.025,0.01))
    savefig(figdir+'MixedLayer/paperfigs/Twopan_Fig3.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Twopan_Fig3.pdf',bbox_inches='tight')

plot_thickhov_eddy_twopan()


def plot_thickhov_eddy_onepan():
    f,axi=subplots(1,1,figsize=(10,2.5),constrained_layout=True)
    # cont_thick(WM['d_sm'],WM['d_sm_out'],axi[1],1500,[500,1000],'Deep ISIW: $\sigma_{\Theta}=$ 27.73 - 27.77')
    cont_thick(WM['u_sm'],WM['u_sm_out'],axi,750,[500],'Upper ISIW: $\sigma_{\Theta}=$ 27.65 - 27.73')
    # plot_eddyvar(axi[2],'Variance in cross-stream velocity at 500m instrument')
    axx=axi
    axx.set_xlim([datetime.datetime(2014,9,5),datetime.datetime(2016,7,15)])
    axx.xaxis.set_major_locator(years)
    axx.xaxis.set_minor_locator(threemonth)
    axx.xaxis.set_minor_formatter(monthFMT)
    axx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter(''))
    for yy in [2014,2015]:
        axx.axvline(datetime.datetime(yy,12,1),color='k',linewidth=3)
    for yy in [2015,2016]:
        axx.axvline(datetime.datetime(yy,4,30),color='k',linewidth=3)

    axi.xaxis.set_minor_formatter(monthFMT)
    axi.xaxis.set_major_formatter(yearFMT)
    # axi[-1].set_yticks(arange(0,0.025,0.01))
    savefig(figdir+'MixedLayer/paperfigs/Onepan_Fig3.png',bbox_inches='tight',dpi=300)
    savefig(figdir+'MixedLayer/paperfigs/Onepan_Fig3.pdf',bbox_inches='tight')

plot_thickhov_eddy_onepan()
