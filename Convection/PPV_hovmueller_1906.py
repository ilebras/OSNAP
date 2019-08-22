
###############################################################################################################################
###############################################################################################################################
############################################## PV timeseries with MLDs   ####################################################
###############################################################################################################################
###############################################################################################################################
from aux_funcs import *

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))
ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb'))
oom=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_denmerged_xray.pickle','rb'))

months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')

d1=27.67
d2=27.73
d3=27.77

def plotPV(moornum,axx,tit):
    hh=axx.pcolor(dat.date.resample(date='W'),dat.depth[::10],log10(dat.PV[moornum,::10,:].resample(date='W').mean(dim='date')),cmap=cm.RdYlBu_r,vmin=-11.5,vmax=-10)
    axx.plot(dat.date,dat.ML[:,moornum],'.',color='k',markersize=9)
    axx.plot(dat.date,dat.ML[:,moornum],'.',color='yellow',markersize=5,zorder=101)
    if moornum==4:
        axx.contour(dat.date.resample(date='W').values,dat.depth.values,dat['potential density'][moornum,:,:].resample(date='W').mean(dim='date').values,levels=[d1,d2],colors='k',zorder=100)
    else:
        axx.contour(dat.date.resample(date='W').values,dat.depth.values,dat['potential density'][moornum,:,:].resample(date='W').mean(dim='date').values,levels=[d1,d2,d3],colors='k',zorder=100)
    axx.set_title(tit)
    return hh


def plotPV_ooi(axx,tit):
    axx.pcolor(ooi.date.sel(date=slice('2014-9-1','2015-8-14')),ooi.depth,log10(ooi.PV.sel(date=slice('2014-9-1','2015-8-14'))),cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)
    axx.pcolor(ooi.date.sel(date=slice('2015-8-15','2016-8-1')),ooi.depth,log10(ooi.PV.sel(date=slice('2015-8-15','2016-8-1'))),cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)
    axx.contour(oom.date.values,oom.depth.values,oom['pden'].values.T,levels=[d1,d2,d3],colors='grey',zorder=99)
    axx.contour(ooi.date.values,ooi.depth.values,ooi['potential density'].values,levels=[d1,d2,d3],colors='k',zorder=100)
    # axx.contourf(ooi.date.values,ooi.depth.values,ooi['salinity'].values,levels=arange(34.9,35.05,0.01),zorder=80,cmap=cm.rainbow)

    axx.set_title(tit)
    # axx.plot(loco_date,-loco_mld,color='white',linewidth=5)
    # axx.plot(loco_date,-loco_mld,color='k',linewidth=3)
    axx.plot(loco_date,-loco_mld,'.',color='k',markersize=9,zorder=101)
    axx.plot(loco_date,-loco_mld,'.',color='yellow',markersize=5,zorder=101)
    # axx.plot(ooi_max_ml.date,ooi_max_ml,'k.'

def PV_presplot():
    f,(ax1,ax2,ax3)=subplots(3,1,figsize=(11,7),sharex=True)
    plotPV_ooi(ax1,'Irminger gyre interior, OOI')
    hh=plotPV(7,ax2,'Offshore of the boundary current, M1')
    hh=plotPV(4,ax3,'Boundary current maximum, CF5')
    cbar_ax = f.add_axes([0.93, 0.25, 0.02, 0.6])
    cbar=colorbar(hh, cax=cbar_ax,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]')
    ax1.set_ylim(1500,0)
    ax2.set_ylim(1500,0)
    ax1.set_yticks([1500,750,0])
    ax2.set_yticks([1500,750,0])
    ax3.set_ylim(1300,0)
    f.text(-0.0025, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=12)
    ax2.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,20)])
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    savefig(figdir+'MixedLayer/paperfigs/PPVcomp.png',bbox_inches='tight',dpi=300)


PV_presplot()
