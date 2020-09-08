###############################################################################################################################
###############################################################################################################################
############################################## PV timeseries with MLDs   ####################################################
###############################################################################################################################
###############################################################################################################################
from firstfuncs_1618 import *

dat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/CFall_gridded_bymoor_2m.nc')
dat

dat=dat.transpose('distance','depth','date','lat','lon')
###############################################################################################################################
############################################## Calc and add PV ####################################################
###############################################################################################################################
PVmat=NaN*ones((len(dat.distance),len(dat.depth)-1,len(dat.date)))
for ii in range(len(dat.distance)):
    PVmat[ii,:,:]=sw.f(dat.lat[ii])*dat['PDEN'][ii,:,:].diff(dim='depth')/mean(diff(dat.depth))/(dat['PDEN'].values[ii,:-1,:]+dat['PDEN'][ii,:,:].diff(dim='depth').values/2+1e3)

dat=dat.assign_coords(mid_depth=(dat.depth.values[:-1]+diff(dat.depth.values)/2))
# PVmat[PVmat<=1e-15]=1e-15
dat['PV']=(('distance','mid_depth','date'),PVmat)
dat.PV.sel(distance=distvec[4]).plot()

dat=dat.rename({'PDEN':'pden'})
###############################################################################################################################
############################################## Quick MLD calc (as in Convection work)  ##########################################
###############################################################################################################################
denthresh=0.01
MLmat=NaN*zeros((len(dat.distance),len(dat.date)))
for ii in range(len(dat.distance)):
    for jj in range(len(dat.date)):
        mind=150
        denvec=dat.pden[ii,dat.depth>=mind,jj]
        if ~isnan(denvec[0]):
            diff_from_top=denvec-denvec[0]
            if max(diff_from_top)>denthresh:
                MLmat[ii,jj]=dat.depth[dat.depth>=mind][diff_from_top>denthresh][0]

for ii in range(4,8):
    figure()
    plot(dat.date,MLmat[ii,:].T);
    title('CF'+str(ii+1))
dat['ML']=(('distance','date'),MLmat)

########################################################################################
#### smooth out the PV and density fields for plotting purposes
########################################################################################
dat['PV_sm']=NaN*dat['PV']

dat['pden_sm']=NaN*dat['pden']

### smooth out in depth first
for tt,na in enumerate(dat.date.values):
    for mm,na in enumerate(dat.distance.values):
        nanind=~isnan(dat['PV'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.1, output='ba')
            dat['PV_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['PV'][mm,nanind,tt].values)
        else:
            dat['PV_sm'][mm,:,tt]=NaN
        nanind=~isnan(dat['pden'][mm,:,tt])
        if sum(nanind)>10:
            Z,X = sig.butter(2,0.1, output='ba')
            dat['pden_sm'][mm,nanind,tt]=sig.filtfilt(Z,X,dat['pden'][mm,nanind,tt].values)
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

MLcol='#41ab5d' #'#10b2dd'#
PVcbar=cm.plasma_r

ML_weekly=dat.ML.resample(date='3D').mean(dim='date')
dencol='k'



def plotPV(moornum,axx,tit,moor='na'):
    hh=axx.contourf(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),51,vmin=-11.25,vmax=-10,cmap=PVcbar)
    colorbar(hh,label='log$_{10}$(PPV) [m$^{-1}$ s$^{-1}$]')
    axx.contour(dat.date.values,dat.mid_depth.values,log10(dat.PV_sm[moornum,:,:].values),levels=[-11],colors='red',linewidths=3)
    axx.plot(ML_weekly.date,ML_weekly[moornum,:].T,color=MLcol,linewidth=3)
    # axx.plot(dat.date,dat.ML_[moornum,:].T,color=MLcol,linewidth=3)
    if moornum==4:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2],colors=dencol,zorder=100,linewidths=3)
    else:
        denlab=axx.contour(dat.date.values,dat.depth.values,dat['pden_sm'][moornum,:,:].values,levels=[d1,d2,d3],colors=dencol,zorder=100,linewidths=3)
    # clabel(denlab)
    axx.set_title(tit,fontsize=14)
    axx.set_ylim(1500,150)
    axx.set_ylabel('depth [m]')

figdir

def plotallmoors():
    for mm in range(4,8):
        f,axi=subplots(1,1,figsize=(15,3))
        plotPV(mm,axi,'CF'+str(mm+1))
        savefig(figdir+'ReAnalyze/convection/PPV_hovmueller_CF'+str(mm+1)+'_2m.png',bbox_inches='tight')


d1=27.65
d2=27.73
d3=27.77

plotallmoors()
