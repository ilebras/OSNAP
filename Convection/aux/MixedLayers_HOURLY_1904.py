
27.49+9.73+10.03+11.30+10.71+45.08+4.51+5.95+37.18+17.85+10.31+12.87

#############################################################
########RE-doing this calculation with hourly data  ###############
#############################################################
from aux_funcs import *
datadir=datadir+'OSNAP2016recovery/'

dat=pickle.load(open(datadir+'/pickles/xarray/CF_M_2014-2016_hourlyTSD_1904_noCF7500.pickle','rb'))
# #############################################################
# ############## Mixed layer calc  ######################
# #############################################################
# def getMLDplus():
#     MLD={}
#     MLT={}
#     MLS={}
#     MLP={}
#     ML_psal={}
#     ML_ptmp={}
#     ML_sig0={}
#     dendiff={}
#     PV={}
#     denthresh=0.0075
#     for jj in range(5,9):
#         MLD[jj]=NaN*ones(len(dat[jj].date))
#         MLS[jj]=NaN*ones(len(dat[jj].date))
#         MLT[jj]=NaN*ones(len(dat[jj].date))
#         MLP[jj]=NaN*ones(len(dat[jj].date))
#         ML_psal[jj]=NaN*ones(len(dat[jj].date))
#         ML_ptmp[jj]=NaN*ones(len(dat[jj].date))
#         ML_sig0[jj]=NaN*ones(len(dat[jj].date))
#         PV[jj]=NaN*ones((len(arange(0,nanmin(dat[jj].depth),-2)),len(dat[jj].date)))
#         dendiff[jj]=nancumsum(dat[jj]['sigma0'][::-1,:].diff(dim='dpvec'),axis=0)[::-1,:]
#         for ii in range(len(dat[jj].date)):
#             middepth=dat[jj].depth[:-1,ii]+diff(dat[jj].depth[:,ii])/2
#             f=interp1d(dendiff[jj][:,ii],middepth,fill_value=NaN,bounds_error=False)
#             MLD[jj][ii]=f(denthresh)
#             fsig=interp1d(dat[jj].depth[:,ii],dat[jj].sigma0[:,ii],fill_value=NaN,bounds_error=False)
#             depthvec=arange(0,nanmin(dat[jj].depth),-2)
#             DenInterp=fsig(depthvec)
#             PV[jj][:,ii]=-sw.f(CFlat[jj-1])/(DenInterp+1e3)*gradient(DenInterp)/gradient(arange(0,nanmin(dat[jj].depth),-2))
#             if ~isnan(MLD[jj][ii]):
#                 fsal=interp1d(dat[jj].depth[:,ii],dat[jj].psal[:,ii],fill_value=NaN,bounds_error=False)
#                 MLS[jj][ii]=nanstd(fsal(range(0,int(MLD[jj][ii]),-1)))
#                 ML_psal[jj][ii]=nanmean(fsal(range(0,int(MLD[jj][ii]),-1)))
#                 ftmp=interp1d(dat[jj].depth[:,ii],dat[jj].ptmp[:,ii],fill_value=NaN,bounds_error=False)
#                 MLT[jj][ii]=nanstd(ftmp(range(0,int(MLD[jj][ii]),-1)))
#                 ML_ptmp[jj][ii]=nanmean(ftmp(range(0,int(MLD[jj][ii]),-1)))
#                 MLP[jj][ii]=nanstd(fsig(range(0,int(MLD[jj][ii]),-1)))
#                 ML_sig0[jj][ii]=nanmean(fsig(range(0,int(MLD[jj][ii]),-1)))
#
#     return MLD, MLS, MLT, MLP, PV, dendiff, ML_psal,ML_ptmp,ML_sig0
#
#
# MLD, MLS, MLT, MLP, PV, dendiff, ML_psal,ML_ptmp,ML_sig0 =getMLDplus()
#
# #mask any mixed layers with standard deviations above
# MLD_masked={}
# for ii in MLD:
#     MLD_masked[ii]=MLD[ii].copy()
#     masking=(MLS[ii]>0.005) | (MLT[ii]>0.05) | (MLP[ii] >0.05)
#     MLD_masked[ii][masking]=NaN
#
# dat[5].dpvec
#
#
#
# date=dat[8].date
#
# for ii in MLD:
#     figure()
#     plot(date, MLD[ii],'k.')
#     plot(date, MLD_masked[ii],'r.')
#     title('CF'+str(ii))
#
# newdat={}
# for ii in range(5,9):
#     newXR = xr.Dataset({'ML': (('date'), MLD_masked[ii]),'ML_precorr': (('date'), MLD[ii]),
#                         'ML_psal': (('date'), ML_psal[ii]),'ML_ptmp': (('date'), ML_ptmp[ii]),'ML_sig0': (('date'), ML_sig0[ii]),
#                         'PV': (('depth_reg','date'), -PV[ii]), 'date': date.values,'depth_reg': arange(0,-len(PV[ii])*2,-2)})
#     newdat[ii]=xr.merge([dat[ii], newXR])
#
# pickle.dump(newdat,open(datadir+'pickles/xarray/CF_hrly_wMLPV_1904.pickle','wb'))
#
#
# dat=pickle.load(open(datadir+'pickles/xarray/CF_hrly_wMLPV_1904.pickle','rb'))

###########################################
##### First, compare with OOI TS #####
##########################################
[FLA,FLB]=pickle.load(open(datadir+'pickles/OOI/FLMAB_xrays.pickle','rb'))
FLA_daily=FLA.resample(date='D').mean()
FLB_daily=FLB.resample(date='D').mean()

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat2=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3


def TS_BA(var,axx):
    if var=='before':
        d1='2014-9-15'
        d2='2014-10-15'
    elif var=='after':
        d1='2015-4-15'
        d2='2015-5-15'

    axx.plot(FLA_daily['salinity'].sel(date=slice(d1,d2))[0,0],FLA_daily['temperature'].sel(date=slice(d1,d2))[0,0],'o',color='purple',label='Irminger gyre: OOI')
    axx.plot(FLA_daily['salinity'].sel(date=slice(d1,d2)),FLA_daily['temperature'].sel(date=slice(d1,d2)),'o',color='purple',label='')
    axx.plot(FLB_daily['salinity'].sel(date=slice(d1,d2)),FLB_daily['temperature'].sel(date=slice(d1,d2)),'o',color='purple',label='')
    mm=8
    axx.plot(dat[mm]['psal'].sel(date=slice(d1,d2))[0,0],dat[mm]['ptmp'].sel(date=slice(d1,d2))[0,0],'o',color='C1',label='Just offshore')
    axx.plot(dat[mm]['psal'][:,:].sel(date=slice(d1,d2)),dat[mm]['ptmp'][:,:].sel(date=slice(d1,d2)),'o',color='C1',alpha=0.5,label='')
    mm=5
    axx.plot(dat[mm]['psal'].sel(date=slice(d1,d2))[0,0],dat[mm]['ptmp'].sel(date=slice(d1,d2))[0,0],'o',color='C0',label='Slope Current')
    axx.plot(dat[mm]['psal'][:,:].sel(date=slice(d1,d2)),dat[mm]['ptmp'][:,:].sel(date=slice(d1,d2)),'o',color='C0',alpha=0.1,label='')
    den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.68,27.74,27.8])

def TSBA_pres():
    f,(ax1,ax2)=subplots(1,2,figsize=(10,3),sharex=True,sharey=True)
    TS_BA('before',ax1)
    TS_BA('after',ax2)
    ax1.text(34.98,4.7,'upper IIW',fontsize=12)
    ax1.text(34.98,4.1,'deep IIW',fontsize=12)
    ax1.set_ylim(2.5,6.5)
    ax1.set_xlim(34.85,35.02)
    ax1.set_xticks(arange(34.85,35.02,0.05))
    ax1.set_title('October 2014 (before convection)')
    ax2.set_title('May 2015 (after convection)')
    ax1.legend(loc=(0.3,-0.4),ncol=3)
    ax1.set_ylabel('Temperature ($^\circ$C)')
    f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=12)
    savefig(figdir+'MixedLayer/pres/TS_BA_hrly.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/pres/TS_BA_hrly.png',bbox_inches='tight',dpi=300)


TSBA_pres()

help(hexbin)

# Below was not successful, instead I should make TS ranges and means in density space, use lines and v.light transparent shading. going to try this with daily data first.
# ### Instead, show each of these volumetrically before and after
# def TSvol(salvar,tmpvar,axx,d1,d2):
#     smin=34.85
#     smax=35.02
#     tmin=2.5
#     tmax=6.5
#     vv=axx.hexbin(salvar.where((salvar>smin)&(salvar<smax)).where((tmpvar>tmin)&(tmpvar<tmax)).sel(date=slice(d1,d2)).values.flatten(),
#                tmpvar.where((salvar>smin)&(salvar<smax)).where((tmpvar>tmin)&(tmpvar<tmax)).sel(date=slice(d1,d2)).values.flatten()
#                ,cmap=cm.rainbow,mincnt=1,rasterized='True',gridsize=100,vmin=1,vmax=10)
#     # colorbar(vv)
#     den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.68,27.74,27.8])
#
# def TSvolpans():
#
#     f,axx=subplots(3,2,sharex=True,sharey=True,figsize=(12,10))
#     #Before
#     b1='2014-9-15'
#     b2='2014-10-15'
#     TSvol(FLA_daily['salinity'],FLA_daily['temperature'],axx[0,0],b1,b2)
#     TSvol(dat[8]['psal'],dat[8]['ptmp'],axx[1,0],b1,b2)
#     TSvol(dat[5]['psal'],dat[5]['ptmp'],axx[2,0],b1,b2)
#     #After
#     a1='2015-4-15'
#     a2='2015-5-15'
#     TSvol(FLA_daily['salinity'],FLA_daily['temperature'],axx[0,1],a1,a2)
#     TSvol(dat[8]['psal'],dat[8]['ptmp'],axx[1,1],a1,a2)
#     TSvol(dat[5]['psal'],dat[5]['ptmp'],axx[2,1],a1,a2)
#
#     axx[0,0].set_ylim(2.5,6.5)
#     axx[0,0].set_xlim(34.85,35.02)
#
#
# TSvolpans()

## I don't think the above is actually very successful in showing what I want, need something to compare cold, fresh bump with


XXXXXXXXXXXXX

# def plotPV(ii,axx,fff='yes'):
#     if fff=='yes':
#         f,axx=subplots(1,1,figsize=(12,3))
#     hh=axx.pcolor(dat[ii].date,dat[ii].depth_reg[::10],log10(dat[ii].PV[::10,:]),cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)
#     axx.plot(dat[ii].date,dat[ii].ML,'k.')
#     if fff=='yes':
#         colorbar(hh,label='log10(PPV) [(ms)$^{-1}$]')
#         ylabel('depth [m]')
#     if ii==8:
#         title('M1')
#     else:
#         title('CF'+str(ii))
#
#     if fff=='yes':
#         savefig(figdir+'MixedLayer/PV_CF'+str(ii)+'_fulldepth_hrly.png',bbox_inches='tight')
#         ylim(1300,0)
#         savefig(figdir+'MixedLayer/PV_CF'+str(ii)+'_top1300_hrly.png',bbox_inches='tight')
#         ylim(2000,0)
#         savefig(figdir+'MixedLayer/PV_CF'+str(ii)+'_top2000_hrly.png',bbox_inches='tight')
#     return hh
#
#
#
# def PV_presplot():
#     f,(ax1,ax2)=subplots(2,1,figsize=(9,5),sharex=True)
#     hh=plotPV(5,ax1,fff='no')
#     ax1.set_title('Slope current')
#     hh=plotPV(8,ax2,fff='no')
#     ax2.set_title('Just offshore')
#     cbar_ax = f.add_axes([0.95, 0.15, 0.03, 0.7])
#     cbar=colorbar(hh, cax=cbar_ax,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]')
#     ax1.set_ylim(1300,0)
#     ax2.set_ylim(1700,0)
#     f.text(-0.0025, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=12)
#     savefig(figdir+'MixedLayer/pres/PV_hov_CF5_M1_hrly.png',bbox_inches='tight',dpi=300)
#
# PV_presplot()
#
# plot(PV[8])
#
# univec['PV']=['PV',arange(-11.5,-10,0.025),cm.rainbow_r,arange(-11.5,-10,0.1),NaN]
#
# denvec=arange(27.5,27.8,0.01)
# middenvec=denvec[:-1]+diff(denvec)/2
# varvec=['PV','sal','tmp','dpth','thick']
# denmats={}
# for var in varvec:
#     denmats[var]=zeros((len(dat.distance),len(middenvec),len(dat.date)))
#     for moornum in range(8):
#         for dd in range(len(middenvec)):
#             if var=='dpth':
#                 denmats[var][moornum,dd,:]=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');
#             elif var=='thick':
#                 alldpths=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1])
#                 denmats[var][moornum,dd,:]=alldpths.max(dim='depth')-alldpths.min(dim='depth');
#             else:
#                 denmats[var][moornum,dd,:]=dat[univec[var][0]][moornum,:,:].where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');
#
#
# def makedendat():
#     dendat=xr.Dataset({'PV': (['distance','den','date'],  denmats['PV']),
#                        'depth': (['distance','den','date'],  denmats['dpth']),
#                        'thickness': (['distance','den','date'],  denmats['thick']),
#                        'sal': (['distance','den','date'],  denmats['sal']),
#                        'tmp': (['distance','den','date'],  denmats['tmp'])},
#                        coords={'den': middenvec,
#                                'distance': dat.distance.values,
#                                'date': dat.date.values})
#
#     return dendat
#
# dendat=makedendat()
#
# def plotDenThick(moornum,axx):
#     hh=axx.pcolor(dat.date.values,middenvec,log10(dendat['thickness'][moornum,:,:]),vmin=1,vmax=2.6,cmap=cm.rainbow)
#     axx.plot(dat.date.values,[27.67]*712,color='k',linewidth=2)
#     axx.plot(dat.date.values,[27.73]*712,color='k',linewidth=2)
#     axx.plot(dat.date.values,[27.79]*712,color='k',linewidth=2)
#     return hh
#
# def Den_thick_presplot():
#     f,(ax1,ax2)=subplots(2,1,figsize=(9,5),sharex=True,sharey=True)
#     cbar_ax = f.add_axes([0.95, 0.15, 0.03, 0.7])
#     hh=plotDenThick(4,ax1)
#     ax1.set_title('Slope current')
#     hh=plotDenThick(7,ax2)
#     ax2.set_title('Just offshore')
#     cbar=colorbar(hh,label='\n log$_{10}$( Thickness ) [m]',ticks=log10(array([10,30,100,300])),cax=cbar_ax)
#     cbar.ax.set_yticklabels(['10','30','100','300'])
#     ax1.set_ylim(max(middenvec),min(middenvec))
#     f.text(0.02, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=12)
#     ax2.text(datetime.datetime(2014,9,1),27.71,'upper IIW',fontsize=14)
#     ax2.text(datetime.datetime(2014,9,1),27.77,'deep IIW',fontsize=14)
#     savefig(figdir+'MixedLayer/pres/Den_thick_CF5_M1.png',bbox_inches='tight',dpi=300)
#
# Den_thick_presplot()
#
#
#
# middenvec[17]
# middenvec[18]
# middenvec[23]
# middenvec[29]
#
#
# def ThickComp_Tseries(d1,d2,ax1):
#         ddiff='W'
#         ax1.plot(dendat.date,nansum(dendat.thickness[4,d1:d2,:],axis=0),alpha=0.5)
#         ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[4,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),label='Slope current',color='C0',linewidth=2)
#         ax1.plot(dendat.date,nansum(dendat.thickness[-1,d1:d2,:],axis=0),label='',color='C1',alpha=0.5)
#         ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[-1,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),color='C1',linewidth=2,label='Offshore')
#
# dlist=['2014-9-15','2015-2-20','2015-4-10','2015-9-15']
#
# years=matplotlib.dates.YearLocator()
# months=matplotlib.dates.MonthLocator()
# threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
# monthFMT=matplotlib.dates.DateFormatter('%B')
# yearFMT=matplotlib.dates.DateFormatter('\n %Y')
#
# def ThickComp_pres():
#     f,(ax1,ax2)=subplots(2,1,figsize=(9,4),sharex=True,sharey=True)
#     ThickComp_Tseries(18,24,ax1)
#     ThickComp_Tseries(24,30,ax2)
#     ax1.legend(loc=1)
#     # ax1.set_xticks([])
#     ax1.set_title('upper IIW')
#     ax2.set_ylim(0,1200)
#     ax2.set_title('deep IIW')
#     f.text(0.02, 0.5, 'Layer thickness [m]', va='center', rotation='vertical',fontsize=12)
#     ax1.set_xlim(dat.date.values[0],dat.date.values[-1])
#     ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
#     ax2.xaxis.set_major_locator(years)
#     ax2.xaxis.set_minor_locator(threemonth)
#     ax2.xaxis.set_minor_formatter(monthFMT)
#     ax2.xaxis.set_major_formatter(yearFMT)
#     savefig(figdir+'MixedLayer/pres/CompThickness_CF5M1.pdf',bbox_inches='tight')
#     [ax1.axvline(dd,color='C'+str(ii+2),linewidth=3) for ii,dd in enumerate(dlist)]
#     [ax2.axvline(dd,color='C'+str(ii+2),linewidth=3) for ii,dd in enumerate(dlist)]
#
#     savefig(figdir+'MixedLayer/pres/CompThickness_CF5M1_wline.pdf',bbox_inches='tight')
#
# ThickComp_pres()
#
# def vertden_comp():
#             figure()
#             for ii,d1 in enumerate(dlist):
#                 figure(figsize=(2.5,2))
#                 print(d1)
#                 denchunk=dat['potential density'][3:,:,:].sel(date=d1).T
#                 contour(dat.distance[3:],dat.depth,denchunk,levels=[27.68,27.74,27.8],linewidths=3,colors='C'+str(ii+2))
#                 fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
#                 axvline(distvec[4],color='C0',linewidth=4)
#                 axvline(distvec[7],color='C1',linewidth=4)
#                 if ii!=0:
#                     gca().set_yticklabels('')
#                 else:
#                     ylabel('depth [m]',fontsize=12)
#                 ylim(1500,0)
#                 xlim(30,98)
#                 savefig(figdir+'MixedLayer/pres/VertDenSecComp_'+str(ii)+'.pdf',bbox_inches='tight')
#
# vertden_comp()
#
#
# # Potential density for properties at 750 - most relevant for ISW
# SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
# CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
# pdenmat2=zeros((shape(salmat)))
# for ii in range(len(salvec)):
#     for jj in range(len(tmpvec)):
#         pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3
#
#
#
# [FLA,FLB]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/FLMAB_xrays.pickle','rb'))
# FLA_daily=FLA.resample(date='D').mean()
# FLB_daily=FLB.resample(date='D').mean()
#
# def TS_BA(var,axx):
#     if var=='before':
#         d1='2014-9-15'
#         d2='2014-10-15'
#     elif var=='after':
#         d1='2015-4-15'
#         d2='2015-5-15'
#
#     axx.plot(FLA_daily['salinity'].sel(date=slice(d1,d2))[0,0],FLA_daily['temperature'].sel(date=slice(d1,d2))[0,0],'o',color='purple',label='Irminger gyre: OOI')
#     axx.plot(FLA_daily['salinity'].sel(date=slice(d1,d2)),FLA_daily['temperature'].sel(date=slice(d1,d2)),'o',color='purple',label='')
#     axx.plot(FLB_daily['salinity'].sel(date=slice(d1,d2)),FLB_daily['temperature'].sel(date=slice(d1,d2)),'o',color='purple',label='')
#     for mm in [7,4]:
#         if mm==7:
#             axx.plot(dat['salinity'][mm,0,:].sel(date=slice(d1,d2)),dat['temperature'][mm,0,:].sel(date=slice(d1,d2)),'o',color='C1',label='Just offshore')
#             axx.plot(dat['salinity'][mm,:,:].sel(date=slice(d1,d2)),dat['temperature'][mm,:,:].sel(date=slice(d1,d2)),'o',color='C1',alpha=0.05,label='')
#         else:
#             axx.plot(dat['salinity'][mm,0,:].sel(date=slice(d1,d2)),dat['temperature'][mm,0,:].sel(date=slice(d1,d2)),'o',color='C'+str(mm-4),label='Slope current')
#             axx.plot(dat['salinity'][mm,:,:].sel(date=slice(d1,d2)),dat['temperature'][mm,:,:].sel(date=slice(d1,d2)),'o',color='C'+str(mm-4),alpha=0.05,label='')
#     den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.68,27.74,27.8])
#
# def TSBA_pres():
#     f,(ax1,ax2)=subplots(1,2,figsize=(10,3),sharex=True,sharey=True)
#     TS_BA('before',ax1)
#     TS_BA('after',ax2)
#     ax2.text(34.98,4.7,'upper IIW',fontsize=12)
#     ax2.text(34.98,4.2,'deep IIW',fontsize=12)
#     ax1.set_ylim(2.5,6.5)
#     ax1.set_xlim(34.85,35.02)
#     ax1.set_xticks(arange(34.85,35.02,0.05))
#     ax1.set_title('October (before convection)')
#     ax2.set_title('May (after convection)')
#     ax2.legend(loc=(1.05,0.3))
#     ax1.set_ylabel('Temperature ($^\circ$C)')
#     f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=12)
#     savefig(figdir+'MixedLayer/pres/TS_BA.pdf',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/pres/TS_BA.png',bbox_inches='tight',dpi=300)
#
#
# TSBA_pres()
