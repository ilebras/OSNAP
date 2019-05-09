# versname='1809lpfilt_noextrap'
# daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_'+versname+'.pickle','rb'))
#
# #############################################################
# ############## Mixed layer calc  ######################
# #############################################################
# # get mixed layer depth, salinity, temp and pden
# MLD=ones((len(daily.date),len(distvec)))*NaN
# MLT=ones((len(daily.date),len(distvec)))*NaN
# MLS=ones((len(daily.date),len(distvec)))*NaN
# MLP=ones((len(daily.date),len(distvec)))*NaN
# PV=ones((len(distvec),len(daily.depth),len(daily.date)))*NaN
# dendiff={}
# for jj,dist in enumerate(distvec):
#     dd=int(dist)
#     print(dd)
#     dendiff[jj]=vstack((zeros(len(daily.date)),nancumsum(diff(daily['potential density'].sel(distance=dd),axis=0),axis=0)))
#     denthresh=0.0075
#     for ii in range(len(daily.date)):
#         dind=(dendiff[jj][:,int(ii)]<=denthresh) & ~isnan(daily['potential density'][jj,:,ii])
#         den_smooth=run_ave(daily['potential density'][jj,:,ii],25)
#         PV[jj,:,ii]=sw.f(CFlat[jj])/(den_smooth+1e3)*gradient(den_smooth)/gradient(daily.depth)
#         if sum(dind)>10:
#             MLD[ii,jj]=daily.depth[dind][-1].values
#             MLS[ii,jj]=std(daily['salinity'][jj,dind,ii])
#             MLT[ii,jj]=std(daily['temperature'][jj,dind,ii])
#             MLP[ii,jj]=std(daily['potential density'][jj,dind,ii])
#
#
# shape(MLD)
# MLDplus = xr.Dataset({'MLplus': (('date','distance'), MLD), 'date': daily.date.values,'distance':daily.distance.values})
# dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))
#
# newdat=xr.merge([dat, MLDplus])
#
#
# pickle.dump(newdat,open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_'+versname+'_wMLPV.pickle','wb'))
#
# #mask any mixed layers with standard deviations above
# MLD_masked=MLD.copy()
# masking=(MLS>0.005) | (MLT>0.05) | (MLP >0.05)
# MLD_masked[masking]=NaN
# sum(masking)
#
#
# def stdcheck():
#     f,axx=subplots(5,1,sharex=True,figsize=(12,14))
#     axx[0].plot(daily.date,MLS,'.')
#     axx[0].axhline(0.005,color='k')
#     axx[0].set_ylabel('std( ML salinity)')
#     # axx[0].set_ylim(0,0.01)
#     axx[1].plot(daily.date,MLT,'.')
#     axx[1].axhline(0.05,color='k')
#     axx[1].set_ylabel('std( ML temperature)')
#     # axx[1].set_ylim(0,0.1)
#     axx[2].plot(daily.date,MLP,'.')
#     axx[2].axhline(0.05,color='k')
#     # axx[2].set_ylim(0,0.1)
#     axx[2].set_ylabel('std( ML density)')
#     axx[3].plot(daily.date,MLD,'.')
#     axx[3].set_ylabel('ML depth')
#     axx[3].set_ylim(0,750)
#     axx[4].plot(daily.date,MLD_masked,'.')
#     axx[4].set_ylabel('ML depth (w/constraints)')
#     axx[4].set_ylim(0,750)
#
#
# stdcheck()
#
# daily
#
# shape(PV)
#
# newXR = xr.Dataset({'ML': (('date','distance'), MLD_masked), 'PV': (('distance','depth','date'), PV), 'moornum': (('distance'),range(1,9)), 'date': daily.date.values,'distance':daily.distance.values})
# daily=pickle.load(open('../pickles/xarray/CF_xarray_notid_'+versname+'.pickle','rb'))
# newdat=xr.merge([daily, newXR])
#
# pickle.dump(newdat,open('../pickles/xarray/CF_xarray_notid_'+versname+'_wMLPV.pickle','wb'))

from aux_funcs import *

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))
dat['across track velocity']=-1*dat['across track velocity']

def onecont(field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,addcont=0):
    if axx==0:
        fig,axx=subplots(1,1)
    ax1=axx.contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    if addcont==1:
        ax2=contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
        clabel(ax2,fmt='%1.1f')
    axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance [km]',fontsize=16)
    ylabel('depth [m]',fontsize=16)
    xlim([25,100])
    ylim([2200,0])
    title(tit,fontsize=22)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axx.set_yticks(range(500,2100,500))
    return ax1



#plot mean velocity and density levels for intro slides
def makesec():
        fig,axx=subplots(1,1,figsize=(5,4))
        field='uacross'
        axvel=onecont(dat[univec[field][0]].mean(dim='date').T,'',univec[field][1][:-3],cm.GnBu,univec[field][3],axx=axx,nomoorlines=0)
        cbaxes = fig.add_axes([0.95, 0.15, 0.04, 0.7])
        cbar=colorbar(axvel,ticks=univec[field][3],label='[m/s]',cax=cbaxes)
        axden=axx.contour(dat.distance,dat.depth,dat['potential density'].mean(dim='date').T,levels=[27.68,27.74,27.8],colors='k',linewidths=2)
        manual_locations = [(70,500), (70, 1000), (70,1500)]
        clabel(axden,fmt='%1.2f',manual=manual_locations)
        axx.text(distvec[4]-10,-100,'Boundary current max',color='black',fontsize=18,zorder=200)
        axx.axvline(distvec[4],color='white',linewidth=9)
        axx.axvline(distvec[4],color='',linewidth=5)
        axx.axvline(distvec[-1],color='white',linewidth=9)
        axx.text(distvec[-1]-10,-100,'Just offshore',color='black',fontsize=18,zorder=200)
        axx.axvline(distvec[-1],color='C1',linewidth=5)
        savefig(figdir+'MixedLayer/pres/Meansec_velden_v2.png',bbox_inches='tight')
        # savefig(figdir+'MixedLayer/pres/Meansec_velden.pdf',bbox_inches='tight')

makesec()

gridded=dat.copy()


# NOTE: SHOULD UPDATE THIS FOR NEW M1 MICROCAT DATA!!/SHOULD BE USED IN HOURLY ANALYSIS AT LEAST --done for hourly.

#go back to mooring resolution, not fully gridded...
dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))

hrly=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_hrly_wMLPV_1904.pickle','rb'))

def plotPV(moornum,axx,fff='yes',hrlyML='no'):
    if fff=='yes':
        f,axx=subplots(1,1,figsize=(12,3))
    hh=axx.pcolor(dat.date,dat.depth[::10],log10(dat.PV[moornum,::10,:]),cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)

    if hrlyML=='yes':
        if moornum==7:
            axx.plot(hrly[moornum+1].date,-hrly[moornum+1].ML,'k.',alpha=0.3)
        elif moornum==4:
            axx.plot(dat.date,dat.ML[:,moornum],'k.')
    else:
        axx.plot(dat.date,dat.MLplus[:,moornum],'.',color='grey')
        axx.plot(dat.date,dat.ML[:,moornum],'k.')
    if fff=='yes':
        colorbar(hh,label='log10(PPV) [(ms)$^{-1}$]')
        ylabel('depth [m]')
    if moornum==7:
        title('M1')
    else:
        title('CF'+str(moornum+1))

    if fff=='yes':
        savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_fulldepth.png',bbox_inches='tight')
        ylim(1300,0)
        savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_top1300.png',bbox_inches='tight')
        ylim(2000,0)
        savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_top2000.png',bbox_inches='tight')
    return hh

years=matplotlib.dates.YearLocator()
months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')



def PV_presplot(hM):
    f,(ax1,ax2)=subplots(2,1,figsize=(9,5),sharex=True)
    hh=plotPV(4,ax1,fff='no',hrlyML=hM)
    ax1.set_title('Slope current')
    hh=plotPV(7,ax2,fff='no',hrlyML=hM)
    ax2.set_title('Just offshore')
    cbar_ax = f.add_axes([0.95, 0.15, 0.03, 0.7])
    cbar=colorbar(hh, cax=cbar_ax,label='\n log$_{10}$( PPV ) [m$^{-1}$ s$^{-1}$]')
    ax1.set_ylim(1300,0)
    ax2.set_ylim(1700,0)
    f.text(-0.0025, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=12)
    ax2.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,20)])
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    if hM=='yes':
        savefig(figdir+'MixedLayer/pres/PV_hov_CF5_M1_hrly.png',bbox_inches='tight',dpi=300)
    else:
        savefig(figdir+'MixedLayer/pres/PV_hov_CF5_M1.png',bbox_inches='tight',dpi=300)

PV_presplot('yes')



# versname='1808lpfilt'
# [cc,egic,eg,ic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/transdic_'+versname+'.pickle','rb'))


univec['PV']=['PV',arange(-11.5,-10,0.025),cm.rainbow_r,arange(-11.5,-10,0.1),NaN]

denvec=arange(27.5,27.8,0.01)
middenvec=denvec[:-1]+diff(denvec)/2
varvec=['PV','sal','tmp','dpth','thick']
denmats={}
for var in varvec:
    denmats[var]=zeros((len(dat.distance),len(middenvec),len(dat.date)))
    for moornum in range(8):
        for dd in range(len(middenvec)):
            if var=='dpth':
                denmats[var][moornum,dd,:]=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');
            elif var=='thick':
                alldpths=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1])
                denmats[var][moornum,dd,:]=alldpths.max(dim='depth')-alldpths.min(dim='depth');
            else:
                denmats[var][moornum,dd,:]=dat[univec[var][0]][moornum,:,:].where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');


def makedendat():
    dendat=xr.Dataset({'PV': (['distance','den','date'],  denmats['PV']),
                       'depth': (['distance','den','date'],  denmats['dpth']),
                       'thickness': (['distance','den','date'],  denmats['thick']),
                       'sal': (['distance','den','date'],  denmats['sal']),
                       'tmp': (['distance','den','date'],  denmats['tmp'])},
                       coords={'den': middenvec,
                               'distance': dat.distance.values,
                               'date': dat.date.values})

    return dendat

dendat=makedendat()

dendat

def plotDenThick(moornum,axx):
    hh=axx.pcolor(dat.date.values,middenvec,log10(dendat['thickness'][moornum,:,:]),vmin=1,vmax=2.6,cmap=cm.rainbow)
    axx.plot(dat.date.values,[27.68]*712,color='k',linewidth=2)
    axx.plot(dat.date.values,[27.74]*712,color='k',linewidth=2)
    axx.plot(dat.date.values,[27.8]*712,color='k',linewidth=2)
    return hh

def Den_thick_presplot():
    f,(ax1,ax2)=subplots(2,1,figsize=(9,5),sharex=True,sharey=True)
    cbar_ax = f.add_axes([0.95, 0.15, 0.03, 0.7])
    hh=plotDenThick(4,ax1)
    ax1.set_title('Slope current')
    hh=plotDenThick(7,ax2)
    ax2.set_title('Just offshore')
    cbar=colorbar(hh,label='\n log$_{10}$( Thickness ) [m]',ticks=log10(array([10,30,100,300])),cax=cbar_ax)
    cbar.ax.set_yticklabels(['10','30','100','300'])
    ax1.set_ylim(max(middenvec),min(middenvec))
    f.text(0.02, 0.5, 'potential density, $\sigma_0$ [kg m$^{-3}$]', va='center', rotation='vertical',fontsize=12)
    ax2.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,20)])
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    ax2.text(datetime.datetime(2014,9,20),27.72,'upper IIW',fontsize=14)
    ax2.text(datetime.datetime(2014,9,20),27.78,'deep IIW',fontsize=14)
    savefig(figdir+'MixedLayer/pres/Den_thick_CF5_M1.png',bbox_inches='tight',dpi=300)

Den_thick_presplot()



middenvec[17]
middenvec[18]
middenvec[23]
middenvec[29]


def ThickComp_Tseries(d1,d2,ax1):
        ddiff='W'
        ax1.plot(dendat.date,nansum(dendat.thickness[4,d1:d2,:],axis=0),alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[4,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),label='Slope current',color='C0',linewidth=2)
        ax1.plot(dendat.date,nansum(dendat.thickness[-1,d1:d2,:],axis=0),label='',color='C1',alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[-1,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),color='C1',linewidth=2,label='Offshore')

dlist=['2014-9-15','2015-2-20','2015-4-10','2015-9-15']

def ThickComp_pres():
    f,(ax1,ax2)=subplots(2,1,figsize=(9,4),sharex=True,sharey=True)
    ThickComp_Tseries(18,24,ax1)
    ThickComp_Tseries(24,30,ax2)
    ax1.legend(loc=1)
    # ax1.set_xticks([])
    ax1.set_title('upper IIW')
    ax2.set_ylim(0,1200)
    ax2.set_title('deep IIW')
    f.text(0.02, 0.5, 'Layer thickness [m]', va='center', rotation='vertical',fontsize=12)
    # ax1.set_xlim(dat.date.values[0],dat.date.values[-1])
    ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)
    savefig(figdir+'MixedLayer/pres/CompThickness_CF5M1.pdf',bbox_inches='tight')
    [ax1.axvline(dd,color='C'+str(ii+2),linewidth=3) for ii,dd in enumerate(dlist)]
    [ax2.axvline(dd,color='C'+str(ii+2),linewidth=3) for ii,dd in enumerate(dlist)]

    savefig(figdir+'MixedLayer/pres/CompThickness_CF5M1_wline.pdf',bbox_inches='tight')

ThickComp_pres()



def vertden_comp():
            figure()
            for ii,d1 in enumerate(dlist):
                figure(figsize=(2.5,2))
                print(d1)
                denchunk=dat['potential density'][3:,:,:].sel(date=d1).T
                contour(dat.distance[3:],dat.depth,denchunk,levels=[27.68,27.74,27.8],linewidths=3,colors='C'+str(ii+2))
                fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
                axvline(distvec[4],color='C0',linewidth=4)
                axvline(distvec[7],color='C1',linewidth=4)
                if ii!=0:
                    gca().set_yticklabels('')
                else:
                    ylabel('depth [m]',fontsize=12)
                ylim(1500,0)
                xlim(30,98)
                savefig(figdir+'MixedLayer/pres/VertDenSecComp_'+str(ii)+'.pdf',bbox_inches='tight')

vertden_comp()

def plotden(ii):
    figure(figsize=(12,3))
    plot(dat.date,dat['potential density'][:,ii,:].T);
    title(str(ii*2)+'m')
### Look at the horizontal density gradient evolution in time, geostrophic velocity change vs. velocity pattern?
plotden(250)


versname='1810JHcal'
[cc,egic,eg,ic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/transdic_'+versname+'.pickle','rb'))


def comp_dengrad():
    f,ax1=subplots(1,1,figsize=(12,3))
    ax1.plot(dat.date,egic['trans'])
    ax2=ax1.twinx()
    ax2.plot(dat.date,(dat['potential density'][-1,150,:].T)-(dat['potential density'][3,150,:].T),color='C1',alpha=0.8)
    # ax2.plot(dat.date,(dat['potential density'][-1,500,:].T)-(dat['potential density'][4,500,:].T),color='C2',alpha=0.5)

comp_dengrad()

### Plot some rolling standard deviations of velocity data to see if there's an eddy pattern
## Looks fine but this should really be done with hourly data...

def rolling_std(var,binstd):
    rollin=zeros(len(dat.date))
    for ii in range(binstd):
        rollin[ii]=std(var[:ii+binstd])
    for ii in range(binstd,len(dat.date)-binstd):
        rollin[ii]=std(var[ii-binstd:ii+binstd])
    for ii in range(len(dat.date)-binstd,len(dat.date)):
            rollin[ii]=std(var[ii-binstd:])
    figure(figsize=(12,3))
    plot(dat.date,rollin)


rolling_std(egic['trans'],20)


rolling_std(dat['across track velocity'][4,150,:],25)


rolling_std(dat['along track velocity'][4,150,:],20)
botvel={}
botvel[4]=xr.open_dataset(datadir+'OSNAP2016recovery/AQD_Data_CF/OS_OSNAP-CF4_201408_CM_411m.nc')
botvel[5]=xr.open_dataset(datadir+'OSNAP2016recovery/AQD_Data_CF/OS_OSNAP-CF5_201408_CM_1266.nc')

def roll_std_bot(var,binstd):
    rollin=zeros(len(var.TIME))
    for ii in range(binstd):
        rollin[ii]=std(var[:ii+binstd])
    for ii in range(binstd,len(var.TIME)-binstd):
        rollin[ii]=std(var[ii-binstd:ii+binstd])
    for ii in range(len(var.TIME)-binstd,len(var.TIME)):
            rollin[ii]=std(var[ii-binstd:])
    figure(figsize=(12,3))
    plot(var.TIME,rollin)


## Right thing to do might be, use a smaller window and then smooth it
## OR do a band-pass or the like, get some more objective/quantitative measure of the strength of high-freq.
## OR do a wavelet...
roll_std_bot((botvel[4].UCUR**2+botvel[4].VCUR**2),100)
roll_std_bot((botvel[5].UCUR**2+botvel[5].VCUR**2),100)

roll_std_bot((botvel[4].UCUR**2+botvel[4].VCUR**2),500)
roll_std_bot((botvel[5].UCUR**2+botvel[5].VCUR**2),500)


roll_std_bot((botvel[4].UCUR**2+botvel[4].VCUR**2),1000)
roll_std_bot((botvel[5].UCUR**2+botvel[5].VCUR**2),1000)

XXXXXXXXXXXXXXXXX

## Now doing this part with hourly data

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
#     elif var=='oct2015':
#         d1='2015-9-15'
#         d2='2015-10-15'
#     elif var=='may2016':
#         d1='2016-4-15'
#         d2='2016-5-15'
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
#     ax2.text(34.98,4.1,'deep IIW',fontsize=12)
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
#
#
# def TSBA_later():
#     f,(ax1,ax2)=subplots(1,2,figsize=(10,3),sharex=True,sharey=True)
#     TS_BA('oct2015',ax1)
#     TS_BA('may2016',ax2)
#     ax1.text(34.98,4.7,'upper IIW',fontsize=12)
#     ax1.text(34.98,4.1,'deep IIW',fontsize=12)
#     ax1.set_ylim(2.5,6.5)
#     ax1.set_xlim(34.85,35.02)
#     ax1.set_xticks(arange(34.85,35.02,0.05))
#     ax1.set_title('October 2015')
#     ax2.set_title('May 2016')
#     ax1.legend(loc=(0.3,-0.4),ncol=3)
#     ax1.set_ylabel('Temperature ($^\circ$C)')
#     f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=12)
#
# TSBA_later()
