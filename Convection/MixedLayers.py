# versname='1809lpfilt_noextrap'
# daily=pickle.load(open('../pickles/xarray/CF_xarray_notid_'+versname+'.pickle','rb'))

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


dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))

# versname='1808lpfilt'
# [cc,egic,eg,ic]=pickle.load(open('../pickles/transdic_'+versname+'.pickle','rb'))
#
# def plotML(moornum):
#     figure(figsize=(12,3))
#     plot(dat.date,dat.ML[:,moornum],'.-')
#     if moornum==7:
#         title('M1')
#     else:
#         title('CF'+str(moornum+1))
#     ylim(0,550)
#     ylabel('mixed layer depth [m]')
#     savefig(figdir+'MixedLayer/ML_CF'+str(moornum+1)+'.png',bbox_inches='tight')
#
# for rr in range(3,8):
#     plotML(rr)
#
#
#
# seltimes=[(datetime.datetime(2015,3,1),datetime.datetime(2015,3,11)),
#           (datetime.datetime(2016,1,1),datetime.datetime(2016,1,11)),
#           (datetime.datetime(2016,1,11),datetime.datetime(2016,1,23)),
#           (datetime.datetime(2016,1,25),datetime.datetime(2016,2,5)),
#           (datetime.datetime(2016,3,1),datetime.datetime(2016,3,11)),]
#
# snaptimes=[datetime.datetime(2014,10,1),
#            datetime.datetime(2015,1,1),
#            datetime.datetime(2015,4,1),
#            datetime.datetime(2015,7,1),
#            datetime.datetime(2015,10,1),
#            datetime.datetime(2016,1,1),
#            datetime.datetime(2016,4,1),
#            datetime.datetime(2016,7,1),]
#
#
#
# def transplot():
#     egic['trans'].plot(figsize=(12,3))
#
#     ylabel('[Sv]')
#     for ss in seltimes:
#         axvline(ss[0],color='k')
#         axvline(ss[1],color='k')
#     axhline(-27,color='k')
#     xlim(datetime.datetime(2015,8,1),datetime.datetime(2016,8,1))
#
# # transplot()
#
# def plotPVtseries():
#     figure(figsize=(12,3))
#     for moornum in range(3,8):
#         plot(dat.date,log10(dat.PV[moornum,:,:].where(dat.depth<300).mean(dim='depth')),label='CF'+str(moornum+1))
#     xlim(datetime.datetime(2015,8,1),datetime.datetime(2016,8,1))
#     for ss in seltimes:
#             axvline(ss[0],color='k')
#             axvline(ss[1],color='k')
#     grid('on')
#     legend()
#
# # plotPVtseries()
#
# def stack_tseries():
#     f,axx=subplots(3,1,sharex=True,figsize=(14,12))
#     for moornum in range(4,8):
#         axx[0].plot(dat.date,log10(dat.PV[moornum,:,:].where(dat.depth<300).mean(dim='depth')),label='CF'+str(moornum+1))
#         # axx[0].plot(dat.date,dat.PV[moornum,:,:].where(dat.depth<300).mean(dim='depth'),label='CF'+str(moornum+1))
#     axx[0].grid('on')
#     axx[0].set_ylabel('log10(PV) above 300m')
#     for moornum in range(4,8):
#         axx[1].plot(dat.date,log10(dat.PV[moornum,:,:].where((dat.depth>500) & (dat.depth<1e3)).mean(dim='depth')),label='CF'+str(moornum+1))
#         # axx[1].plot(dat.date,dat.PV[moornum,:,:].where((dat.depth>500) & (dat.depth<1e3)).mean(dim='depth'),label='CF'+str(moornum+1))
#         axx[1].grid('on')
#     axx[1].legend()
#     axx[1].set_ylabel('log10(PV) between 500 and 1000m')
#     egic['trans'].plot(ax=axx[2])
#     axx[2].grid('on')
#     axx[2].set_ylabel('EGIC transport [Sv]')
#     # for ss in snaptimes:
#     #     for ii in range(3):
#     #         axx[ii].axvline(ss,color='grey')
#     savefig(figdir+'MixedLayer/PV_xport_tseries.png',bbox_inches='tight')
#
# # stack_tseries()
#
def plotPV(moornum):
    figure(figsize=(12,3))
    pcolor(dat.date,dat.depth[::10],log10(dat.PV[moornum,::10,:]),cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)
    plot(dat.date,dat.ML[:,moornum],'k.')
    colorbar(label='log10(PV)')
    if moornum==7:
        title('M1')
    else:
        title('CF'+str(moornum+1))
    gca().invert_yaxis()
    ylabel('depth [m]')
    axvline('2014-10-1',color='k',linewidth=3)
    axvline('2015-5-1',color='k',linewidth=3)
    savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_fulldepth_wline.png',bbox_inches='tight')
    ylim(1300,0)
    savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_top1300_wline.png',bbox_inches='tight')
    ylim(2000,0)
    savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_top2000_wline.png',bbox_inches='tight')
    # xlim(datetime.datetime(2016,1,1),datetime.datetime(2016,4,15))
    # for ss in seltimes:
    #     axvline(ss[0],color='k')
    #     axvline(ss[1],color='k')
    # savefig(figdir+'MixedLayer/PV_CF'+str(moornum+1)+'_winterzoom.png',bbox_inches='tight')

plotPV(4)
plotPV(7)
# for rr in range(3,8):
    # plotPV(rr)

#### Get hovmueller of PV within density classes:
##Build xarray dataset binned by density instead of depths
#
# def DenHov(d1,d2):
#     PVspec=dat.PV.where(dat['potential density']>d1).where(dat['potential density']<d2).mean(dim='depth')
#     figure(figsize=(12,3))
#     pcolor(dat.date,dat.distance,log10(PVspec),cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)
#     ylabel('distance [km]')
#     colorbar(label='log10(PV)')
#     ylim(45,90)
#
# DenHov(27.74,27.8)
#
# DenHov(27.68,27.74)
denvec

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

denmats['thick']

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

# Look at the evolution of properties in a density layer
def plotvarsonden(densel):
    f,axx=subplots(3,1,sharex=True,figsize=(25,8))
    for ii,var in enumerate(['sal','tmp','depth']):
        for rr in range(4,8):
            axx[ii].plot(dendat.date,dendat[var][rr,densel,:].T,label='CF'+str(rr+1));
        axx[ii].set_ylabel(var)
    axx[0].set_title(str(int(middenvec[densel]*1e3)/1e3)+' isopycnal')
    axx[1].legend(loc=(1.01,0.2))

plotvarsonden(-6)

plotvarsonden(-7)

plotvarsonden(-8)

# Look at the evolution of properties in a density layer
def plot_vardiff_onden(densel):
    f,axx=subplots(3,1,sharex=True,figsize=(15,8))
    for ii,var in enumerate(['sal','tmp','depth']):
        for rr in range(5,8):
            axx[ii].plot(dendat.date,(dendat[var][4,densel,:]-dendat[var][rr,densel,:]).T,label='CF5-CF'+str(rr+1));
        ddiff='W'
        axx[ii].plot(dendat.resample(date=ddiff).mean().date,(dendat[var][4,densel,:]-dendat[var][-2,densel,:]).resample(date=ddiff).mean().T,label='CF5-CF7, smoothed',color='k',linewidth='3');
        axx[ii].set_ylabel(var)
        axx[ii].axhline(0,color='grey')
    axx[0].set_title(str(int(middenvec[densel]*1e3)/1e3)+' isopycnal')
    axx[1].legend(loc=(1.01,0.2))
    # print('correlation between saldiff and dendiff CF5-M1')
    # print(corrcoef(dendat['sal'][4,densel,:]-dendat['sal'][-1,densel,:],dendat['depth'][4,densel,:]-dendat['depth'][-1,densel,:] )[1,0])

plot_vardiff_onden(-10)
plot_vardiff_onden(-7)

help(twinx)

def comp_diffs(densel):
    f,ax1=subplots(1,1,figsize=(15,3))
    ddiff='M'
    var='sal'
    ax1.plot(dendat.resample(date='3D').mean().date,(dendat[var][4,densel,:]-dendat[var][-2,densel,:]).resample(date='3D').mean().T,color='C2',alpha=0.5);
    ax1.plot(dendat.resample(date=ddiff).mean().date,(dendat[var][4,densel,:]-dendat[var][-2,densel,:]).resample(date=ddiff).mean().T,color='C2',linewidth='3');
    ax1.set_ylabel('salinity difference',color='C2')
    ax1.axhline(0,color='C2',alpha=0.5,linewidth=2)
    ax1.set_ylim(-0.02,0.03)
    ax2=twinx(ax=ax1)
    ax2.set_ylim(100,700)
    var='depth'
    ax2.plot(dendat.resample(date='3D').mean().date,(dendat[var][4,densel,:]-dendat[var][-2,densel,:]).resample(date='3D').mean().T,color='C1',alpha=0.5);
    ax2.plot(dendat.resample(date=ddiff).mean().date,(dendat[var][4,densel,:]-dendat[var][-2,densel,:]).resample(date=ddiff).mean().T,color='C1',linewidth='3');
    ax2.set_ylabel('depth difference',color='C1')
    ax1.set_title(str(int(middenvec[densel]*1e3)/1e3)+' isopycnal (CF5-CF7)')

comp_diffs(-11)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_695.png',bbox_inches='tight')

comp_diffs(-10)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_705.png',bbox_inches='tight')

comp_diffs(-9)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_715.png',bbox_inches='tight')

comp_diffs(-8)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_725.png',bbox_inches='tight')

comp_diffs(-7)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_735.png',bbox_inches='tight')

comp_diffs(-6)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_745.png',bbox_inches='tight')

comp_diffs(-5)
savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_month3d_755.png',bbox_inches='tight')

#
# comp_diffs(-7)
# savefig(figdir+'MixedLayer/HorizGrads/SalDpthDiff_CF5-CF7_3day.png',bbox_inches='tight')


winterchunks=['2014-8-1','2014-9-1','2014-11-15','2014-12-15','2015-1-15','2015-2-15','2015-3-15','2015-4-15','2015-7-15','2015-11-15','2015-12-15','2016-1-15','2016-2-15','2016-3-15','2016-4-15','2016-8-1']

def vertden_evol():
    chunks=winterchunks
    for ii in range(len(chunks)-1):
            figure()
            denchunk=dat['potential density'][4:,:,:].sel(date=slice(chunks[ii],chunks[ii+1])).T
            [contour(dat.distance[4:],dat.depth,denchunk.isel(date=tt),levels=[27.7,27.74,27.8],colors='grey') for tt in range(len(denchunk.date))]
            contour(dat.distance[4:],dat.depth,denchunk.mean(dim='date'),levels=[27.7,27.74,27.8],colors='r',linewidths=3)
            title(chunks[ii]+' to '+chunks[ii+1])
            ylim(1500,0)
            savefig(figdir+'MixedLayer/HorizGrads/VertDenSecFrom'+chunks[ii]+'_to_'+chunks[ii+1]+'.png',bbox_inches='tight')
vertden_evol()


middenvec[-11]
#Plot evolution of depth 27.7 and 27.74 isopycnals at CF5 and M1
def compdendepth():
    f,axx=subplots(2,1,sharex=True,sharey=True,figsize=(12,5))
    axx[1].plot(dendat.date,dendat.depth[4,-7,:],label='CF5')
    axx[1].plot(dendat.date,dendat.depth[-1,-7,:],label='M1')
    axx[1].set_title('27.74 isopycnal')
    axx[0].plot(dendat.date,dendat.depth[4,-11,:],label='CF5')
    axx[0].plot(dendat.date,dendat.depth[-1,-11,:],label='M1')
    axx[0].set_title('27.7 isopycnal')
    axx[0].legend(loc=(1.05,-0.2))
    ylim(1300,0)
    savefig(figdir+'MixedLayer/HorizGrads/CompIsopycDepth_CF5M1.png',bbox_inches='tight')

compdendepth()

int((middenvec[-5]+middenvec[-6])/2*100)/100

middenvec




def plotdenthickness():
    for ii in range(-5,-11,-1):
        figure(figsize=(15,3))
        plot(dendat.date,dendat.depth[4,ii,:]-dendat.depth[4,ii-1,:],label='CF5')
        plot(dendat.date,dendat.depth[-1,ii,:]-dendat.depth[-1,ii-1,:],label='M1')
        legend(loc=(1.05,0.2))
        title('isopyncal thickness of '+str(int((middenvec[ii]+middenvec[ii-1])/2*100)/100)+'$\pm$ 0.005')
        savefig(figdir+'MixedLayer/HorizGrads/CompIsopycThickness_CF5M1_'+str(ii)+'.png',bbox_inches='tight')

plotdenthickness()

def plotdenthickness_axtwin():
    for ii in range(-5,-11,-1):
        f,ax1=subplots(1,1,figsize=(15,3))
        ax1.plot(dendat.date,dendat.depth[4,ii,:]-dendat.depth[4,ii-1,:],label='CF5')
        ax1.set_ylabel('CF5',color='C0')
        ax2=twinx(ax=ax1)
        ax2.plot(dendat.date,dendat.depth[-1,ii,:]-dendat.depth[-1,ii-1,:],label='M1',color='C1')
        ax2.set_ylabel('M1',color='C1')
        ax2.set_ylim(0,500)
        # legend(loc=(1.05,0.2))
        title('isopyncal thickness of '+str(int((middenvec[ii]+middenvec[ii-1])/2*100)/100)+'$\pm$ 0.005')
        savefig(figdir+'MixedLayer/HorizGrads/CompIsopycThickness_CF5M1_axtwin_'+str(ii)+'.png',bbox_inches='tight')

plotdenthickness_axtwin()


dat['surfden']=dat.ML.copy()
dat['surftmp']=dat.ML.copy()
dat['surfsal']=dat.ML.copy()
dat['minPVden']=dat.ML.copy()
dat['minPVtmp']=dat.ML.copy()
dat['minPVsal']=dat.ML.copy()
dat['minPV']=dat.ML.copy()
dat['minPVdepth']=dat.ML.copy()
dmaxind=[500,500,500,500,500,750,750,750]
for mm in range(8):
    for dd in range(len(dat.date)):
        pden=dat['potential density'][mm,:,dd]
        if sum(~isnan(pden))>0:
            dat['surfden'][dd,mm]=pden[~isnan(pden)][0]
            dat['surfsal'][dd,mm]=dat['salinity'][mm,~isnan(pden),dd][0]
            dat['surftmp'][dd,mm]=dat['temperature'][mm,~isnan(pden),dd][0]
        if sum(~isnan(dat.PV[mm,:,dd]))>0:
            minPVind=where(nanmin(dat.PV[mm,:dmaxind[mm],dd])==dat.PV[mm,:dmaxind[mm],dd])[0][0]
            dat['minPVden'][dd,mm]=dat['potential density'][mm,minPVind,dd]
            dat['minPVsal'][dd,mm]=dat['salinity'][mm,minPVind,dd]
            dat['minPVtmp'][dd,mm]=dat['temperature'][mm,minPVind,dd]
            dat['minPVdepth'][dd,mm]=dat.depth[minPVind]
            dat['minPV'][dd,mm]=nanmin(dat.PV[mm,:dmaxind[mm],dd])


def MLvsminPV():
    for ii in range(4,8):
        figure(figsize=(12,3))
        plot(dat.date,dat['minPVden'][:,ii],label='minimum PV density')
        plot(dat.date,dat['surfden'][:,ii],label='surface density')
        axvline(datetime.datetime(2014,11,15))
        axvline(datetime.datetime(2015,4,15))
        axvline(datetime.datetime(2015,11,15))
        axvline(datetime.datetime(2016,4,15))
        title('CF'+str(ii+1))
        ylim(27.85,27.4)
        axhline(27.74,color='C2')
        legend(loc=(1.05,0.3))
        savefig(figdir+'MixedLayer/DenComp_MLPV_CF'+str(moornum+1)+'.png',bbox_inches='tight')

# MLvsminPV()

# figure(figsize=(12,3))
# plot(dat.date,dat['minPVden'][:,4])
# plot(dat.date,dat['surfden'][:,-1])
#
# ylim(27.8,27.4)

# dat['sal_whenML']=dat['salinity'].copy()
# dat['tmp_whenML']=dat['temperature'].copy()
#
# def NaNout():
#     MLind=isnan(dat.ML[:,-1])
#     dat['minPVden'][MLind,4]=NaN
#     LSind=dat['minPVden'][:,-1]<27.74
#     dat['minPVden'][LSind,-1]=NaN
#     dat['minPVsal'][LSind,-1]=NaN
#     dat['minPVtmp'][LSind,-1]=NaN
#     for ii in range(8):
#         dat['surfden'][MLind,ii]=NaN
#         dat['surfsal'][MLind,ii]=NaN
#         dat['surftmp'][MLind,ii]=NaN
#
#         dat['sal_whenML'][ii,:,MLind]=NaN
#         dat['tmp_whenML'][ii,:,MLind]=NaN

# NaNout()



def DenCompPlot():
    figure(figsize=(12,3))
    plot(dat.date,dat['minPVden'][:,4],label='Density of PV minimum at CF5')
    plot(dat.date,dat['surfden'][:,-1],label='Mixed layer density at M1')
    # plot(dat.date,dat['surfden'][:,-2],label='Mixed layer density at CF7')
    # plot(dat.date,dat['surfden'][:,-3],label='Mixed layer density at CF6')
    # plot(dat.date,dat['surfden'][:,4],label='Mixed layer density at CF5')
    title('Comparing CF5 and M1 densities at all times with analysis chunks')
    axvline(datetime.datetime(2014,11,15))
    axvline(datetime.datetime(2015,4,15))
    axvline(datetime.datetime(2015,11,15))
    axvline(datetime.datetime(2016,4,15))
    axhline(27.74,color='grey')
    legend(loc=(1.01,0.2))
    ylim(27.8,27.4)
    ylabel('potential density')
    savefig(figdir+'MixedLayer/DenComp_MLPV_alltimes.png',bbox_inches='tight')


DenCompPlot()

# Potential density for properties at 750 - most relevant for ISW
SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat2=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3

[FLA,FLB]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/FLMAB_xrays.pickle','rb'))


FLA_daily=FLA.resample(date='D').mean()
FLB_daily=FLB.resample(date='D').mean()



chunks=['2014-8-1','2014-11-15','2015-4-15','2015-11-15','2016-4-15','2016-8-1']



def OOI_TS():
    for ii in range(len(chunks)-1):
            dlim=7e3
            figure(figsize=(10,8))
            plot(FLA_daily['salinity'].where(FLA.depth<=dlim),FLA_daily['temperature'].where(FLA.depth<=dlim),'o',color='grey',label='')
            plot(FLB_daily['salinity'].where(FLB.depth<=dlim),FLB_daily['temperature'].where(FLB.depth<=dlim),'o',color='grey',label='')
            plot(FLA_daily['salinity'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=dlim),FLA_daily['temperature'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=dlim),'ro',label='')
            plot(FLB_daily['salinity'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLB.depth<=dlim),FLB_daily['temperature'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLB.depth<=dlim),'o',color='purple',label='')
            plot(34.885,3.6,'yo',markersize=20,label='Irminger properties (Piron et al. 2017)')
            plot(34.84,3.3,'co',markersize=20,label='Labrador properties (Piron et al. 2017)')
            den=contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.7,27.74,27.8])
            title('OOI flanking moorings from '+chunks[ii]+' to '+chunks[ii+1])
            legend(loc=(1.05,0.5))
            xlim(34.6,35.1)
            ylim(0,10)
            savefig(figdir+'MixedLayer/TSchunks/OOIallTSfrom'+chunks[ii]+'_to_'+chunks[ii+1]+'.png',bbox_inches='tight')

OOI_TS()

#override here to examine in a bit more detail
chunks=winterchunks

def TScomp_OOI():
    for ii in range(len(chunks)-1):
        for mm in range(4,8):
            figure(figsize=(10,8))
            plot(dat['salinity'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),dat['temperature'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='grey',alpha=0.2,label='')
            plot(FLA_daily['salinity'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=50),FLA_daily['temperature'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=50),'ro',label='')
            plot(FLB_daily['salinity'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=50),FLB_daily['temperature'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLB.depth<=50),'o',color='purple',label='')
            plot(FLA_daily['salinity'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=50)[0,0],FLA_daily['temperature'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=50)[0,0],'ro',label='OOI FLA above 50m')
            plot(FLB_daily['salinity'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLA.depth<=50)[0,0],FLB_daily['temperature'].sel(date=slice(chunks[ii],chunks[ii+1])).where(FLB.depth<=50)[0,0],'o',color='purple',label='OOI FLB above 50m')
            plot(dat['minPVsal'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),dat['minPVtmp'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='orange',label='CF'+str(mm+1)+' PV minimum properties')
            den=contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.7,27.74,27.8])
            title('CF'+str(mm+1)+' from '+chunks[ii]+' to '+chunks[ii+1])
            legend(loc=(1.05,0.5))
            xlim(34.6,35.1)
            ylim(0,10)
            savefig(figdir+'MixedLayer/TSchunks/OOIsurfcomp_TSfrom'+chunks[ii]+'_to_'+chunks[ii+1]+'CF'+str(mm+1)+'.png',bbox_inches='tight')

TScomp_OOI()

chunks




def TScomplot():
    for ii in range(len(chunks)-1):
        for mm in range(4,8):
            figure(figsize=(10,8))
            plot(dat['salinity'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),dat['temperature'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='grey',alpha=0.2,label='')
            plot(dat['surfsal'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),dat['surftmp'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),'o',label='Surface properties')
            plot(dat['minPVsal'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),dat['minPVtmp'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),'o',label='PV minimum properties')
            den=contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.7,27.74,27.8])
            title('CF'+str(mm+1)+' from '+chunks[ii]+' to '+chunks[ii+1])
            legend(loc=(1.05,0.5))
            xlim(34.6,35.1)
            ylim(0,10)
            savefig(figdir+'MixedLayer/TSchunks/CF'+str(mm+1)+'_TSfrom'+chunks[ii]+'_to_'+chunks[ii+1]+'.png',bbox_inches='tight')


TScomplot()


def TScomplot():
    for ii in range(len(chunks)-1):
            figure(figsize=(10,8))
            mm=-1
            plot(dat['salinity'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),dat['temperature'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='grey',alpha=0.2,label='')
            mm=4
            plot(dat['salinity'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),dat['temperature'][mm,:,:].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='grey',alpha=0.2,label='')
            plot(dat['minPVsal'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),dat['minPVtmp'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='C1',label='PV minimum properties at CF5')
            mm=-1
            plot(dat['surfsal'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),dat['surftmp'][:,mm].sel(date=slice(chunks[ii],chunks[ii+1])),'o',color='C2',label='Surface properties at M1')

            den=contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.7,27.74,27.8])
            title(chunks[ii]+' to '+chunks[ii+1])
            legend(loc=(1.05,0.5))
            xlim(34.6,35.1)
            ylim(0,10)
            savefig(figdir+'MixedLayer/TSchunks/CompM1CF5_TSfrom'+chunks[ii]+'_to_'+chunks[ii+1]+'.png',bbox_inches='tight')

TScomplot()

def DenHovMoor(moornum,var):
    figure(figsize=(12,3))
    if var=='PV':
        pcolor(dat.date,middenvec,log10(dendat.PV[moornum,:,:]),vmin=min(univec[var][1]),vmax=max(univec[var][1]),cmap=univec[var][2])
        if moornum==4:
            plot(dat.date,dat.surfden[:,-1],color='k')
        colorbar(label='log10(PV)')
    else:
        pcolor(dat.date.values,middenvec,dendat[var][moornum,:,:],cmap=univec[var][2],vmin=min(univec[var][1]),vmax=max(univec[var][1]))
        colorbar(label=univec[var][0])
    if var=='sal':
        contour(dat.date.values,middenvec,dendat[var][moornum,:,:],levels=[34.9],colors='k')
    ylabel('potential density')
    title('CF'+str(moornum+1))
    ylim(max(middenvec),min(middenvec))
    axhline(27.74,color='grey')
    savefig(figdir+'MixedLayer/Hovden_'+var+'CF'+str(moornum+1)+'.png',bbox_inches='tight')

DenHovMoor(4,'PV')

for moornum in range(4,8):
    figure(figsize=(12,3))
    var='thickness'
    pcolor(dat.date.values,middenvec,log10(dendat[var][moornum,:,:]),vmin=1,vmax=2.6,cmap=cm.rainbow)
    axhline(27.67,color='grey')
    axhline(27.73,color='grey')
    axhline(27.79,color='grey')
    cbar=colorbar(ticks=log10(array([10,30,100,300])))
    cbar.ax.set_yticklabels(['10','30','100','300'])
    cbar.set_label('Thickness [m]')
    ylim(max(middenvec),min(middenvec))
    savefig(figdir+'MixedLayer/HovdenThick_CF'+str(moornum+1)+'.png',bbox_inches='tight')
#Plot thickness within these two ranges: 26.67-27.73-27.79
middenvec[17]
middenvec[23]
middenvec[29]

dlist=['2014-9-15','2015-2-20','2015-4-20','2015-9-15']

def plotthickness_range(d1,d2):
        ddiff='W'
        f,ax1=subplots(1,1,figsize=(15,3),sharex=True)
        ax1.plot(dendat.date,nansum(dendat.thickness[4,d1:d2,:],axis=0),alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[4,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),label='CF5',color='C0',linewidth=2)
        ax1.set_ylabel('Layer thickness [m]')

        ax1.plot(dendat.date,nansum(dendat.thickness[-1,d1:d2,:],axis=0),label='',color='C1',alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[-1,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),color='C1',linewidth=2,label='M1')
        ax1.set_ylim(0,1200)
        # [ax1.axvline(dd,color='C'+str(ii+2),linewidth=3) for ii,dd in enumerate(dlist)]
        ax1.legend()
        savefig(figdir+'MixedLayer/HorizGrads/CompThicknessRange_CF5M1_'+str(d1)+'-'+str(d2)+'.png',bbox_inches='tight')

plotthickness_range(17,23)

plotthickness_range(23,29)

distvec[4]

def vertden_comp():
            figure()
            for ii,d1 in enumerate(dlist):
                figure(figsize=(4,3))
                print(d1)
                denchunk=dat['potential density'][3:,:,:].sel(date=d1).T
                contour(dat.distance[3:],dat.depth,denchunk,levels=[27.67,27.73,27.79],linewidths=3,colors='C'+str(ii+2))
                fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
                axvline(distvec[4],color='C0',linewidth=4)
                axvline(distvec[7],color='C1',linewidth=4)
                if ii!=0:
                    gca().set_yticklabels('')
                ylim(1500,0)
                xlim(30,98)
                savefig(figdir+'MixedLayer/pres/VertDenSecComp_'+str(ii)+'.png',bbox_inches='tight')
distvec[7]

vertden_comp()
vertden_check((dlist[1])
vertden_check(dlist[2])

vertden_check(dlist[3])
def plotthickness_range_WSAL(d1,d2):
        ddiff='W'
        f,(ax1,ax2)=subplots(2,1,figsize=(15,5),sharex=True)
        ax1.plot(dendat.date,nansum(dendat.thickness[4,d1:d2,:],axis=0),alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[4,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),label='CF5',color='C0',linewidth=2)
        ax1.set_ylabel('Layer thickness [m]')

        ax1.plot(dendat.date,nansum(dendat.thickness[-1,d1:d2,:],axis=0),label='',color='C1',alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[-1,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),color='C1',linewidth=2,label='M1')
        ax1.set_ylim(0,1200)
        ax1.legend()
        # ax2=twinx(ax=ax1)
        ax2.plot(dendat.date,nanmean(dendat.sal[4,d1:d2,:]-dendat.sal[-1,d1:d2,:],axis=0),color='C2',alpha=0.5)
        ax2.plot(dendat.resample(date=ddiff).mean(dim='date').date,nanmean((dendat.sal[4,d1:d2,:]-dendat.sal[-1,d1:d2,:]).resample(date=ddiff).mean(dim='date'),axis=0),color='C2',linewidth=2)
        ax2.axhline(0,color='grey')
        ax2.set_title('Mean salinity anomaly in isopycnal range')
        # ax2.set_ylim(-0.2,0.5)
        savefig(figdir+'MixedLayer/HorizGrads/CompThicknessRange_PlusSal_CF5M1_'+str(d1)+'-'+str(d2)+'.png',bbox_inches='tight')

plotthickness_range_WSAL(17,23)

plotthickness_range_WSAL(23,29)

plotthickness_range(23,29)

def plotthickness_range(d1,d2):
        ddiff='W'
        f,ax1=subplots(1,1,figsize=(15,3))
        ax1.plot(dendat.date,nansum(dendat.thickness[4,d1:d2,:],axis=0),alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[4,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),label='CF5',color='C0',linewidth=2)
        ax1.set_ylabel('CF5',color='C0')
        ax1.plot(dendat.date,nansum(dendat.thickness[-1,d1:d2,:],axis=0),label='M1',color='C1',alpha=0.5)
        ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[-1,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),color='C1',linewidth=2)
        ax1.set_ylabel('M1',color='C1')
        ax1.set_ylim(0,1200)
        savefig(figdir+'MixedLayer/HorizGrads/CompIsopycThickness_CF5M1_'+str(d1)+'-'+str(d2)+'.png',bbox_inches='tight')


dendat
d1='2014-9-15'
d2='2014-10-15'
minsalind=argmin(dat['salinity'][4,:,:].sel(date=slice(d1,d2)),axis=0)

dat['salinity'][4,:,:].sel(date=slice(d1,d2))[minsalind,:]

dendat

dendat['sal'][4,:,:].sel(date=slice(d1,d2)),

def TS_BA(var):
    if var=='before':
        d1='2014-9-15'
        d2='2014-10-15'
    elif var=='after':
        d1='2015-4-15'
        d2='2015-5-15'
    figure()
    plot(FLA_daily['salinity'].sel(date=slice(d1,d2))[0,0],FLA_daily['temperature'].sel(date=slice(d1,d2))[0,0],'o',color='purple',label='OOI')
    plot(FLA_daily['salinity'].sel(date=slice(d1,d2)),FLA_daily['temperature'].sel(date=slice(d1,d2)),'o',color='purple',label='')
    plot(FLB_daily['salinity'].sel(date=slice(d1,d2)),FLB_daily['temperature'].sel(date=slice(d1,d2)),'o',color='purple',label='')
    for mm in range(4,8)[::-1]:
        if mm==7:
            plot(dat['salinity'][mm,0,:].sel(date=slice(d1,d2)),dat['temperature'][mm,0,:].sel(date=slice(d1,d2)),'o',color='C'+str(mm-4),label='M1')
        else:
            plot(dat['salinity'][mm,0,:].sel(date=slice(d1,d2)),dat['temperature'][mm,0,:].sel(date=slice(d1,d2)),'o',color='C'+str(mm-4),label='CF'+str(mm+1))
        plot(dat['salinity'][mm,:,:].sel(date=slice(d1,d2)),dat['temperature'][mm,:,:].sel(date=slice(d1,d2)),'o',color='C'+str(mm-4),alpha=0.5,label='')
        # if var=='before':
        #     d1='2014-8-15'
        #     d2='2014-9-15'
        #     maxsalind=argmax(dendat['sal'][4,:-1,:].sel(date=slice(d1,d2)),axis=1)
        #     salref=dendat['sal'][4,:-1,:].sel(date=slice(d1,d2))[:,maxsalind]
        #     plot(dendat['sal'][4,:-1,:].sel(date=slice(d1,d2))[:,maxsalind],dendat['tmp'][4,:-1,:].sel(date=slice(d1,d2))[:,maxsalind],'o',color='k',label='')
    den=contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.67,27.73,27.79])
    legend()
    xlim(34.85,35.0)
    ylim(2.5,6.5)
    savefig(figdir+'MixedLayer/TSchunks/TS'+var+'.png',bbox_inches='tight')
    return salref

salref=TS_BA('before')



TS_BA('after')


salanom=dendat['sal'][:,:-1,:]-salref

def plotsalanom_range(d1,d2):
        ddiff='W'
        f,ax1=subplots(1,1,figsize=(15,3))
        ax1.plot(dendat.date,nanmean(salanom[4,d1:d2,:],axis=0),alpha=0.5)
        axvline('2014-9-15')
        axhline(0)
        # ax1.plot(dendat.resample(date=ddiff).mean().date,nansum(dendat.thickness[4,d1:d2,:].resample(date=ddiff).mean(dim='date'),axis=0),label='CF5',color='C0',linewidth=2)
        # ax1.set_ylabel('CF5',color='C0')
        # savefig(figdir+'MixedLayer/HorizGrads/SalAnomCF5_'+str(d1)+'-'+str(d2)+'.png',bbox_inches='tight')


plotsalanom_range(17,23)

plotsalanom_range(23,29)


for rr in range(5,8):
    DenHovMoor(rr,'PV')


for rr in range(4,8):
    DenHovMoor(rr,'tmp')


for rr in range(3,8):
    DenHovMoor(rr,'sal')


skipno=20

def PVcont(axi,timech):
    ax1=axi.contourf(dat.distance[3:],dat.depth[::skipno],log10(dat.PV)[3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,50,vmin=-11.5,vmax=-10,cmap=cm.rainbow_r,extend='both')
    colorbar(ax1,ax=axi,label='log10(PV)')
    ax2=axi.contour(dat.distance[3:],dat.depth[::skipno],log10(dat.PV)[3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,levels=[-11.4,-11],colors='k')
    gca().invert_yaxis()
    fstr='%1.1f'
    clabels=clabel(ax2, fmt=fstr, fontsize=10)
    # for dd in distvec[3:]:
    #     axi.axvline(dd,color='k')


def contelse(var,axx,timech):
    ax1=axx.contourf(dat.distance[3:],dat.depth[::skipno],dat[univec[var][0]][3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,univec[var][1],cmap=univec[var][2],extend='both')
    colorbar(ax1,ax=axx,label=univec[var][0])
    ax2=axx.contour(dat.distance[3:],dat.depth[::skipno],dat[univec[var][0]][3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,levels=univec[var][3],colors='k')
    tit=var
    if ('den' in tit) | ('sal' in tit):
        fstr='%1.2f'
    elif 'tmp' in tit:
        fstr='%1.0f'
    elif  ('uac' in tit):
        fstr='%1.1f'
    clabels=clabel(ax2, fmt=fstr, fontsize=10)

def plotsecs(timech):
    #plot mean sections of PV, salinity,temperature,density and velocity for a select period of time
    f,axx=subplots(5,1,figsize=(8,24),sharex=True,sharey=True)
    PVcont(axx[0],timech)
    contelse('pden',axx[1],timech)
    contelse('uacross',axx[2],timech)
    contelse('sal',axx[3],timech)
    contelse('tmp',axx[4],timech)
    axx[2].set_ylabel('depth [m]')
    axx[4].set_xlabel('distance [km]')
    axx[0].set_title(str(timech[0].year)+'-'+str(timech[0].month)+'-'+str(timech[0].day)+' to '+
             str(timech[1].year)+'-'+str(timech[1].month)+'-'+str(timech[1].day),fontsize=16)
    savefig(figdir+'MixedLayer/Sections_bigtrans_'+str(timech[0].year)+'-'+str(timech[0].month)+'-'+str(timech[0].day)+' to '+
             str(timech[1].year)+'-'+str(timech[1].month)+'-'+str(timech[1].day)+'.png',bbox_inches='tight')

plotsecs(seltimes[0])

def plotPVdenvel(timech):
    contourf(dat.distance[3:],dat.depth[::skipno],log10(dat.PV)[3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,50,vmin=-11.5,vmax=-10,cmap=cm.rainbow_r,extend='both')
    var='uacross'
    contour(dat.distance[3:],dat.depth[::skipno],dat[univec[var][0]][3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,levels=univec[var][3][::2],colors='k')
    var='pden'
    contour(dat.distance[3:],dat.depth[::skipno],dat[univec[var][0]][3:,::skipno,:].sel(date=slice(timech[0],timech[1])).mean(dim='date').T,levels=univec[var][3],colors='purple',linewidths=4)
    gca().invert_yaxis()
    xlabel('distance [km]')
    ylabel('depth [m]')
    savefig(figdir+'MixedLayer/PV_wdenvel_'+str(timech[0].year)+'-'+str(timech[0].month)+'-'+str(timech[0].day)+' to '+
             str(timech[1].year)+'-'+str(timech[1].month)+'-'+str(timech[1].day)+'.png',bbox_inches='tight')

plotPVdenvel(seltimes[1])


def PVcont_daily(axi,timech):
    ax1=axi.contourf(dat.distance[3:],dat.depth[::skipno],log10(dat.PV)[3:,::skipno,:].sel(date=dat.date[timech]).T,50,vmin=-11.5,vmax=-10,cmap=cm.rainbow_r,extend='both')
    colorbar(ax1,ax=axi,label='log10(PV)')
    ax2=axi.contour(dat.distance[3:],dat.depth[::skipno],log10(dat.PV)[3:,::skipno,:].sel(date=dat.date[timech]).T,levels=[-11.4,-11],colors='k')
    gca().invert_yaxis()
    fstr='%1.1f'
    clabels=clabel(ax2, fmt=fstr, fontsize=10)
    # for dd in distvec[3:]:
    #     axi.axvline(dd,color='k')


def contelse_daily(var,axx,timech):
    ax1=axx.contourf(dat.distance[3:],dat.depth[::skipno],dat[univec[var][0]][3:,::skipno,:].sel(date=dat.date[timech]).T,univec[var][1],cmap=univec[var][2],extend='both')
    colorbar(ax1,ax=axx,label=univec[var][0])
    ax2=axx.contour(dat.distance[3:],dat.depth[::skipno],dat[univec[var][0]][3:,::skipno,:].sel(date=dat.date[timech]).T,levels=univec[var][3],colors='k')
    tit=var
    if ('den' in tit) | ('sal' in tit):
        fstr='%1.2f'
    elif 'tmp' in tit:
        fstr='%1.0f'
    elif  ('uac' in tit):
        fstr='%1.1f'
    clabels=clabel(ax2, fmt=fstr, fontsize=10)

univec['uacross']=['across track velocity',arange(-0.6,0.6,0.025),cm.RdBu_r,arange(-0.8,0.8,0.1),'[m/s]']

def plot_dailysecs():
    #plot daily mean sections of PV, salinity,temperature,density and velocity
    for timech in range(100):#len(dat.date)):
        f,axx=subplots(5,1,figsize=(8,24),sharex=True,sharey=True)
        PVcont_daily(axx[0],timech)
        contelse_daily('pden',axx[1],timech)
        contelse_daily('uacross',axx[2],timech)
        contelse_daily('sal',axx[3],timech)
        contelse_daily('tmp',axx[4],timech)
        axx[2].set_ylabel('depth [m]')
        axx[4].set_xlabel('distance [km]')
        axx[0].set_title(str(dat.date[timech].values)[:10],fontsize=16)
        savefig(figdir+'MixedLayer/daily/Sections_'+str(timech).zfill(3)+'.png',bbox_inches='tight')

plot_dailysecs()

def plotTSevolve():
    for ii in range(len(dat.date)):
        figure()
        mlim=4
        p1=0#300m
        p2=750#1500m
        for mm in range(4,8):
            if mm==7:
                labit='M1'
            else:
                labit='CF'+str(mm+1)
            plot(dat.salinity[mm,p1:p2,ii].values.flatten(),dat.temperature[mm,p1:p2,ii].values.flatten(),'.',label=labit);
        den=contour(salvec[76:81],tmpvec,pdenmat[:,76:81],colors='grey',levels=[27.5,27.68,27.74,27.8])
        clabel(den,fmt='%1.2f')
        legend()
        xlim(34.85,35.05)
        xticks(arange(34.85,35.06,0.05))
        ylim(2,6)
        xlabel('salinity')
        ylabel('potential temperature ($^\circ$C)')
        yticks(range(2,7,1))
        title(str(dat.date[ii].values)[:10],fontsize=16)
        savefig(figdir+'MixedLayer/TSdaily/TS_'+str(ii).zfill(3)+'.png',bbox_inches='tight')

plotTSevolve()

def plotTSPV():
    for ii in range(len(dat.date)):
        figure()
        mlim=4
        p1=0#300m
        p2=750#1500m
        for mm in range(4,8):
            if mm==7:
                labit='M1'
            else:
                labit='CF'+str(mm+1)
            scatter(dat.salinity[mm,p1:p2,ii].values.flatten(),dat.temperature[mm,p1:p2,ii].values.flatten(),
                    c=log10(dat.PV[mm,p1:p2,ii].values.flatten()),vmin=-11.5,vmax=-10,cmap=cm.rainbow_r);
        colorbar(label='log10(PV)')
        xlabel('salinity')
        ylabel('potential temperature ($^\circ$C)')
        den=contour(salvec[76:81],tmpvec,pdenmat[:,76:81],colors='grey',levels=[27.5,27.68,27.74,27.8])
        clabel(den,fmt='%1.2f')
        xlim(34.85,35.05)
        xticks(arange(34.85,35.06,0.05))
        ylim(2,6)
        yticks(range(2,7,1))
        title(str(dat.date[ii].values)[:10],fontsize=16)
        savefig(figdir+'MixedLayer/TSdaily/TS_PV_'+str(ii).zfill(3)+'.png',bbox_inches='tight')

plotTSPV()
help(hexbin)

def plotTSPVhex_30day():
    for ii in range(25):
        figure()
        mlim=4
        p1=0#300m
        p2=750#1500m
        for mm in range(4,8):
            if mm==7:
                labit='M1'
            else:
                labit='CF'+str(mm+1)
            hexbin(dat.salinity[mm,p1:p2,30*ii:30*(ii+1)].values.flatten(),dat.temperature[mm,p1:p2,30*ii:30*(ii+1)].values.flatten(),
                    C=log10(dat.PV[mm,p1:p2,30*ii:30*(ii+1)].values.flatten()),vmin=-11.5,vmax=-10,cmap=cm.rainbow_r);
        colorbar(label='log10(PV)')
        xlabel('salinity')
        ylabel('potential temperature ($^\circ$C)')
        den=contour(salvec[76:81],tmpvec,pdenmat[:,76:81],colors='grey',levels=[27.5,27.68,27.74,27.8])
        clabel(den,fmt='%1.2f')
        xlim(34.85,35.05)
        xticks(arange(34.85,35.06,0.05))
        ylim(2,6)
        yticks(range(2,7,1))
        title(str(dat.date[30*ii].values)[:7],fontsize=16)
        savefig(figdir+'MixedLayer/monthly/TS_PV30day_'+str(dat.date[30*ii].values)[:7]+'.png',bbox_inches='tight')

plotTSPVhex_30day()

def plotTSPVhex_30day():
    for ii in range(25):
        figure()
        mlim=4
        p1=0#300m
        p2=750#1500m
        for mm in range(4,8):
            if mm==7:
                labit='M1'
            else:
                labit='CF'+str(mm+1)
            hexbin(dat.salinity[mm,p1:p2,30*ii:30*(ii+1)].values.flatten(),dat.temperature[mm,p1:p2,30*ii:30*(ii+1)].values.flatten(),
                    C=log10(dat.PV[mm,p1:p2,30*ii:30*(ii+1)].values.flatten()),vmin=-11.5,vmax=-10,cmap=cm.rainbow_r);
        colorbar(label='log10(PV)')
        xlabel('salinity')
        ylabel('potential temperature ($^\circ$C)')
        den=contour(salvec[76:81],tmpvec,pdenmat[:,76:81],colors='grey',levels=[27.5,27.68,27.74,27.8])
        clabel(den,fmt='%1.2f')
        xlim(34.85,35.05)
        xticks(arange(34.85,35.06,0.05))
        ylim(2,6)
        yticks(range(2,7,1))
        title(str(dat.date[30*ii].values)[:7],fontsize=16)
        savefig(figdir+'MixedLayer/monthly/TS_PV30day_'+str(dat.date[30*ii].values)[:7]+'.png',bbox_inches='tight')

plotTSPVhex_30day()



def plotTSPVhex_10day():
    for ii in range(75):
        figure()
        mlim=4
        p1=0#300m
        p2=750#1500m
        for mm in range(4,8):
            if mm==7:
                labit='M1'
            else:
                labit='CF'+str(mm+1)
            hexbin(dat.salinity[mm,p1:p2,10*ii:10*(ii+1)].values.flatten(),dat.temperature[mm,p1:p2,10*ii:10*(ii+1)].values.flatten(),
                    C=log10(dat.PV[mm,p1:p2,10*ii:10*(ii+1)].values.flatten()),vmin=-11.5,vmax=-10,cmap=cm.rainbow_r);
        colorbar(label='log10(PV)')
        xlabel('salinity')
        ylabel('potential temperature ($^\circ$C)')
        den=contour(salvec[76:81],tmpvec,pdenmat[:,76:81],colors='grey',levels=[27.5,27.68,27.74,27.8])
        clabel(den,fmt='%1.2f')
        xlim(34.85,35.05)
        xticks(arange(34.85,35.06,0.05))
        ylim(2,6)
        yticks(range(2,7,1))
        title(str(dat.date[10*ii].values)[:10],fontsize=16)
        savefig(figdir+'MixedLayer/tenday/TS_PV10day_'+str(dat.date[10*ii].values)[:10]+'.png',bbox_inches='tight')

plotTSPVhex_10day()

XXXXXXXX


## Tried looking at these for signs of coincidence/preceeding wind events/ heat flux anomalies
## Want to try using ERA5 instead. These didn't quite cut it and Bob says NARR is no longer recommended.
#
# dat=xr.open_dataset('../data/aux_data/ERA_1804/tau_hflux_180411.nc')
#
# dat.history
# era=dat.rename({'inss':'tauy','iews':'taux','longitude':'lon','latitude':'lat','time':'date'})
# era['lon']=era['lon']-360
#
#
# narr=pickle.load(open('../pickles/wind/NARR_linet_newtheta.pickle','rb'))
# narr
#
# era['tau along']=era['taux']*cos(theta)+era['tauy']*sin(theta)
#
# def plot_narrfield(var):
#     figure(figsize=(12,3))
#     plot(narr.date,narr[var].mean(dim='distance'))
#     plot(era.date,era['tau along'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon')*30)
#     for ss in seltimes:
#             axvline(ss[0],color='k')
#             axvline(ss[1],color='k')
#     xlim(datetime.datetime(2015,8,1),datetime.datetime(2016,8,1))
#     title(var)
#
# plot_narrfield('along flow wind')
#
# era.lat
#
# def plot_erafield(var):
#     figure(figsize=(12,3))
#     plot(era.date,era[var].sel(lat=60).sel(lon=slice(-45,-25)).mean(dim='lon'))
#     for ss in seltimes:
#             axvline(ss[0],color='k')
#             axvline(ss[1],color='k')
#     xlim(datetime.datetime(2015,8,1),datetime.datetime(2016,8,1))
#     title(var)
#
# plot_erafield('taux')
#
# plot_erafield('tau along')
#
# plot_erafield('sshf')
#
# plot_erafield('slhf')
