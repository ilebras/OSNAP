from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'
# import cmocean

################################################################################################################################
###########################################    Load NorESM and obs    #######################################################
suffix='_sigma_UVcorr'
osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_18yrs_2004'+suffix+'.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_18yrs_2004'+suffix+'.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_18yrs_2004'+suffix+'.nc')

tot_trans=(osnap.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')-bso.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')-fs.TRANS.sum(dim='DZ').sum(dim='LONGITUDE'))
tot_trans.mean()

################################################################################################################################
################################################################################################################################
###########################################      GET overturning    #######################################################
################################################################################################################################
################################################################################################################################
denstep=0.01
denvec=arange(25.9,28.3,denstep)
def getpsi(east,dpdim):
    psimat=NaN*ones((len(denvec),len(east.TIME)))
    for ii,dd in enumerate(denvec):
            psimat[ii,:]=east['TRANS'].where(east['PDEN']>dd-denstep/2).where(east['PDEN']<=dd+denstep/2).sum(dim=dpdim).sum(dim='LONGITUDE').values
    east=east.assign_coords(DENLEV=(denvec))
    east['TRANSDEN']=((['DENLEV','TIME']),psimat)
    east['TRANSDEN_adj']=(east.TRANSDEN-(east.TRANSDEN.mean(dim='TIME').sum(dim='DENLEV')/len(east.DENLEV))).cumsum(dim='DENLEV')
    east['PSI']=east['TRANSDEN'].cumsum(dim='DENLEV')
    east['PSI_adj']=east['TRANSDEN'].cumsum(dim='DENLEV')
    east['SIGMAX']=(('TIME'),east.DENLEV[east.PSI_adj.argmax(dim='DENLEV').values])
    east['SIGMEAN']=east.DENLEV[east.PSI_adj.mean(dim='TIME').argmax(dim='DENLEV').values]
    # east['MOC']=east.PSI.max(dim='DENLEV')
    # east['MOCMEAN']=east.PSI[east.PSI.mean(dim='TIME').argmax(dim='DENLEV').values,:]

    return east

osnap=getpsi(osnap,'DZ')
bso=getpsi(bso,'DZ')
fs=getpsi(fs,'DZ')

# def plot_comp_psi(which,whichobs,name):
#     figure()
#     # plot(osnap_obs.PSI,osnap_obs.DENLEV,'C0',alpha=0.2,label='')
#     # plot(osnap.PSI,osnap.DENLEV,'C1',alpha=0.2,label='')
#     plot(whichobs.PSI.mean(dim='TIME'),whichobs.DENLEV,linewidth=3,label=name+' (2014-2016)',color='limegreen')
#     plot(which.PSI.mean(dim='TIME'),which.DENLEV,linewidth=3,label='NorESM (2010-2018)',color='purple')
#     if name=='OSNAP':
#         plot(which.PSI.sel(TIME=slice('2014-8-1','2016-9-1')).mean(dim='TIME'),which.DENLEV,linewidth=3,label='NorESM during OSNAP',color='purple',linestyle='--')
#     xlabel('Streamfunction [Sv]')
#     ylabel('pot. density anomaly [kg m$^{-3}$]')
#     legend(fontsize=16)
#     ylim(28.3,26)
#     savefig(figdir+name+'_PsiComp.png',bbox_inches='tight')
#     savefig(figdir+name+'_PsiComp.pdf',bbox_inches='tight')
#
# plot_comp_psi(osnap,osnap_obs,'OSNAP')
# plot_comp_psi(fs,fs_obs,'Fram Strait')
# plot_comp_psi(bso,bso_obs,'Barents Sea Opening')
#
# def MOC_comp():
#     osnap.MOC.plot()
#     osnap.MOCMEAN.plot()
#     osnap_obs.MOC.plot()
#     osnap_obs.MOCMEAN.plot()
#     xlim(datetime.datetime(2014,1,1),datetime.datetime(2017,1,1))
#     figure()
#     osnap.SIGMAX.plot()
#     axhline(osnap.SIGMEAN)
#     osnap_obs.SIGMAX.plot()
#     axhline(osnap_obs.SIGMEAN,color='C1')
#     xlim(datetime.datetime(2014,1,1),datetime.datetime(2017,1,1))
#     figure()
#     osnap.TRANS.sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
#     osnap_obs.TRANS.sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
#
# MOC_comp()

################################################################################################################################
################################################################################################################################
########################################      PLOT MEAN SECTION PANELS    #######################################################
################################################################################################################################
################################################################################################################################


# oslim=-42.4 #shelf only
oslim=-40 #just a little more
# oslim=-38 #generous

fslim=2
############################################################################################
################  NORESM  WATER MASS PARTITIONING   #############################################

AWS={}
PWS={}
DWS={}
AWN={}
PWN={}

# pws_depth=100
sigmax=osnap.SIGMEAN.values

sigmax
# xray=osnap
# for var in ['TRANS','PSAL','PTMP','PDEN']:
#     if var=='TRANS':
#         DWS[var]=xray['TRANS'].where(xray.LONGITUDE>=-42.4).where(xray.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
#         PWS[var]=xray['TRANS'].where(xray.LONGITUDE<-42.4).sum(dim='DZ').sum(dim='LONGITUDE')+xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
#         AWS[var]=xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
#
#     else:
#         DWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=-42.4).where(xray.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=-42.4).where(xray.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
#         PWS[var]=((xray[var]*xray['TRANS']).where(xray.LONGITUDE<-42.4).sum(dim='DZ').sum(dim='LONGITUDE')+(xray[var]*xray['TRANS']).where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE'))/PWS['TRANS']
#         AWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
#
xray=osnap
for var in ['TRANS','PSAL','PTMP','PDEN']:
    if var=='TRANS':
        DWS[var]=xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
        PWS[var]=xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
        AWS[var]=xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')

    else:
        DWS[var]=(xray[var]*xray['TRANS']).where(xray.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
        PWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE<oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
        AWS[var]=(xray[var]*xray['TRANS']).where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')/xray['TRANS'].where(xray.LONGITUDE>=oslim).where(xray.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')

#split up FS into PW and AW using definition in Tsubouchi et al. 2018, which is boundary between "EGC" + "Middle", 2W
# Here I'm using 4E, as it better maximizes transport...? Check on this...
#call all BSO water AW.

xray=fs
for var in ['TRANS','PSAL','PTMP','PDEN']:
    if var=='TRANS':
        PWN[var]=-fs['TRANS'].where(fs.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')
        AWN[var]=-fs['TRANS'].where(fs.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DZ').sum(dim='LONGITUDE')
    else:
        PWN[var]=-(fs[var]*fs['TRANS']).where(fs.PDEN<sigmax).sum(dim='DZ').sum(dim='LONGITUDE')/PWN['TRANS']
        AWN[var]=-((fs[var]*fs['TRANS']).where(fs.PDEN>=sigmax).sum(dim='DZ').sum(dim='LONGITUDE')+(bso[var]*bso['TRANS']).sum(dim='DZ').sum(dim='LONGITUDE'))/AWN['TRANS']


def plot_WMN():
    scatter(PWN['PSAL'].groupby('TIME.year').mean('TIME').values,PWN['PTMP'].groupby('TIME.year').mean('TIME').values,color='purple')
    scatter(AWN['PSAL'].groupby('TIME.year').mean('TIME').values,AWN['PTMP'].groupby('TIME.year').mean('TIME').values,color='orange')
    scatter(PWS['PSAL'].groupby('TIME.year').mean('TIME').values,PWS['PTMP'].groupby('TIME.year').mean('TIME').values,color='blue')
    scatter(AWS['PSAL'].groupby('TIME.year').mean('TIME').values,AWS['PTMP'].groupby('TIME.year').mean('TIME').values,color='red')
    scatter(DWS['PSAL'].groupby('TIME.year').mean('TIME').values,DWS['PTMP'].groupby('TIME.year').mean('TIME').values,color='grey')

plot_WMN()

PWN['TRANS'].mean()

PWS['TRANS'].mean()

PWN['TRANS'].plot()
PWS['TRANS'].plot()

##### Compare with diagnostics quickly
startyear=2000
startmonth=1
endyear=2018
endmonth=12

def load_allyears(varlab):
    tmp1=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_200001-200912.nc')[0])
    tmp2=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_201001-201812.nc')[0])
    varout=xr.concat([tmp1,tmp2],dim='time')
    varout=varout.rename({'time':'TIME'})
    varout['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return varout

hts=load_allyears('heattransports')
cp=3850
rhow=1025

(PWN['PTMP']*PWN['TRANS']+AWN['PTMP']*AWN['TRANS']).plot()
(-(hts.net_ht_BSO+hts.net_ht_FS)/cp/rhow/1e6).plot()


(AWS['PTMP']*AWS['TRANS']+DWS['PTMP']*DWS['TRANS']+PWS['PTMP']*PWS['TRANS']).plot()
((hts.net_ht_OSNAP)/cp/rhow/1e6).plot()


# ################################################################################################################################
# ################################################################################################################################
# ###########################################      SAVE WM props    #######################################################
# ################################################################################################################################
# ################################################################################################################################


WM={}
for ii,xrw in enumerate([AWS,PWS,DWS,AWN,PWN]):
    WM[ii]=xr.concat([xrw['PDEN'],xrw['PSAL'],xrw['PTMP'],xrw['TRANS']],pd.Index(['PDEN','PSAL','PTMP','TRANS'],name='var'))

WMall=xr.concat([WM[ww] for ww in WM],pd.Index(['AWS','PWS','DWS','AWN','PWN'],name='WM'))
WMall=WMall.to_dataset('var')

WMall.TRANS.sum(dim='WM').mean()

WMall.to_netcdf(datadir+'NorESM/NorESM_WMs_18yrs_2004_sigma_UVcorr.nc','w')
