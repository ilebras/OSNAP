from aux_funcs import *

datadir=datadir+'OSNAP2016recovery/'

datadir

#####################################################################
# Load all T and S into large dictionaries
#####################################################################
sal={}
tmp={}
date={}
prs={}
time={}
for moornum in range(1,8):
    sal[moornum]={}
    tmp[moornum]={}
    date[moornum]={}
    time[moornum]={}
    prs[moornum]={}

    #choose ctdlist for each mooring
    if moornum==7:
        ctdlist=glob.glob(datadir+'MCTD_Data_CF/MAT/CF'+str(moornum)+'*mat_ilebras.mat')
    else:
        ctdlist=glob.glob(datadir+'MCTD_Data_CF/NetCDF/*CF'+str(moornum)+'*.nc')
    if moornum==7:
        ctdlist=hstack((ctdlist,glob.glob(datadir+'RBR/CF'+str(moornum)+'*xr420*_ilebras.mat')))

    #load in each sal,tmp set
    for dd in ctdlist:
        if moornum==7:
            dat = io.loadmat(dd)
            tmp_hrly=array([float(tt) for tt in dat['temp'][:].flatten()])
            sal_hrly=array([float(ss) for ss in dat['psal'][:].flatten()])
            prs_hrly=array(dat['pres'][:].flatten())
            time_hrly=array(dat['dtnum'][:].flatten())
            date_hrly=array([datetime.datetime(1,1,1)+datetime.timedelta(days=tt-366) for tt in time_hrly])

        else:
            dat = Dataset(dd, 'r')
            time_hrly=array(dat.variables['TIME'][:])
            tmp_hrly=array([float(tt) for tt in dat.variables['TEMP'][:].flatten()])
            sal_hrly=array([float(ss) for ss in dat.variables['PSAL'][:].flatten()])
            prs_hrly=array(list(dat.variables['PRES'][:].flatten()))
            date_hrly=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=tt) for tt in time_hrly])


        prskey=int(round(nanmean(prs_hrly)/50.0)*50.0)
        timekey=mean(time_hrly)
        sal[moornum][prskey]=sal_hrly
        tmp[moornum][prskey]=tmp_hrly
        prs[moornum][prskey]=prs_hrly
        date[moornum][prskey]=date_hrly
        time[moornum][prskey]=time_hrly




########################################################################
##### Replace salinity with Jamie's and Kim's calibrated values ########
#######################################################################

KM=sort(glob.glob(datadir+'/MCTD_Data_CF/KMartini/*'))

sal[1][150]=io.loadmat(datadir+'MC_cal_JH_1810/CF1/cf1_171_2014_correctedsal.mat')['salc'].flatten()

sal[2][100]=io.loadmat(datadir+'MC_cal_JH_1810/CF2/cf2_100_2014_correctedsal.mat')['salc'].flatten()
sal[2][200]=io.loadmat(datadir+'MC_cal_JH_1810/CF2/cf2_176_2014_correctedsal.mat')['salc'].flatten()

sal[3][200]=io.loadmat(datadir+'MC_cal_JH_1810/CF3/cf3_181_2014_correctedsal.mat')['salc'].flatten()

sal[4][100]=xr.open_dataset(KM[0])['s0'].values
sal[4][200]=io.loadmat(datadir+'MC_cal_JH_1810/CF4/cf4_200_2014_correctedsal.mat')['salc'].flatten()
sal[4][350]=io.loadmat(datadir+'MC_cal_JH_1810/CF4/cf4_350_2014_correctedsal.mat')['salc'].flatten()
sal[4][400]=io.loadmat(datadir+'MC_cal_JH_1810/CF4/cf4_500_2014_correctedsal.mat')['salc'].flatten()

sal[5][100]=io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_100_2014_correctedsal.mat')['salc'].flatten()
sal[5][250]=io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_250_2014_correctedsal.mat')['salc'].flatten()
sal[5][500]=io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_500_2014_correctedsal.mat')['salc'].flatten()
sal[5][750]=io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_750_2014_correctedsal.mat')['salc'].flatten()
sal[5][1000]=io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_1000_2014_correctedsal.mat')['salc'].flatten()
sal[5][1300]=io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_1275_2014_correctedsal.mat')['salc'].flatten()

#remove the spike at CF5 250m
plot(sal[5][250])
date5=datetime.datetime(1950,1,1)+datetime.timedelta(days=int(24305))
sal[5][250][date[5][250]>date5]=nan

sal[6][150]=xr.open_dataset(KM[1])['s0'].values[:-4]
sal[6][300]=io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_250_2014_correctedsal.mat')['salc'].flatten()
sal[6][550]=io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_500_2014_correctedsal.mat')['salc'].flatten()
sal[6][800]=io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_750_2014_correctedsal.mat')['salc'].flatten()
sal[6][1050]=io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_1000_2014_correctedsal.mat')['salc'].flatten()
sal[6][1550]=io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_1500_2014_correctedsal.mat')['salc'].flatten()
sal[6][1850]=io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_1830_2014_correctedsal.mat')['salc'].flatten()

# sal[7][50]=sal[7][50]-0.025 could re-implement this...
sal[7][100]=io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_100_2014_correctedsal.mat')['salc'].flatten()
sal[7][250]=io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_250_2014_correctedsal.mat')['salc'].flatten()
sal[7][500][:]=NaN
sal[7][750]=io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_750_2014_correctedsal.mat')['salc'].flatten()
sal[7][1000]=io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_1000_2014_correctedsal.mat')['salc'].flatten()
sal[7][1500]=io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_1500_2014_correctedsal.mat')['salc'].flatten()
sal[7][1900]=io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_1900_2014_correctedsal.mat')['salc'].flatten()

###################################################################
################## Some useful functions #########################
###################################################################
# set up basic ts plotting framework

def TSplot(moor,prs,alph=0.3):
    plot(sal[moor][prs],tmp[moor][prs],'o',
         markersize=5,label='CF'+str(moor)+', '+str(prs)+' db',
         alpha=alph,mew=0)

def TSplotafter(moor,prs,alph=0.3):
    plot(salfinal[moor][prs],tmp[moor][prs],'o',
         markersize=5,label='CF'+str(moor)+', '+str(prs)+' db',
         alpha=alph,mew=0)

def spruce_TS():
    dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k')
    ylabel('potential temperature ($^\circ\ C$)')
    xlabel('salinity')
    legend(numpoints=1,markerscale=2)


def plotnfilt4(what):
    figure(figsize=(14,3))
    divit=30
    plot(date[4][350],what)
    plot(date[4][350][::divit],hrly_ave(what,divit),'ro-')
    plot(date[4][350][::divit*6],hrly_ave(what,divit*6),'yo-')


pden={}
ptmp={}
for mm in range(1,8):
    pden[mm]={}
    ptmp[mm]={}
    for key in sort(list(sal[mm].keys())):
        SAvec=gsw.SA_from_SP(sal[mm][key],prs[mm][key],CFlon[mm-1],CFlat[mm-1])
        ptmp[mm][key]=gsw.pt0_from_t(SAvec,tmp[mm][key],prs[mm][key])
        pden[mm][key]=gsw.sigma0(SAvec,gsw.CT_from_pt(SAvec,ptmp[mm][key]))


# for mm in range(2,8):
#     prslist=list(sort(list(sal[mm].keys())))
#     # if mm==3:
#     #     prslist.remove(100)
#     if mm==7:
#         prslist.remove(500)
#     for ii,kk in enumerate(prslist[:-1]):
#         if len(pden_precorr[mm][kk])>len(pden_precorr[mm][prslist[ii+1]]):
#             lesslen=len(pden_precorr[mm][prslist[ii+1]])
#         else:
#             lesslen=len(pden_precorr[mm][kk])
#         figure(figsize=(12,3))
#         plot(date[mm][kk][:lesslen],pden_precorr[mm][prslist[ii+1]][:lesslen]-pden_precorr[mm][kk][:lesslen])
#         axhline(0)
#         # plot(date[mm][kk][:lesslen],pden[mm][prslist[ii+1]][:lesslen]-pden[mm][kk][:lesslen])
#         title('CF'+str(mm)+': '+str(prslist[ii+1])+'db-'+str(kk)+'db; '+str(sum((pden_precorr[mm][prslist[ii+1]][:lesslen]-pden_precorr[mm][kk][:lesslen])<0))+' inversions')
#         ylabel('Density difference (kg m$^{-3}$)')
        # savefig('../figures/salcalib_15min/JH/CF'+str(mm)+'_'+str(prslist[ii+1])+'db-'+str(kk)+'db_JH.png',bbox_inches='tight')

#####################################################################
######### Plot salinity time series before and after ################
#####################################################################

for mm in range(1,8):
    figure(figsize=(14,3))
    for key in sort(list(sal[mm].keys())):
        plot(date[mm][key],sal[mm][key],label=str(key)+' db')
    title('Salinity at CF'+str(mm)+' dip calibrated')
    ylabel('Salinity')
    legend()
    # savefig('../figures/salcalib_15min/JH/CF'+str(mm)+'_sal_JH.png',bbox_inches='tight')

#####################################################################
######### Plot TS space before and after ############################
#####################################################################

#
# for mm in range(1,8):
#     figure()
#     for key in sort(list(sal[mm].keys())):
#         TSplot(mm,key)
#     spruce_TS()
#     title('CF'+str(mm)+': $\Theta$-S space dip calibrated')
#     if mm>=6:
#         xlim([34,35.2])
    # savefig('../figures/salcalib_15min/JH/CF'+str(mm)+'_TS_JH.png',bbox_inches='tight')
month={}
for key in date:
    month[key]={}
    for key2 in date[key]:
        month[key][key2]=[dd.month for dd in date[key][key2]]


#just to be consistent with nomenclature as the below is copied from another script.
#(want to save as netcdf from this file directly and circumvent pickle nonsense.)
date_all=date.copy()
month_all=month.copy()
prs_all=prs.copy()
sal_all=sal.copy()
tmp_all=tmp.copy()

#confirm that length/start date lines up
for ii in range(1,8):
    for kk in prs_all[ii]:
        print(ii,shape(prs_all[ii][kk]),date_all[ii][kk][0],date_all[ii][kk][-1])
#lots of different time starts here... will have to rethink strategy

#make an xarray for each instrument.
#then resample at 1H, and merge to get one for each mooring
CF_16_xrays={}
for ii in range(2,8):
    dpths=list(date_all[ii].keys())
    CF_16_xrays[ii]={}
    for dd in dpths:
        CF_16_xrays[ii][dd]=xr.Dataset({'PRES': (['TIME'], prs_all[ii][dd] ),
                                        'PSAL': (['TIME'], sal_all[ii][dd]),
                                        'PTMP': (['TIME'],  tmp_all[ii][dd]),},
                                        coords={'TIME': date_all[ii][dd],})


def add_SA_CT_PT(xray):
    if 'PRES' in list(xray.data_vars):
            PRES_out=xray['PRES']
    else:
            PRES_out=-gsw.p_from_z(xray['DEPTH'])
    SA_out=gsw.SA_from_SP(xray['PSAL'],PRES_out,xray.LONGITUDE,xray.LATITUDE)
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['ASAL']=(('TIME','DEPTH'),SA_out)
    xray['PTMP']=(('TIME','DEPTH'),PT_out)
    xray['CTMP']=(('TIME','DEPTH'),CT_out)
    xray['PDEN']=(('TIME','DEPTH'),PD_out)





#note: not going to do CF1... because Astrid doesn't need it.
CF_16_hourly={}
for ii in range(2,8):
    dpths=list(date_all[ii].keys())
    CF_16_hourly[ii]=xr.concat([CF_16_xrays[ii][dd] for dd in sort(dpths)],dim='DEPTH')
    CF_16_hourly[ii]=CF_16_hourly[ii].assign_coords({'DEPTH':sort(dpths),
                                                     'LATITUDE': CFlat[ii-1],
                                                     'LONGITUDE': CFlon[ii-1]})
    CF_16_hourly[ii]=CF_16_hourly[ii].resample(TIME='1H').mean()
    add_SA_CT_PT(CF_16_hourly[ii])
    CF_16_hourly[ii].to_netcdf(datadir+'Hourly_netcdf/CF'+str(ii)+'_mcat_2016recovery_hourly.nc','w',format='netCDF4')


# pickle.dump([date,month,prs,sal,ptmp,pden],
#             open(datadir+'/pickles/TSdailydic/TS_15min_dic_wJHdipcal.pickle','wb'))

XXXXXXXXXXX
#
# #####################################################################
# ######### Plot TS space correspondence between moorings #############
# #####################################################################
#
# dlim=[0,100,300,750,2000]
#
# for ii,dd in enumerate(dlim[:-1]):
#
#     figure()
#     for mm in range(4,8):
#         for key in sort(list(sal[mm].keys())):
#             if (key>dd) & (key<=dlim[ii+1]):
#                 TSplotafter(mm,key)
#     spruce_TS()
#     # legend(loc=(1.05,0.2))
#     title('CF4-7: '+str(dd)+'db-'+str(dlim[ii+1])+'db $\Theta$-S space AFTER corrections')
#     if dd>=100:
#         xlim([33,35.5])
#     if dd>=300:
#         xlim([34,35.2])
#     if dd>=500:
#         xlim([34.8,35.1])
#     savefig('../figures/salcalib_15min/CF4-7_TSafter_'+str(dd)+'-'+str(dlim[ii+1])+savename+'.png',bbox_inches='tight')
#
#     figure()
#     for mm in range(4,8):
#         for key in sort(list(sal[mm].keys())):
#             if (key>dd) & (key<=dlim[ii+1]):
#                 TSplot(mm,key)
#     spruce_TS()
#     # legend(loc=(1.05,0.2))
#     title('CF4-7: '+str(dd)+'db-'+str(dlim[ii+1])+'db $\Theta$-S space BEFORE corrections')
#     if dd>=100:
#         xlim([33,35.5])
#     if dd>=300:
#         xlim([34,35.2])
#     if dd>=500:
#         xlim([34.8,35.1])
#     savefig('../figures/salcalib_15min/CF4-7_TSbefore'+str(dd)+'-'+str(dlim[ii+1])+savename+'.png',bbox_inches='tight')
#
#
# shelflim=200
#
# figure()
# for mm in range(1,5):
#     for key in sort(list(sal[mm].keys())):
#         if (key<=shelflim):
#             TSplotafter(mm,key)
# spruce_TS()
# # legend(loc=(1.05,0.2))
# title('CF1-4 $\Theta$-S space AFTER corrections')
# savefig('../figures/salcalib_15min/CF1-4_TSafter'+savename+'.png',bbox_inches='tight')
#
# figure()
# for mm in range(1,5):
#     for key in sort(list(sal[mm].keys())):
#         if (key<=shelflim):
#             TSplot(mm,key)
# spruce_TS()
# # legend(loc=(1.05,0.2))
# title('CF1-4 $\Theta$-S space BEFORE corrections')
# savefig('../figures/salcalib_15min/CF1-4_TSbefore'+savename+'.png',bbox_inches='tight')
#
# for key in date:
#     month[key]={}
#     for key2 in date[key]:
#         month[key][key2]=[dd.month for dd in date[key][key2]]

#
# max_sdiff={}
# mean_sdiff={}
# max_pdendiff={}
# mean_pdendiff={}
# for key in sal:
#     if 'corr' in str(key):
#         max_sdiff[key]={}
#         mean_sdiff[key]={}
#         max_pdendiff[key]={}
#         mean_pdendiff[key]={}
#         for dpkey in sal[key]:
#             max_sdiff[key][dpkey]=round(nanmax(abs(sal[key][dpkey]-sal[int(key[0])][dpkey])),3)
#             mean_sdiff[key][dpkey]=round(nanmean(sal[key][dpkey]-sal[int(key[0])][dpkey]),3)
#             max_pdendiff[key][dpkey]=round(nanmax(abs(pden[int(key[0])][dpkey]-pden_precorr[int(key[0])][dpkey])),3)
#             mean_pdendiff[key][dpkey]=round(nanmean(pden[int(key[0])][dpkey]-pden_precorr[int(key[0])][dpkey]),3)
#
#
# max_pdendiff
