from aux_funcs import *

#####################################################################
# Load all T and S into large dictionaries
#####################################################################

datadir

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
        ctdlist=glob.glob(datadir+'OSNAP2016recovery/MCTD_Data_CF/MAT/CF'+str(moornum)+'*mat_ilebras.mat')
    else:
        ctdlist=glob.glob(datadir+'OSNAP2016recovery/MCTD_Data_CF/NetCDF/*CF'+str(moornum)+'*.nc')
    if moornum==7:
        ctdlist=hstack((ctdlist,glob.glob(datadir+'OSNAP2016recovery/RBR/CF'+str(moornum)+'*xr420*_ilebras.mat')))
    #load in each sal,tmp set
    for dd in ctdlist:
        if moornum==7:
            dat = io.loadmat(dd)
            tmp_hrly=hrly_ave([float(tt) for tt in dat['temp'][:].flatten()],aveconst)
            sal_hrly=hrly_ave([float(ss) for ss in dat['psal'][:].flatten()],aveconst)
            prs_hrly=hrly_ave(list(dat['pres'][:].flatten()),aveconst)
            time_hrly=dat['dtnum'][:].flatten()[::aveconst][:len(prs_hrly)]
            date_hrly=array([datetime.datetime(1,1,1)+datetime.timedelta(days=int(tt-366)) for tt in time_hrly])

        else:
            dat = Dataset(dd, 'r')
            time_hrly=array(dat.variables['TIME'][:])[::aveconst][:-1]
            tmp_hrly=hrly_ave([float(tt) for tt in dat.variables['TEMP'][:].flatten()],aveconst)[:len(time_hrly)]
            sal_hrly=hrly_ave([float(ss) for ss in dat.variables['PSAL'][:].flatten()],aveconst)[:len(time_hrly)]
            prs_hrly=hrly_ave(list(dat.variables['PRES'][:].flatten()),aveconst)[:len(time_hrly)]
            date_hrly=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=int(tt)) for tt in time_hrly])

        prskey=int(round(nanmean(prs_hrly)/50.0)*50.0)
        timekey=mean(time_hrly)

        sal[moornum][prskey]=sal_hrly
        tmp[moornum][prskey]=gsw.pt0_from_t(sal_hrly,tmp_hrly,prs_hrly)
        prs[moornum][prskey]=prs_hrly
        date[moornum][prskey]=date_hrly
        time[moornum][prskey]=time_hrly

figure(figsize=(12,3))
plot(date[7][500],tmp[7][500])

###################################################################
#####Replace salinity with (daily averaged!) dip calibrated data ##########
###################################################################

KM=sort(glob.glob(datadir+'/MCTD_Data_CF/KMartini/*'))

sal[1][150]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF1/cf1_171_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[1][150])]

sal[2][100]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF2/cf2_100_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[2][100])]
sal[2][200]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF2/cf2_176_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[2][200])]

sal[3][200]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF3/cf3_181_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[3][200])]

sal[4][100]=hrly_ave(xr.open_dataset(KM[0])['s0'].values,aveconst)[:len(time[4][100])]
sal[4][200]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF4/cf4_200_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[4][200])]
sal[4][350]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF4/cf4_350_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[4][350])]
sal[4][400]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF4/cf4_500_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[4][400])]

#add a correction to 400m instrument
sal[4][400]=sal[4][400]+0.02

sal[5][100]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_100_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[5][100])]
sal[5][250]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_250_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[5][250])]
sal[5][500]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_500_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[5][500])]
sal[5][750]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_750_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[5][750])]
sal[5][1000]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_1000_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[5][1000])]
sal[5][1300]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF5/cf5_1275_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[5][1300])]

#remove the spike at CF5 250m

date5=datetime.datetime(1950,1,1)+datetime.timedelta(days=int(24305))
sal[5][250][date[5][250]>date5]=nan

def removeline(tvec,svec,anchor=0):
    linefunc=poly1d(polyfit(tvec,svec, 1))
    liney=linefunc(tvec)
    corrected=svec-liney+liney[anchor]
    return corrected




# sal[6][150]=hrly_ave(xr.open_dataset(KM[1])['s0'].values[:-4],aveconst)[:len(time[6][150])]
#replacing with my calibration:
sal[6][150]=removeline(time[6][150],sal[6][150])-0.0122
sal[6][100]=sal[6][100]-0.024

sal[6][300]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_250_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[6][300])]
sal[6][550]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_500_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[6][550])]
sal[6][800]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_750_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[6][800])]
sal[6][1050]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_1000_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[6][1050])]
sal[6][1550]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_1500_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[6][1550])]
sal[6][1850]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF6/cf6_1830_2014_correctedsal.mat')['salc'].flatten(),aveconst)[:len(time[6][1850])]

sal[7][100]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_100_2014_correctedsal.mat')['salc'].flatten(),aveconst)
sal[7][250]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_250_2014_correctedsal.mat')['salc'].flatten(),aveconst)
sal[7][500][:]=NaN
sal[7][750]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_750_2014_correctedsal.mat')['salc'].flatten(),aveconst)
sal[7][1000]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_1000_2014_correctedsal.mat')['salc'].flatten(),aveconst)
sal[7][1500]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_1500_2014_correctedsal.mat')['salc'].flatten(),aveconst)
sal[7][1900]=hrly_ave(io.loadmat(datadir+'MC_cal_JH_1810/CF7/cf7_1900_2014_correctedsal.mat')['salc'].flatten(),aveconst)
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


pden_precorr={}

for mm in range(1,8):
    print(mm)
    pden_precorr[mm]={}
    for key in sort(list(sal[mm].keys())):
        print(key)
        pden_precorr[mm][key]=gsw.sigma0(sal[mm][key],gsw.CT_from_pt(sal[mm][key],tmp[mm][key]))

for mm in range(2,8):
    prslist=list(sort(list(sal[mm].keys())))
    # if mm==3:
    #     prslist.remove(100)
    if mm==7:
        prslist.remove(500)
    for ii,kk in enumerate(prslist[:-1]):
        if len(pden_precorr[mm][kk])>len(pden_precorr[mm][prslist[ii+1]]):
            lesslen=len(pden_precorr[mm][prslist[ii+1]])
        else:
            lesslen=len(pden_precorr[mm][kk])
        figure(figsize=(12,3))
        plot(date[mm][kk][:lesslen],pden_precorr[mm][prslist[ii+1]][:lesslen]-pden_precorr[mm][kk][:lesslen])
        axhline(0)
        # plot(date[mm][kk][:lesslen],pden[mm][prslist[ii+1]][:lesslen]-pden[mm][kk][:lesslen])
        title('CF'+str(mm)+': '+str(prslist[ii+1])+'db-'+str(kk)+'db; '+str(sum((pden_precorr[mm][prslist[ii+1]][:lesslen]-pden_precorr[mm][kk][:lesslen])<0))+' inversions')
        ylabel('Density difference (kg m$^{-3}$)')
        # savefig('../figures/salcalib/JH/CF'+str(mm)+'_'+str(prslist[ii+1])+'db-'+str(kk)+'db_JH.png',bbox_inches='tight')

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


for mm in range(1,8):
    figure()
    for key in sort(list(sal[mm].keys())):
        TSplot(mm,key)
    spruce_TS()
    title('CF'+str(mm)+': $\Theta$-S space dip calibrated')
    if mm>=6:
        xlim([34,35.2])
    # savefig('../figures/salcalib_15min/JH/CF'+str(mm)+'_TS_JH.png',bbox_inches='tight')
month={}
for key in date:
    month[key]={}
    for key2 in date[key]:
        month[key][key2]=[dd.month for dd in date[key][key2]]


pickle.dump([date,month,prs,sal,tmp],
            open('../pickles/TSdailydic/TS_daily_dic_wJHIL.pickle','wb'))

XXXXXXXXXXX
