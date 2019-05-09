from aux_funcs import *

datadir

savename='_15min'

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
            date_hrly=array([datetime.datetime(1,1,1)+datetime.timedelta(days=int(tt-366)) for tt in time_hrly])

        else:
            dat = Dataset(dd, 'r')
            time_hrly=array(dat.variables['TIME'][:])
            tmp_hrly=array([float(tt) for tt in dat.variables['TEMP'][:].flatten()])
            sal_hrly=array([float(ss) for ss in dat.variables['PSAL'][:].flatten()])
            prs_hrly=array(list(dat.variables['PRES'][:].flatten()))
            date_hrly=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=int(tt)) for tt in time_hrly])


        prskey=int(round(nanmean(prs_hrly)/50.0)*50.0)
        timekey=mean(time_hrly)
#         if moornum==8:
#             sal[moornum][str(prskey)+'_'+str(timekey)]=sal_hrly
#             tmp[moornum][str(prskey)+'_'+str(timekey)]=tmp_hrly
#             prs[moornum][str(prskey)+'_'+str(timekey)]=prs_hrly
#             date[moornum][str(prskey)+'_'+str(timekey)]=date_hrly
#             time[moornum][str(prskey)+'_'+str(timekey)]=time_hrly
#         else:
        sal[moornum][prskey]=sal_hrly
        tmp[moornum][prskey]=tmp_hrly
        prs[moornum][prskey]=prs_hrly
        date[moornum][prskey]=date_hrly
        time[moornum][prskey]=time_hrly




###################################################################
################## Some useful functions #########################
###################################################################


def removeline(tvec,svec,anchor=0):
    nanind=~isnan(svec)
    linefunc=poly1d(polyfit(tvec[nanind],svec[nanind], 1))
    liney=svec.copy()
    liney[nanind]=linefunc(tvec[nanind])
    corrected=svec.copy()
    corrected[nanind]=svec[nanind]-liney[nanind]+liney[anchor]
    return liney,corrected

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


# def TSplot_seasonal(moor,prs):
#     scatter(sal[moor][prs],tmp[moor][prs],c=month[moor][prs],
#             cmap=cm.brg_r,lw=0,alpha=0.5,vmin=1,vmax=12)
#
# def spruce_seasonal():
#     cbar=colorbar(ticks=range(1,13,2))
#     cbar.set_ticklabels(['January','March','May','July','September','November'])
#     ylabel('potential temperature ($^\circ\ C$)')
#     xlabel('salinity')
#     dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k')
#     clabel(dens)


def plotnfilt4(what):
    figure(figsize=(14,3))
    divit=30
    plot(date[4][350],what)
    plot(date[4][350][::divit],hrly_ave(what,divit),'ro-')
    plot(date[4][350][::divit*6],hrly_ave(what,divit*6),'yo-')


## Some code is hidden below for a rainy day

#####################################################################
### Inspect CF1 TS space and how it relates to CF2

#####################################################################
#
# for key in sal[1]:
#     TSplot(1,key)
#
# TSplot(2,50)
# TSplot(2,100)
# spruce_TS()
#
#
# for key in sal['1recon']:
#     TSplot('1recon',key)
# TSplot(1,150)
# TSplot(2,50)
# TSplot(2,100)
#
# spruce_TS()
#
#
# for key in sal['1recon']:
#     TSplot_seasonal('1recon',key)
# spruce_seasonal()
#
#
# for key in sal[1]:
#     TSplot_seasonal(1,key)
# for key in sal[2]:
#     TSplot_seasonal(2,key)
# spruce_seasonal()
#
#
# for key in sal[2]:
#     TSplot_seasonal(2,key)
# spruce_seasonal()
#
#
# TSplot_seasonal(1,50)
# TSplot_seasonal(2,50)
# spruce_seasonal()
#
#
# TSplot_seasonal('1recon',50)
# TSplot_seasonal(2,50)
# spruce_seasonal()
#
#
# TSplot_seasonal(1,100)
# TSplot_seasonal(2,100)
# spruce_seasonal()
#
#
# TSplot_seasonal('1recon',100)
# TSplot_seasonal(2,100)
# spruce_seasonal()

#####################################################################
#### Look at the end of record for CF1 and CF4 and compare to neighbors

#####################################################################
# end4=datetime.datetime(1950,1,1)+datetime.timedelta(days=valitime)
#
#
#
#
# const=100
#
#
#
#
# sal[2].keys()
#
#
#
# figure(figsize=(12,4))
# keyvec=[150,200]
# for moor in range(1,4):
#     for key in keyvec:
#         if key in sal[moor].keys():
#             plot(date[moor][key][date[moor][key]>datetime.datetime(2016, 4, 1)],
#              sal[moor][key][date[moor][key]>datetime.datetime(2016, 4,1)],
#             label='CF'+str(moor)+' '+str(key)+'db')
# axvline(end4,color='k')
# legend(loc=(1.05,0.3))
# savefig('../figures/salcorr/CF1_salatend.png',bbox_inches='tight')
#
#
# colorvec=['purple','grey','orange']
#
# figure(figsize=(12,4))
# keyvec=[150,200]
# ii=0
# for moor in range(1,4):
#     for key in keyvec:
#         if key in sal[moor].keys():
#             plot(sal[moor][key],
#              tmp[moor][key],'.',color=colorvec[ii],
#             label='CF'+str(moor)+' '+str(key)+'db',
#                 alpha=0.5)
#
#             plot(sal[moor][key][date[moor][key]>end4],
#              tmp[moor][key][date[moor][key]>end4],'o',color=colorvec[ii],
#             label='CF'+str(moor)+' '+str(key)+'db',zorder=10)
#     ii+=1
# legend(loc=(1.05,0.3),numpoints=1)
# savefig('../figures/salcorr/CF1_salatend_TS.png',bbox_inches='tight')
#
#
#
# figure(figsize=(12,4))
# keyvec=[50,100]
# for moor in range(3,6):
#     for key in keyvec:
#         if key in sal[moor].keys():
#             plot(date[moor][key][date[moor][key]>datetime.datetime(2016, 4, 1)],
#              sal[moor][key][date[moor][key]>datetime.datetime(2016, 4,1)],
#             label='CF'+str(moor)+' '+str(key)+'db')
# axvline(end4,color='k')
# legend(loc=(1.05,0.3))
# savefig('../figures/salcorr/CF4_salatend.png',bbox_inches='tight')
#
#
# figure(figsize=(12,4))
# keyvec=[50,100]
# for moor in range(3,6):
#     for key in keyvec:
#         if key in sal[moor].keys():
#             plot(sal[moor][key][date[moor][key]>end4],
#              tmp[moor][key][date[moor][key]>end4],'o',
#             label='CF'+str(moor)+' '+str(key)+'db')
# legend(loc=(1.05,0.3),numpoints=1)
# savefig('../figures/salcorr/CF4_salatend_TS.png',bbox_inches='tight')
# savefig('../figures/salcorr/CF4_salatend_TS.png',bbox_inches='tight')

#####################################################################
#####################################################################
#####################################################################
###########                 CF1                     #################
#####################################################################
#####################################################################
#####################################################################
liney={}

liney[1]={}
sal['1corr']={}
liney[1][50],sal['1corr'][50]=removeline(time[1][50],sal[1][50])

sal['1corr'][100]=sal[1][100].copy()
sal['1corr'][100][-2:]=NaN

sal['1corr'][150]=sal[1][150].copy()


#####################################################################
#####################################################################
#####################################################################
###########                 CF2                     #################
#####################################################################
#####################################################################
#####################################################################

liney[2]={}
sal['2corr']={}
for key in sal[2]:
    liney[2][key],sal['2corr'][key]=removeline(time[2][key],sal[2][key])


sal['2corr'][50]=sal['2corr'][50]-0.065
sal['2corr'][100]=sal['2corr'][100]
sal['2corr'][200]=sal['2corr'][200]+0.01


#####################################################################
#####################################################################
#####################################################################
###########                 CF3                     #################
#####################################################################
#####################################################################
#####################################################################

# Going to leave out 100db instrument's salinity as it is buggy and not long -- not worth keeping around.
sal['3corr']={}
sal['3corr'][100]=sal[3][100].copy()
# sal['3corr'][100][:]=NaN #Actually, leave it in!!

#####################################################################
#####################################################################
#####################################################################
###########                 CF4                     #################
#####################################################################
#####################################################################
#####################################################################


#####################################################################
###### Drift spotted in CF4 @ 100db #########################
#####################################################################


driftstartdate=datetime.datetime(2015,11,1)
driftenddate=datetime.datetime(2016,7,1)

# detrend 50db instrument first
sal['4corr']={}
liney[4]={}
liney[4][50],sal['4corr'][50]=removeline(time[4][50],sal[4][50])

saldif4=sal[4][100]-sal['4corr'][50]


# Then 100db inst
timeind4=(date[4][100]>driftstartdate) & (date[4][100]<driftenddate)
salcorr4=saldif4[timeind4]
tcorr4=time[4][100][timeind4]
liney4,corr4=removeline(tcorr4,salcorr4)

sal['4corr'][100]=sal[4][100].copy()
sal['4corr'][100][date[4][100]>driftenddate]=NaN
sal['4corr'][100][timeind4]=corr4+sal[4][50][timeind4]
sal['4corr'][100]=sal['4corr'][100]

# 100db inst fig

figure(figsize=(14,3))
plot(date[4][50],saldif4,label='Original difference')
axvline(driftstartdate,color='k')
axvline(driftenddate,color='k')
plot(date[4][100][timeind4],liney4,'b')
plot(date[4][50][timeind4],corr4,label='Corrected difference')
grid('on')
title('CF4 Salinity$_{100\mathrm{db}}$- Salinity$_{50\mathrm{db}}$')
ylabel('salinity difference')
legend()
savefig('../figures/salcalib_15min/CF4_100_scorr'+savename+'.png',bbox_inches='tight')

#####################################################################
###### Another drift at CF4 @ 350db #########################
#####################################################################

#Remove linear fit to salinity at 350db

# 350 inst fig (and detrending/adjusting)
figure(figsize=(14,3))
plot(date[4][350],sal[4][350],label='Original record')
liney[4][350],sal['4corr'][350]=removeline(time[4][350],sal[4][350])
plot(date[4][350],liney[4][350],'b')
plot(date[4][350],sal['4corr'][350],label='Detrended record')
axvline(datetime.datetime(2016,1,1),color='k')
sal['4corr'][350][date[4][350]>datetime.datetime(2016,1,1)]=sal['4corr'][350][date[4][350]>datetime.datetime(2016,1,1)]+0.025
plot(date[4][350],sal['4corr'][350],label='Final record')
title('Salinity time series at 350db instrument on CF4')
ylabel('salinity')
legend()
grid('on')
savefig('../figures/salcalib_15min/CF4_350_scorr_'+savename+'.png',bbox_inches='tight')

#Remove linear fit to difference between 400 and 350db instruments
dsplit400=datetime.datetime(2015,10,1)
dateind400=[date[4][350]<dsplit400]
saldiff434=(sal[4][400]-sal['4corr'][350])[dateind400]
liney[4][400],saldiff434_corr=removeline(time[4][350][dateind400],saldiff434)
newdiff434=(sal[4][400]-sal['4corr'][350]).copy()
newdiff434[dateind400]=saldiff434_corr
newdiff434[date[4][350]>=dsplit400]=(sal[4][400]-sal['4corr'][350])[date[4][350]>=dsplit400]-0.9*liney[4][400][-1]

sal['4corr'][400]=sal['4corr'][350]+newdiff434+0.022

tester=sal['4corr'][400].copy()

# 400db inst fig
figure(figsize=(14,3))
plot(date[4][350],sal[4][400]-sal['4corr'][350],label='Original difference')
plot(date[4][350][dateind400],liney[4][400],'b')
plot(date[4][350],newdiff434,label='Detrended difference')
plot(date[4][350],newdiff434+0.022,label='Final difference')
title('CF4: Salinity$_{400\mathrm{db}}$- Salinity$_{350\mathrm{db}}$')
axvline(dsplit400,color='k')
axhline(0,color='k')
grid('on')
legend()
ylabel('salinity difference')
savefig('../figures/salcalib_15min/CF4_400_scorr_'+savename+'.png',bbox_inches='tight')
#####################################################################
#####################################################################
#####################################################################
###########                 CF5                    #################
#####################################################################
#####################################################################
#####################################################################

#####################################################################
### Examine an apparent spike at the end of 250m instrument on CF5
#####################################################################

#spike is after this date
date5=datetime.datetime(1950,1,1)+datetime.timedelta(days=int(24305))

sal['5corr']={}
tmp['5corr']={}
sal['5corr'][250]=sal[5][250].copy()
tmp['5corr'][250]=tmp[5][250].copy()
sal['5corr'][250][date[5][250]>date5]=nan


#####################################################################
### Density inversions/drift on CF5, 500-1300db
#####################################################################


# Remove linear trend in all records

liney[5]={}
for key in sal[5]:
    if key!=250:
        liney[5][key],sal['5corr'][key]=removeline(time[5][key],sal[5][key])
    else:
        sal5_250=sal['5corr'][250].copy()
        liney525,sal['5corr'][250][~isnan(sal5_250)]=removeline(time[5][250][~isnan(sal5_250)],sal5_250[~isnan(sal5_250)])

sal['5corr'][1300]=sal['5corr'][1300]+0.003



#####################################################################
#####################################################################
#####################################################################
###########                 CF6                    #################
#####################################################################
#####################################################################
#####################################################################

#####################################################################
### Next, deal with a giant drift in salinity at 150m at CF6 #####
#####################################################################


liney[6]={}
sal['6corr']={}

liney[6][150],sal['6corr'][150]=removeline(time[6][150],sal[6][150])

sal['6corr'][150]=sal['6corr'][150]-0.0122
sal['6corr'][100]=sal[6][100]-0.024

figure(figsize=(14,3))
plot(date[6][150],sal[6][150],label='Original record')
plot(date[6][150],liney[6][150],'b')
plot(date[6][150],sal['6corr'][150]+0.0122,label='Detrended record')
plot(date[6][150],sal['6corr'][150],label='Final record')
title('Salinity time series at 150db instrument on CF6')
ylabel('salinity')
legend()
grid('on')
savefig('../figures/salcalib_15min/CF6_150_scorr'+savename+'.png',bbox_inches='tight')


#### Some valuable thoughts on TS space below, keeping around for now.
#
#
# tmp[str(key)+'5corr']={}
# for key2 in sal['5corr']:
#         tmp['5corr'][key2]=tmp[5][key2]
#
#
#
# sal[5].keys()
#
# sal[7].keys()
#
# figure()
# # TSplot(6,300)
# # TSplot(7,250)
# TSplot('6corr',150)
# TSplot('6corr',100)
# # TSplot('5corr',100)
# TSplot(7,50)
# spruce_TS()
# xlim([34.4,35.2])
#
#
#
#
# figure()
# TSplot('6corr',150)
# TSplot('6corr',100)
# # TSplot('5corr',100)
# TSplot(7,100)
#
# figure()
# TSplot(6,100)
# TSplot('6corr',150)
# # TSplot('5corr',100)
# TSplot(7,100)
#
# figure()
# TSplot(6,300)
# TSplot(6,550)
#
# sal[7].keys()
#
# figure()
# TSplot(6,150)
# TSplot(6,300)
# TSplot(5,250)
# TSplot(7,250)
# TSplot(5,500)
# TSplot(6,550)
# TSplot(6,800)
# spruce_TS()
# xlim([34.5,35.2])
#
#
#
#
#
# figure()
# TSplot('6corr',150)
# TSplot('6corr',300)
# TSplot('5corr',250)
# TSplot(7,250)
# TSplot('5corr',500)
# TSplot(6,550)
# TSplot(6,800)
# spruce_TS()
# xlim([34.9,35.1])
#
#
#
# sal['5corr'].keys()
#
# figure()
# for key in sort(list(sal['5corr'].keys()))[-2:]:
#     TSplot(5,key)
#
# figure()
# for key in sort(list(sal['5corr'].keys()))[-2:]:
#     TSplot('5corr',key)
#
# sal['5corr'].keys()
#
#
#
#
# figure()
# TSplot('5corr',250)
# TSplot(6,300)
# TSplot(7,250)
# spruce_TS()
# xlim([34.5,35.2])
#
# sal['5corr'].keys()
#
# TSplot('4corr',400)
# TSplot('5corr',250)
# TSplot('5corr',500)
# TSplot(6,550)
# TSplot('5corr',1000)
# TSplot(6,1050)
# TSplot('5corr',1300)
# TSplot(6,1550)
# TSplot(6,1850)
# spruce_TS()
# xlim([34.5,35.2])
#
# figure()
# TSplot(6,550)
# TSplot(6,800)
#
#
# # TSplot(6,150)
# spruce_TS()
# xlim([33.5,35.5])
# savefig('../notes/figures/CF6_TS.png',bbox_inches='tight')
# figure()
# TSplot(6,100)
# TSplot('6corr',300)
# TSplot('6corr',150)
# spruce_TS()
# xlim([33.5,35.5])
# savefig('../notes/figures/CF6_TScorr.png',bbox_inches='tight')
# figure(figsize=(12,3))
# plot(date[6][100],(sal[6][300]-sal['6corr'][150])/0.1)
# plot(date[6][100],tmp[6][300]-tmp[6][150])
# axhline(0)
#
# figure(figsize=(12,3))
#
# axhline(0)
#
# figure(figsize=(12,3))
# plot(date[6][100],-sal['6corr'][150]+sal[6][300])
# axhline(0)
#
#
# figure(figsize=(12,3))
# plot(date[6][100],-tmp['6corr'][150]+tmp[6][300])
# axhline(0)
#
#
#
# figure(figsize=(12,3))
# plot(date[6][100],sal['6corr'][150]-sal['6corr'][100])
# axhline(0)
#
#
# figure(figsize=(12,3))
# plot(date[6][100],tmp['6corr'][150]-tmp['6corr'][100])
# axhline(0)
#
#
#
# figure(figsize=(12,3))
# plot(date[6][100],sal['6corr'][100],label='CF6, 100db')
# plot(date[6][100],sal['6corr'][150],label='CF6, 150db corr')
# plot(date[6][100],sal[6][300],label='CF6, 300db')
# legend(loc=(1.05,0.5),title='nominal pressure')
# xlabel('time (days)')
# ylabel('salinity')
# title('CF6 salinity correction')
# ylim([34.5,35.1])
# savefig('../notes/figures/CF6_salcorr.png',bbox_inches='tight')
#
# figure(figsize=(12,3))
# for key in sort(list(sal[6].keys()))[:3]:
#     plot(time[6][key],sal[6][key],label=key)
# plot(t100,scorr6_tot,label='150, corrected')
# plot(t100,line6_tot,color='orange')
# legend(loc=(1.05,0.5),title='nominal pressure')
# xlabel('time (days)')
# ylabel('salinity')
# title('CF6 salinity correction')
# xlim([23570,24350])
# ylim([34.5,35.1])
# savefig('../notes/figures/CF6_salcorr.png',bbox_inches='tight')


#####################################################################
#####################################################################
#####################################################################
###########                 CF7                    #################
#####################################################################
#####################################################################
#####################################################################

#This plot reassured me that it makes sense for the 250 instrument is fresher than the one above it in summer.
# for key in sort(list(salfinal[6].keys()))[:5]:
#     plot(date[6][key],salfinal[6][key],label=key)
# legend(loc=(1.05,0.5))

offstime=170+7.358e5

date[7][100][(time[7][100]>=offstime-1) &(time[7][100]<offstime)]


sal['7corr']={}
liney[7]={}
figure(figsize=(14,3))
plot(date[7][100],sal[7][100],label='Original record')
sal['7corr'][100]=sal[7][100].copy()
liney[7][100],sal['7corr'][100][time[7][100]<=offstime]=removeline(time[7][100][time[7][100]<=offstime],sal[7][100][time[7][100]<=offstime])
plot(date[7][100][time[7][100]<=offstime],liney[7][100],'b')
sal['7corr'][100][time[7][100]<=offstime]=sal['7corr'][100][time[7][100]<=offstime]-0.5*(liney[7][100][0]-liney[7][100][-1])
sal['7corr'][100][time[7][100]>offstime]=sal[7][100].copy()[time[7][100]>offstime]+0.5*(liney[7][100][0]-liney[7][100][-1])
plot(date[7][100],sal['7corr'][100],label='First Detrend')
liney[7][100],sal['7corr'][100]=removeline(time[7][100],sal['7corr'][100])
plot(date[7][100],liney[7][100],color='orange')
plot(date[7][100],sal['7corr'][100],label='Second Detrend')
sal['7corr'][100]=sal['7corr'][100]-0.0085
plot(date[7][100],sal['7corr'][100],label='Final record')
title('Salinity time series at 100db instrument on CF7')
ylabel('salinity')
legend()
# grid('on')
savefig('../figures/salcalib_15min/CF7_100_scorr'+savename+'.png',bbox_inches='tight')


# offset 50db instrument
sal['7corr'][50]=sal[7][50]-0.025


#####################################################################
#####################################################################
#####################################################################
#####################################################################
##########              Final analysis               ################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

#####################################################################
# Create corrected sal and tmp dictionaries to reflect changes
#####################################################################

# form basis
salfinal={}
for key in range(1,8):
    salfinal[key]={}
    for key2 in sal[key]:
        salfinal[key][key2]=sal[key][key2]


# replace relevant cases with corrected values
for key in range(1,8):
    for key2 in sal[str(key)+'corr']:
        salfinal[key][key2]=sal[str(key)+'corr'][key2]


#####################################################################
######### Calculate and plot density before and after ###############
#####################################################################


pden_precorr={}
pden={}
for mm in range(1,8):
    pden_precorr[mm]={}
    pden[mm]={}
    for key in sort(list(sal[mm].keys())):
        pden_precorr[mm][key]=gsw.sigma0(sal[mm][key],gsw.CT_from_pt(sal[mm][key],tmp[mm][key]))
        pden[mm][key]=gsw.sigma0(salfinal[mm][key],gsw.CT_from_pt(salfinal[mm][key],tmp[mm][key]))



for mm in range(2,8):
    prslist=list(sort(list(salfinal[mm].keys())))
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
        plot(date[mm][kk][:lesslen],pden[mm][prslist[ii+1]][:lesslen]-pden[mm][kk][:lesslen])
        title('CF'+str(mm)+': '+str(prslist[ii+1])+'db-'+str(kk)+'db; from '+str(sum((pden_precorr[mm][prslist[ii+1]][:lesslen]-pden_precorr[mm][kk][:lesslen])<0))+' to '+str(sum((pden[mm][prslist[ii+1]][:lesslen]-pden[mm][kk][:lesslen])<0))+' inversions')
        ylabel('Density difference (kg m$^{-3}$)')
        savefig('../figures/salcalib_15min/perinst/CF'+str(mm)+'_'+str(prslist[ii+1])+'db-'+str(kk)+'db'+savename+'.png',bbox_inches='tight')

#####################################################################
######### Plot salinity time series before and after ################
#####################################################################

for mm in range(1,8):
    figure(figsize=(14,3))
    for key in sort(list(sal[mm].keys())):
        plot(date[mm][key],sal[mm][key],label=str(key)+' db')
    title('Salinity at CF'+str(mm)+' BEFORE corrections')
    ylabel('Salinity')
    legend()
    savefig('../figures/salcalib_15min/CF'+str(mm)+'_salbefore'+savename+'.png',bbox_inches='tight')

for mm in range(1,8):
    figure(figsize=(14,3))
    for key in sort(list(sal[mm].keys())):
        plot(date[mm][key],salfinal[mm][key],label=str(key)+' db')
    title('Salinity at CF'+str(mm)+' AFTER corrections')
    ylabel('Salinity')
    legend()
    savefig('../figures/salcalib_15min/CF'+str(mm)+'_salafter'+savename+'.png',bbox_inches='tight')

#####################################################################
######### Plot TS space before and after ############################
#####################################################################

for mm in range(1,8):
    figure()
    for key in sort(list(sal[mm].keys())):
        TSplotafter(mm,key)
    spruce_TS()
    title('CF'+str(mm)+': $\Theta$-S space AFTER corrections')
    if mm>=6:
        xlim([34,35.2])
    savefig('../figures/salcalib_15min/CF'+str(mm)+'_TSafter'+savename+'.png',bbox_inches='tight')

for mm in range(1,8):
    figure()
    for key in sort(list(sal[mm].keys())):
        TSplot(mm,key)
    spruce_TS()
    title('CF'+str(mm)+': $\Theta$-S space BEFORE corrections')
    if mm>=6:
        xlim([34,35.2])
    savefig('../figures/salcalib_15min/CF'+str(mm)+'_TSbefore'+savename+'.png',bbox_inches='tight')


#####################################################################
######### Plot TS space correspondence between moorings #############
#####################################################################

dlim=[0,100,300,750,2000]

for ii,dd in enumerate(dlim[:-1]):

    figure()
    for mm in range(4,8):
        for key in sort(list(sal[mm].keys())):
            if (key>dd) & (key<=dlim[ii+1]):
                TSplotafter(mm,key)
    spruce_TS()
    # legend(loc=(1.05,0.2))
    title('CF4-7: '+str(dd)+'db-'+str(dlim[ii+1])+'db $\Theta$-S space AFTER corrections')
    if dd>=100:
        xlim([33,35.5])
    if dd>=300:
        xlim([34,35.2])
    if dd>=500:
        xlim([34.8,35.1])
    savefig('../figures/salcalib_15min/CF4-7_TSafter_'+str(dd)+'-'+str(dlim[ii+1])+savename+'.png',bbox_inches='tight')

    figure()
    for mm in range(4,8):
        for key in sort(list(sal[mm].keys())):
            if (key>dd) & (key<=dlim[ii+1]):
                TSplot(mm,key)
    spruce_TS()
    # legend(loc=(1.05,0.2))
    title('CF4-7: '+str(dd)+'db-'+str(dlim[ii+1])+'db $\Theta$-S space BEFORE corrections')
    if dd>=100:
        xlim([33,35.5])
    if dd>=300:
        xlim([34,35.2])
    if dd>=500:
        xlim([34.8,35.1])
    savefig('../figures/salcalib_15min/CF4-7_TSbefore'+str(dd)+'-'+str(dlim[ii+1])+savename+'.png',bbox_inches='tight')


shelflim=200

figure()
for mm in range(1,5):
    for key in sort(list(sal[mm].keys())):
        if (key<=shelflim):
            TSplotafter(mm,key)
spruce_TS()
# legend(loc=(1.05,0.2))
title('CF1-4 $\Theta$-S space AFTER corrections')
savefig('../figures/salcalib_15min/CF1-4_TSafter'+savename+'.png',bbox_inches='tight')

figure()
for mm in range(1,5):
    for key in sort(list(sal[mm].keys())):
        if (key<=shelflim):
            TSplot(mm,key)
spruce_TS()
# legend(loc=(1.05,0.2))
title('CF1-4 $\Theta$-S space BEFORE corrections')
savefig('../figures/salcalib_15min/CF1-4_TSbefore'+savename+'.png',bbox_inches='tight')

for key in date:
    month[key]={}
    for key2 in date[key]:
        month[key][key2]=[dd.month for dd in date[key][key2]]

pickle.dump([date,month,prs,salfinal,tmp],
            open('../pickles/TSdailydic/TS_daily_dic_wcorr'+savename+'.pickle','wb'))


max_sdiff={}
mean_sdiff={}
max_pdendiff={}
mean_pdendiff={}
for key in sal:
    if 'corr' in str(key):
        max_sdiff[key]={}
        mean_sdiff[key]={}
        max_pdendiff[key]={}
        mean_pdendiff[key]={}
        for dpkey in sal[key]:
            max_sdiff[key][dpkey]=round(nanmax(abs(sal[key][dpkey]-sal[int(key[0])][dpkey])),3)
            mean_sdiff[key][dpkey]=round(nanmean(sal[key][dpkey]-sal[int(key[0])][dpkey]),3)
            max_pdendiff[key][dpkey]=round(nanmax(abs(pden[int(key[0])][dpkey]-pden_precorr[int(key[0])][dpkey])),3)
            mean_pdendiff[key][dpkey]=round(nanmean(pden[int(key[0])][dpkey]-pden_precorr[int(key[0])][dpkey]),3)


max_pdendiff
