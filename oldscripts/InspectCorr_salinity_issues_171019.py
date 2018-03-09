from aux_funcs import *


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

        if moornum==1:
            [cf1date,cf1time,cf1prs,cf1mnprs,cf1sal,cf1tmp]=pd.read_pickle(open('../pickles/TSinterp/CF1_recon.pickle',
                                                                             'rb'))


        prskey=int(round(nanmean(prs_hrly)/50.0)*50.0)
        timekey=mean(time_hrly)

        sal[moornum][prskey]=sal_hrly
        tmp[moornum][prskey]=tmp_hrly
        prs[moornum][prskey]=prs_hrly
        date[moornum][prskey]=date_hrly
        time[moornum][prskey]=time_hrly


# Keeping reconstruction out of this part -- do recon a couple different ways after I've calibrated.
# sal['1recon']={}
# tmp['1recon']={}
# date['1recon']={}
# prs['1recon']={}
# sal['1recon'][50]=cf1sal['0']
# sal['1recon'][100]=cf1sal['3']
# tmp['1recon'][50]=cf1tmp['0']
# tmp['1recon'][100]=cf1tmp['3']
# date['1recon'][50]=cf1date['0']
# date['1recon'][100]=cf1date['3']
# prs['1recon'][50]=cf1prs['0']
# prs['1recon'][100]=cf1prs['3']


month={}
for key in date:
    month[key]={}
    for key2 in date[key]:
        month[key][key2]=[dd.month for dd in date[key][key2]]

###################################################################
################## Some useful functions #########################
###################################################################


def removeline(tvec,svec,anchor=0):
    linefunc=poly1d(polyfit(tvec,svec, 1))
    liney=linefunc(tvec)
    corrected=svec-liney+liney[anchor]
    return liney,corrected

# set up basic ts plotting framework

def TSplot(moor,prs,alph=0.3):
    plot(sal[moor][prs],tmp[moor][prs],'o',
         markersize=5,label='CF'+str(moor)+', '+str(prs)+' db',
         alpha=alph,mew=0)

def spruce_TS():
    dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k')
    ylabel('temperature')
    xlabel('salinity')
    legend(loc=(1.05,0.2),numpoints=1,markerscale=2)


def TSplot_seasonal(moor,prs):
    scatter(sal[moor][prs],tmp[moor][prs],c=month[moor][prs],
            cmap=cm.brg_r,lw=0,alpha=0.5,vmin=1,vmax=12)

def spruce_seasonal():
    cbar=colorbar(ticks=range(1,13,2))
    cbar.set_ticklabels(['January','March','May','July','September','November'])
    ylabel('temperature')
    xlabel('salinity')
    dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k')
    clabel(dens)


def plotnfilt4(what):
    figure(figsize=(14,3))
    divit=30
    plot(date[4][350],what)
    plot(date[4][350][::divit],hrly_ave(what,divit),'ro-')
    plot(date[4][350][::divit*6],hrly_ave(what,divit*6),'yo-')


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
# savefig('../figures/salcorr/CF4_salatend_TS.pdf',bbox_inches='tight')

#####################################################################
#####################################################################
#####################################################################
###########                 CF1                     #################
#####################################################################
#####################################################################
#####################################################################

liney[1]={}
sal['1corr']={}
liney[1][50],sal['1corr'][50]=removeline(time[1][50],sal[1][50])

sal['1corr'][100]=sal[1][100].copy()

sal['1corr'][100][-2:]=NaN

sal['1corr'][100]

sal['1corr'][150]=sal[1][150].copy()

figure(figsize=(14,3))
for key in sal[1]:
    plot(date[1][key],sal['1corr'][key],label=key)
legend(loc=(1.05,0.2))

figure(figsize=(14,3))
for key in sal[1]:
    plot(date[1][key],tmp[1][key],label=key)
legend(loc=(1.05,0.2))

figure(figsize=(14,3))
for key in sal[1]:
    plot(date[1][key],sal[1][key],label=key)
legend(loc=(1.05,0.2))


pden={}
mm=1
pden[mm]={}
for key in sort(list(sal[mm].keys())):
    pden[mm][key]=gsw.sigma0(sal[1][key],gsw.CT_from_t(sal[1][key],tmp[mm][key],prs[mm][key]))

(pden[1][150][:len100]-pden[1][100])[-2]

plot(pden[1][150][:len100]-pden[1][100])

len100=len(sal[mm][100])
len50=len(sal[mm][50])

figure(figsize=(12,3))
plot(date[mm][100],pden[mm][150][:len100]-pden[mm][100])
axhline(0)
title('CF'+str(mm)+': 150db-100db; '+str(sum((pden[mm][150][:len100]-pden[mm][100])<0))+' inv profiles')
savefig('../figures/invercheck/perinst/CF'+str(mm)+'_200db-100db'+testname+'.png')

figure(figsize=(12,3))
plot(date[mm][50],pden[mm][100][:len50]-pden[mm][50])
axhline(0)
title('CF'+str(mm)+': 100db-50db; '+str(sum((pden[mm][100][:len50]-pden[mm][50])<0))+' inv profiles')
savefig('../figures/invercheck/perinst/CF'+str(mm)+'_100db-50db'+testname+'.png')


figure(figsize=(12,3))
plot(date[mm][50],pden[mm][150][:len50]-pden[mm][50])
axhline(0)
title('CF'+str(mm)+': 150db-50db; '+str(sum((pden[mm][150][:len50]-pden[mm][50])<0))+' inv profiles')
savefig('../figures/invercheck/perinst/CF'+str(mm)+'_150db-50db_'+testname+'.png')

#####################################################################
#####################################################################
#####################################################################
###########                 CF2                     #################
#####################################################################
#####################################################################
#####################################################################


salfinal={}
tmpfinal={}
datefinal={}
for key in range(2,3):
    salfinal[key]={}
    tmpfinal[key]={}
    for key2 in sal[key]:
        salfinal[key][key2]=sal[key][key2]
        tmpfinal[key][key2]=tmp[key][key2]

liney[2]={}
sal['2corr']={}
for key in sal[2]:
    liney[2][key],sal['2corr'][key]=removeline(time[2][key],sal[2][key])


sal['2corr'][50]=sal['2corr'][50]-0.065
sal['2corr'][100]=sal['2corr'][100]
sal['2corr'][200]=sal['2corr'][200]+0.01

tmp['2corr']={}

for key in sal[str(mm)+'corr']:
        salfinal[mm][key]=sal[str(mm)+'corr'][key]
        tmp['2corr'][key]=tmp[2][key]

# for key in sal[2]:
#     TSplot(2,key)
#
# for key in sal[2]:
#     TSplot('2corr',key)
#
# figure(figsize=(14,3))
# for key in sal[2]:
#     plot(date[2][key],sal['2corr'][key],label=key)
# legend(loc=(1.05,0.2))
#
# figure(figsize=(14,3))
# for key in sal[2]:
#     plot(date[2][key],sal[2][key],label=key)
# legend(loc=(1.05,0.2))


testname='orig'

pden={}
mm=2
pden[mm]={}
for key in sort(list(salfinal[mm].keys())):
    pden[mm][key]=gsw.sigma0(sal[mm][key],gsw.CT_from_t(sal[mm][key],tmp[mm][key],prs[mm][key]))


figure(figsize=(12,3))
plot(date[mm][kk],pden[mm][200]-pden[mm][100])
axhline(0)
title('CF'+str(mm)+': 200db-100db; '+str(sum((pden[mm][200]-pden[mm][100])<0))+' inv profiles')
savefig('../figures/invercheck/perinst/CF'+str(mm)+'_200db-100db'+testname+'.png')

figure(figsize=(12,3))
plot(date[mm][kk],pden[mm][100]-pden[mm][50])
axhline(0)
title('CF'+str(mm)+': 100db-50db; '+str(sum((pden[mm][100]-pden[mm][50])<0))+' inv profiles')
savefig('../figures/invercheck/perinst/CF'+str(mm)+'_100db-50db'+testname+'.png')


figure(figsize=(12,3))
plot(date[mm][kk],pden[mm][200]-pden[mm][50])
axhline(0)
title('CF'+str(mm)+': 200db-50db; '+str(sum((pden[mm][200]-pden[mm][50])<0))+' inv profiles')
savefig('../figures/invercheck/perinst/CF'+str(mm)+'_200db-50db_'+testname+'.png')


#####################################################################
#####################################################################
#####################################################################
###########                 CF3                     #################
#####################################################################
#####################################################################
#####################################################################

# Going to leave out 100db instrument's salinity as it is buggy and not long -- not worth keeping around.

sal['3corr'][100]=sal[3][100].copy()
sal['3corr'][100][:]=NaN

figure(figsize=(14,3))
for key in sal[3]:
    plot(date[3][key],tmp[3][key],label=key)
legend(loc=(1.05,0.2))

figure(figsize=(14,3))
for key in sal[3]:
    plot(date[3][key],sal[3][key],label=key)
legend(loc=(1.05,0.2))


figure(figsize=(14,3))
for key in sal[3]:
    plot(date[3][key],sal[3][key],label=key)
legend(loc=(1.05,0.2))
xlim([datetime.datetime(2014,8,10),datetime.datetime(2014,12,1)])




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

klist=[50,100,200]

figure(figsize=(10,3))
for key in klist:
    plot(date[4][key],sal[4][key],'-',label=key)
legend(loc=(1.05,0.5))
grid('on')


saldif4=sal[4][100]-sal[4][50]

driftime=24020
valitime=24285


#### Drift is in several phases

# Apply correction to first part
salcorr4_first=saldif4[time[4][100]<=driftime]
tcorr4_first=time[4][100][time[4][100]<=driftime]
liney4_first,corr4_first=removeline(tcorr4_first,salcorr4_first)

# Then second part
timeind4=(time[4][100]>driftime) & (time[4][100]<valitime)
salcorr4=saldif4[timeind4]
tcorr4=time[4][100][timeind4]
liney4,corr4=removeline(tcorr4,salcorr4)
newdiff=saldif4.copy()
newdiff[timeind4]=corr4
newdiff[time[4][100]>valitime]=NaN

figure(figsize=(12,3))
plot(time[4][50],saldif4)
axvline(driftime,color='k')
axvline(valitime,color='k')
plot(tcorr4_first,liney4_first,'k')
plot(tcorr4,liney4,'k')
grid('on')
title('CF4 100db salinity correction')
plot(time[4][50],newdiff)
savefig('../notes/figures/CF4_sdiff.pdf',bbox_inches='tight')



sal['4corr']={}
sal['4corr'][100]=sal[4][100].copy()
sal['4corr'][100][time[4][100]>valitime]=NaN
sal['4corr'][100][time[4][100]<=driftime]=corr4_first+sal[4][50][time[4][50]<=driftime]+0.07
sal['4corr'][100][timeind4]=corr4+sal[4][50][timeind4]+liney4_first[-1]
sal['4corr'][100]=sal['4corr'][100]-0.05

liney[4]={}
liney[4][50],sal['4corr'][50]=removeline(time[4][50],sal[4][50])


tmp['4corr'][50]=tmp[4][50]
tmp['4corr'][100]=tmp[4][100]


#####################################################################
###### Another drift at CF4 @ 350db #########################
#####################################################################

#Remove linear fit to salinity at 350db


liney[4][350],sal['4corr'][350]=removeline(time[4][350],sal[4][350])
sal['4corr'][350][date[4][350]>datetime.datetime(2016,1,1)]=sal['4corr'][350][date[4][350]>datetime.datetime(2016,1,1)]+0.025



#Remove linear fit to difference between 400 and 350db instruments
dsplit400=datetime.datetime(2015,10,1)
dateind400=[date[4][350]<dsplit400]
saldiff434=(sal[4][400]-sal['4corr'][350])[dateind400]
liney434,saldiff434_corr=removeline(time[4][350][dateind400],saldiff434)
newdiff434=(sal[4][400]-sal['4corr'][350]).copy()
newdiff434[dateind400]=saldiff434_corr
newdiff434[date[4][350]>=dsplit400]=(sal[4][400]-sal['4corr'][350])[date[4][350]>=dsplit400]-0.9*liney434[-1]

sal['4corr'][400]=sal['4corr'][350]+newdiff434+0.022


# Extraneous plotting things that I can refer back to
# figure(figsize=(14,3))
# plot(time[4][50],sal[4][50])
# plot(time[4][50],liney[4][50])
# plot(time[4][50],sal['4corr'][50])
#
#
#
# TSplot(4,100)
# TSplot(4,50)
#
#
# TSplot('4corr',100)
# TSplot('4corr',50)
#
# plotnfilt4(newdiff434)
# plot(date[4][350],sal[4][400]-sal['4corr'][350])
# plot(date[4][350][dateind400],liney434)
# plot(date[4][350][dateind400],saldiff434_corr)
# axhline(0)
# savefig('../figures/invercheck/perinst/CF4_saldiff_400-350_corroct.png')
# min(newdiff434)
# figure(figsize=(14,3))
# for key in [50,100,200]:
#     plot(time[4][key],sal[4][key],label=key)
# plot(time[4][100],sal['4corr'][100],label='100, corrected')
# axvline(driftime,color='k')
# axvline(valitime,color='k')
# legend(loc=(1.05,0.5),title='nominal pressure')
# xlabel('time (days)')
# ylabel('salinity')
# title('CF4 100db salinity correction ')
# savefig('../notes/figures/CF4_salcorr.pdf',bbox_inches='tight')
#
 # figure()
# TSplot_seasonal(4,50)
# TSplot_seasonal(4,100)
# TSplot_seasonal(4,200)
# spruce_seasonal()
# title('Before correction')
# savefig('../notes/figures/CF4_seasonalTS.png',bbox_inches='tight')
#
# figure()
# TSplot_seasonal(4,50)
# TSplot_seasonal('4corr',100)
# TSplot_seasonal(4,200)
# spruce_seasonal()
# title('After correction')
# savefig('../notes/figures/CF4_seasonalTScorr.png',bbox_inches='tight')
#
#
#
# figure()
# TSplot_seasonal(4,100)
# spruce_seasonal()
#
#
# figure()
# TSplot_seasonal('4corr',100)
# spruce_seasonal()
#
#
# figure()
# TSplot(4,50)
# TSplot(4,200)
# TSplot(4,100)
# spruce_TS()
# xlim([32.5,35.5])
# # ylim([0,4])
# savefig('../notes/figures/CF4_TS.png',bbox_inches='tight')
# figure()
# TSplot(4,50)
# TSplot(4,200)
# TSplot('4corr',100)
# spruce_TS()
# xlim([32.5,35.5])
# # ylim([0,4])
# savefig('../notes/figures/CF4_TScorr.png',bbox_inches='tight')
# figure(figsize=(14,3))
# plot(time[4][50],saldif4)
# axvline(driftime,color='r')
# plot(tcorr4,corr4)
# plot(tcorr4,liney4)
# plot(tcorr4_first,liney4_first)
# axvline(valitime,color='r')
# grid('on')
#
# figure(figsize=(14,3))
# plot(time[4][100],sal[4][100])
# axvline(driftime,color='r')
# axvline(valitime,color='r')
# plot(tcorr4,corr4+sal[4][50][timeind4])
# plot(tcorr4_first,corr4_first+sal[4][50][time[4][50]<=driftime])

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

date5=datetime.datetime(1950,1,1)+datetime.timedelta(days=int(24305))


figure(figsize=(14,3))
for key in sal[5]:
    plot(date[5][key],sal[5][key],label=key)
legend(loc=(1.05,0.25),title='nominal pressure')
axvline(date5,color='k')
xlabel('time (days)')
ylabel('salinity')
title('CF5 salinity correction')
savefig('../notes/figures/CF5_salcorr.pdf',bbox_inches='tight')


sal['5corr']={}
tmp['5corr']={}
sal['5corr'][250]=sal[5][250].copy()
tmp['5corr'][250]=tmp[5][250].copy()
sal['5corr'][250][date[5][250]>date5]=nan

sal[5].keys()

figure()
TSplot(5,100)
TSplot(5,250)
TSplot(5,500)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF5_TScomp.png',bbox_inches='tight')

figure()
TSplot(5,250,alph=1)
TSplot('5corr',250,alph=1)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF5_TScorr.png',bbox_inches='tight')


#####################################################################
### Density inversions/drift on CF5, 500-1300db
#####################################################################

keylist5=[500,750,1000,1300]

TSplot(5,500)
TSplot(5,750)
TSplot(5,1000)
TSplot(5,1300)
spruce_TS()
xlim([34.8,35.1])

figure(figsize=(14,3))
for key in keylist5:
    plot(date[5][key],sal[5][key],label=key)
legend()

# Remove linear trend in all records

sal[5].keys()


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



mm=6
sal[str(mm)+'corr']={}
klist=[150]
for key in klist:
    liney[mm][key],sal[str(mm)+'corr'][key]=removeline(time[mm][key],sal[mm][key])

sal['6corr'][150]=sal['6corr'][150]-0.0122
sal['6corr'][100]=sal[6][100]-0.024


## Some valuable thoughts on TS space below, keeping around for now.


tmp[str(key)+'5corr']={}
for key2 in sal['5corr']:
        tmp['5corr'][key2]=tmp[5][key2]



sal[5].keys()

sal[7].keys()

figure()
# TSplot(6,300)
# TSplot(7,250)
TSplot('6corr',150)
TSplot('6corr',100)
# TSplot('5corr',100)
TSplot(7,50)
spruce_TS()
xlim([34.4,35.2])




figure()
TSplot('6corr',150)
TSplot('6corr',100)
# TSplot('5corr',100)
TSplot(7,100)

figure()
TSplot(6,100)
TSplot('6corr',150)
# TSplot('5corr',100)
TSplot(7,100)

figure()
TSplot(6,300)
TSplot(6,550)

sal[7].keys()

figure()
TSplot(6,150)
TSplot(6,300)
TSplot(5,250)
TSplot(7,250)
TSplot(5,500)
TSplot(6,550)
TSplot(6,800)
spruce_TS()
xlim([34.5,35.2])





figure()
TSplot('6corr',150)
TSplot('6corr',300)
TSplot('5corr',250)
TSplot(7,250)
TSplot('5corr',500)
TSplot(6,550)
TSplot(6,800)
spruce_TS()
xlim([34.9,35.1])



sal['5corr'].keys()

figure()
for key in sort(list(sal['5corr'].keys()))[-2:]:
    TSplot(5,key)

figure()
for key in sort(list(sal['5corr'].keys()))[-2:]:
    TSplot('5corr',key)

sal['5corr'].keys()




figure()
TSplot('5corr',250)
TSplot(6,300)
TSplot(7,250)
spruce_TS()
xlim([34.5,35.2])

sal['5corr'].keys()

TSplot('4corr',400)
TSplot('5corr',250)
TSplot('5corr',500)
TSplot(6,550)
TSplot('5corr',1000)
TSplot(6,1050)
TSplot('5corr',1300)
TSplot(6,1550)
TSplot(6,1850)
spruce_TS()
xlim([34.5,35.2])

figure()
TSplot(6,550)
TSplot(6,800)


# TSplot(6,150)
spruce_TS()
xlim([33.5,35.5])
savefig('../notes/figures/CF6_TS.png',bbox_inches='tight')
figure()
TSplot(6,100)
TSplot('6corr',300)
TSplot('6corr',150)
spruce_TS()
xlim([33.5,35.5])
savefig('../notes/figures/CF6_TScorr.png',bbox_inches='tight')
figure(figsize=(12,3))
plot(date[6][100],(sal[6][300]-sal['6corr'][150])/0.1)
plot(date[6][100],tmp[6][300]-tmp[6][150])
axhline(0)

figure(figsize=(12,3))

axhline(0)

figure(figsize=(12,3))
plot(date[6][100],-sal['6corr'][150]+sal[6][300])
axhline(0)


figure(figsize=(12,3))
plot(date[6][100],-tmp['6corr'][150]+tmp[6][300])
axhline(0)



figure(figsize=(12,3))
plot(date[6][100],sal['6corr'][150]-sal['6corr'][100])
axhline(0)


figure(figsize=(12,3))
plot(date[6][100],tmp['6corr'][150]-tmp['6corr'][100])
axhline(0)



figure(figsize=(12,3))
plot(date[6][100],sal['6corr'][100],label='CF6, 100db')
plot(date[6][100],sal['6corr'][150],label='CF6, 150db corr')
plot(date[6][100],sal[6][300],label='CF6, 300db')
legend(loc=(1.05,0.5),title='nominal pressure')
xlabel('time (days)')
ylabel('salinity')
title('CF6 salinity correction')
ylim([34.5,35.1])
savefig('../notes/figures/CF6_salcorr.pdf',bbox_inches='tight')


# In[57]:



figure(figsize=(12,3))
for key in sort(list(sal[6].keys()))[:3]:
    plot(time[6][key],sal[6][key],label=key)
plot(t100,scorr6_tot,label='150, corrected')
plot(t100,line6_tot,color='orange')
legend(loc=(1.05,0.5),title='nominal pressure')
xlabel('time (days)')
ylabel('salinity')
title('CF6 salinity correction')
xlim([23570,24350])
ylim([34.5,35.1])
savefig('../notes/figures/CF6_salcorr.pdf',bbox_inches='tight')


#####################################################################
#####################################################################
#####################################################################
###########                 CF7                    #################
#####################################################################
#####################################################################
#####################################################################

#####################################################################
###################### More drift at CF7! #########################
#####################################################################




#This plot reassured me that 250 instrument is fresher than above in summer.
for key in sort(list(salfinal[6].keys()))[:5]:
    plot(date[6][key],salfinal[6][key],label=key)
# axvline(offstime,color='red')
legend(loc=(1.05,0.5))


offstime=170+7.358e5

date[7][100][(time[7][100]>=offstime-1) &(time[7][100]<offstime)]


for key in sort(list(sal[7].keys()))[:5]:
    plot(time[7][key],sal[7][key],label=key)
axvline(offstime,color='red')
xlim([735900,736100])
legend(loc=(1.05,0.5))



sal['7corr']={}
liney[7]={}
sal['7corr'][100]=sal[7][100].copy()
liney[7][100],sal['7corr'][100][time[7][100]<=offstime]=removeline(time[7][100][time[7][100]<=offstime],sal[7][100][time[7][100]<=offstime])
sal['7corr'][100][time[7][100]<=offstime]=sal['7corr'][100][time[7][100]<=offstime]-0.5*(liney[7][100][0]-liney[7][100][-1])
sal['7corr'][100][time[7][100]>offstime]=sal[7][100].copy()[time[7][100]>offstime]+0.5*(liney[7][100][0]-liney[7][100][-1])

liney[7][100],sal['7corr'][100]=removeline(time[7][100],sal['7corr'][100])

sal['7corr'][100]=sal['7corr'][100]-0.0085
sal['7corr'][50]=sal[7][50]-0.025

figure(figsize=(14,3))
plot(date[7][50],sal[7][50])
plot(date[7][100],sal[7][100])
# plot(date[7][100][time[7][100]<=offstime],liney[7][100])
plot(date[7][100],sal['7corr'][100])
plot(date[6][100],sal['6corr'][100])
ylim([34.6,35.2])
# figure(figsize=(14,3))
# plot(date[7][50],tmp[7][50])
# plot(date[7][100],tmp[7][100])
# xlim([datetime.datetime(2014,8,10),datetime.datetime(2015,6,1)])
#
#
# figure(figsize=(14,3))
# plot(date[7][50],sal[7][50])
# plot(date[7][100],sal['7corr'][100])
# ylim([34.6,35.2])
# xlim([datetime.datetime(2014,8,10),datetime.datetime(2015,6,1)])




figure()
TSplot(6,100)
TSplot(7,100)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF7_TS_CF6comp.png',bbox_inches='tight')

figure()
TSplot(6,100)
# TSplot(6,300)
# TSplot('6corr_tot',150)
# TSplot(7,50)
# TSplot(7,250)
TSplot('7corr',100)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF7_TScorr_CF6comp.png',bbox_inches='tight')



figure(figsize=(14,3))
for key in sort(list(sal[7].keys()))[:5]:
    plot(time[7][key],sal[7][key],label=key)
plot(time[7][100],sal['7corr'][sort(list(sal['7corr'].keys()))[-4]],label='100, corrected')
axvline(offstime,color='k')
legend(loc=(1.05,0.25),title='nominal pressure')
xlabel('time (days)')
ylabel('salinity')
title('CF7 salinity correction')
savefig('../notes/figures/CF7_salcorr.pdf',bbox_inches='tight')


figure()
TSplot(6,100)
TSplot('6corr',150)
# TSplot('6corr_tot',150)
TSplot(6,300)
spruce_TS()


figure()
TSplot(7,50)
TSplot(7,250)
TSplot(7,100)
spruce_TS()
xlim([34.5,35.2])
savefig('../notes/figures/CF7_TS.png',bbox_inches='tight')


figure()
TSplot(7,50)
TSplot(7,250)
TSplot('7corr',sort(list(sal['7corr'].keys()))[-4])
spruce_TS()
xlim([34.5,35.2])
savefig('../notes/figures/CF7_TScorr.png',bbox_inches='tight')



#####################################################################
# Create corrected sal and tmp dictionaries to reflect changes
#####################################################################


salfinal={}
tmpfinal={}
datefinal={}
for key in range(1,8):
    salfinal[key]={}
    tmpfinal[key]={}
    for key2 in sal[key]:
        salfinal[key][key2]=sal[key][key2]
        tmpfinal[key][key2]=tmp[key][key2]


salfinal[1][50]=sal['1recon'][50]
salfinal[1][100]=sal['1recon'][100]
tmpfinal[1][50]=tmp['1recon'][50]
tmpfinal[1][100]=tmp['1recon'][100]


date[1][50]=date['1recon'][50]
date[1][100]=date['1recon'][100]
prs[1][50]=prs['1recon'][50]
prs[1][100]=prs['1recon'][100]


salfinal[4][100]=sal['4corr'][100]
tmpfinal[4][100]=tmp['4corr'][100]

salfinal[4][350]=sal['4corr'][350]
salfinal[4][400]=sal['4corr'][400]


salfinal[5][250]=sal['5corr'][250]
tmpfinal[5][250]=tmp['5corr'][250]


salfinal[6][150]=sal['6corr_tot'][150]
tmpfinal[6][150]=tmp['6corr_tot'][150]


salfinal[7][100]=sal['7corr'][sort(list(sal['7corr'].keys()))[-4]]
tmpfinal[7][100]=tmp['7corr'][sort(list(sal['7corr'].keys()))[-4]]

del date['1recon'],date['4corr']
del month['1recon'],month['4corr']


#####################################################################
######### Check out density to see if there are any clues ##########
#####################################################################


# pden_precorr={}
# for mm in range(2,8):
#     pden_precorr[mm]={}
#     figure(figsize=(12,3))
#     for key in sort(list(sal[mm].keys())):
#         pden_precorr[mm][key]=gsw.sigma0(sal[mm][key],gsw.CT_from_t(sal[mm][key],tmp[mm][key],prs[mm][key]))
#         plot(date[mm][key],pden_precorr[mm][key],label=key)
#     title('CF'+str(mm)+' potential density (pre corrections)')
#     legend(loc=(1.05,0.1))
#     savefig('../figures/invercheck/CF'+str(mm)+'_pden_meas_precorr.png',bbox_inches='tight')



pden={}
for mm in range(1,8):
    pden[mm]={}
    figure(figsize=(12,3))
    for key in sort(list(salfinal[mm].keys())):
        pden[mm][key]=gsw.sigma0(salfinal[mm][key],gsw.CT_from_t(salfinal[mm][key],tmpfinal[mm][key],prs[mm][key]))
        plot(date[mm][key],pden[mm][key],label=key)
    title('CF'+str(mm)+' potential density')
    legend(loc=(1.05,0.1))
    savefig('../figures/invercheck/CF'+str(mm)+'_pden_meas.png',bbox_inches='tight')

for mm in range(1,8):
    prslist=sort(list(salfinal[mm].keys()))
    for ii,kk in enumerate(prslist[:-1]):
        figure(figsize=(12,3))
        plot(date[mm][kk],pden[mm][prslist[ii+1]]-pden[mm][kk])
        axhline(0)
        title('CF'+str(mm)+': '+str(prslist[ii+1])+'db-'+str(kk)+'db')
        savefig('../figures/invercheck/perinst/CF'+str(mm)+'_'+str(prslist[ii+1])+'db-'+str(kk)+'db_test.png')



# for mm in range():
#     prslist=sort(list(sal[mm].keys()))
#     for ii,kk in enumerate(prslist[:-1]):
#         figure(figsize=(12,3))
#         plot(date[mm][kk],pden_precorr[mm][prslist[ii+1]]-pden[mm][kk])
#         axhline(0)
#         title('CF'+str(mm)+': '+str(prslist[ii+1])+'db-'+str(kk)+'db')
#         savefig('../figures/invercheck/perinst/CF'+str(mm)+'_'+str(prslist[ii+1])+'db-'+str(kk)+'db_precorr.png')


len(date[1][50])
len(month[1][50])


pickle.dump([date,month,prs,salfinal,tmpfinal],
            open('../pickles/TSinterp/TS_daily_dic_wcorr.pickle','wb'))
