########################################################################################
########################################################################################
# This version is without daily averaging, and focusing on CF4-7
# for the express purpose of creating new fields for Feili +
########################################################################################
########################################################################################

from aux_funcs import *

savename='_nofilter'

############################################
# Load all T and S into large dictionaries
############################################

sal={}
tmp={}
date={}
prs={}
time={}
for moornum in range(4,8):
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



month={}
for key in date:
    month[key]={}
    for key2 in date[key]:
        month[key][key2]=[dd.month for dd in date[key][key2]]


# set up basic ts plotting framework

def TSplot(moor,prs,alph=0.3):
    plot(sal[moor][prs],tmp[moor][prs],'o',
         markersize=5,label='CF'+str(moor)+', '+str(prs)+' db',
         alpha=alph,mew=0)

def spruce_TS():
    dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k')
#     clabel(dens)
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
    savefig('../notes/figures/CF5_salcorr'+savename+'.pdf',bbox_inches='tight')


sal['5corr']={}
tmp['5corr']={}
sal['5corr'][250]=sal[5][250].copy()
tmp['5corr'][250]=tmp[5][250].copy()
sal['5corr'][250][date[5][250]>date5]=nan

sal[5].keys()


TSplot(5,100)
TSplot(5,250)
TSplot(5,500)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF5_TScomp'+savename+'.pdf',bbox_inches='tight')



TSplot(5,250,alph=1)
TSplot('5corr',250,alph=1)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF5_TScorr'+savename+'.pdf',bbox_inches='tight')

#####################################################################
### Next, deal with a giant drift in salinity at 150m at CF6 #####
#####################################################################


figure(figsize=(12,4))
for key in sal[6]:
    plot(date[6][key],sal[6][key],'-',label=key)
legend(loc=(1.05,0.2))

figure(figsize=(10,4))
plot(date[6][150],sal[6][150]-sal[6][100][:-3])


def removeline(tvec,svec,anchor=0):
    linefunc=poly1d(polyfit(tvec,svec, 1))
    liney=linefunc(tvec)
    print(len(svec))
    corrected=svec-liney+liney[anchor]
    return liney,corrected

s100=sal[6][150][~isnan(sal[6][150])]
t100=time[6][150][~isnan(sal[6][150])]

line6_tot,scorr6_tot=removeline(t100,s100)

f=interpolate.interp1d(t100,scorr6_tot)
scorr6_final=f(time[6][150])

plot(time[6][150],scorr6_final)
plot(time[6][150],sal[6][150])
plot(t100,line6_tot)



sal['6corr_tot']={}
sal['6corr_tot'][150]=scorr6_final
tmp['6corr_tot']={}
tmp['6corr_tot'][150]=tmp[6][150]

len(sal['6corr_tot'][150])
len(tmp['6corr_tot'][150])

figure()
TSplot(6,100)
TSplot(6,300)
TSplot(6,150)
spruce_TS()
xlim([33.5,35.5])
savefig('../notes/figures/CF6_TS'+savename+'.pdf',bbox_inches='tight')

figure()
TSplot(6,100)
TSplot(6,300)
TSplot('6corr_tot',150)
spruce_TS()
xlim([33.5,35.5])
savefig('../notes/figures/CF6_TScorr'+savename+'.pdf',bbox_inches='tight')


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
# ylim([34.5,35.1])
savefig('../notes/figures/CF6_salcorr'+savename+'.pdf',bbox_inches='tight')

#####################################################################
###### Similar drift spotted in CF4 ##################################
#####################################################################


figure(figsize=(10,3))
for key in sal[4]:
    plot(date[4][key],sal[4][key],'-',label=key)
legend(loc=(1.05,0.5))


saldif4=sal[4][100]-sal[4][50]
driftime=24050
valitime=24285


# This drift appears after one year, apply linear correction thereafter

timeind4=(time[4][100]>driftime) & (time[4][100]<valitime)

salcorr4=saldif4[timeind4]
tcorr4=time[4][100][timeind4]

timeind4_last=(time[4][100]>=valitime)

liney4,corr4=removeline(tcorr4[~isnan(salcorr4)],salcorr4[~isnan(salcorr4)])


plot(time[4][50],saldif4)
axvline(driftime,color='r')
plot(tcorr4[~isnan(salcorr4)],corr4)
plot(tcorr4[~isnan(salcorr4)],liney4)
axvline(valitime,color='r')
grid('on')



figure(figsize=(14,3))
plot(time[4][100],sal[4][100])
axvline(driftime,color='r')
axvline(valitime,color='r')
plot(tcorr4[~isnan(salcorr4)],corr4+sal[4][50][timeind4][~isnan(salcorr4)])


figure(figsize=(12,3))
plot(time[4][50],saldif4)
axvline(driftime,color='k')
axvline(valitime,color='k')
plot(tcorr4[~isnan(salcorr4)],liney4,'k')
grid('on')
title('CF4 salinity correction')
savefig('../notes/figures/CF4_sdiff'+savename+'.pdf',bbox_inches='tight')


sal['4corr']={}
tmp['4corr']={}
sal['4corr'][100]=sal[4][100].copy()
tmp['4corr'][100]=tmp[4][100].copy()
sal['4corr'][100][time[4][100]>valitime]=NaN

subind4=arange(len(sal[4][100]))[timeind4][~isnan(salcorr4)]

new4corr=(corr4+sal[4][50][subind4])

sal['4corr'][100][subind4]=new4corr



figure(figsize=(14,3))
for key in [50,100]:
    plot(time[4][key],sal[4][key],label=key)
plot(time[4][100],sal['4corr'][100],label='100, corrected')
axvline(driftime,color='k')
axvline(valitime,color='k')
legend(loc=(1.05,0.5),title='nominal pressure')
xlabel('time (days)')
ylabel('salinity')
title('CF4 salinity correction')
savefig('../notes/figures/CF4_salcorr'+savename+'.pdf',bbox_inches='tight')




date['4corr']={}
date['4corr'][100]=date[4][100]
month['4corr']={}
month['4corr'][100]=[dd.month for dd in date['4corr'][100]]


figure()
TSplot_seasonal(4,50)
TSplot_seasonal(4,100)
TSplot_seasonal(4,200)
spruce_seasonal()
title('Before correction')
savefig('../notes/figures/CF4_seasonalTS'+savename+'.pdf',bbox_inches='tight')

figure()
TSplot_seasonal(4,50)
TSplot_seasonal('4corr',100)
TSplot_seasonal(4,200)
spruce_seasonal()
title('After correction')
savefig('../notes/figures/CF4_seasonalTScorr'+savename+'.pdf',bbox_inches='tight')



figure()
TSplot_seasonal(4,100)
spruce_seasonal()

figure()
TSplot_seasonal('4corr',100)
spruce_seasonal()


figure()
TSplot(4,50)
TSplot(4,200)
TSplot(4,100)
spruce_TS()
xlim([32.5,35.5])
# ylim([0,4])
savefig('../notes/figures/CF4_TS'+savename+'.pdf',bbox_inches='tight')


figure()
TSplot(4,50)
TSplot(4,200)
TSplot('4corr',100)
spruce_TS()
xlim([32.5,35.5])
# ylim([0,4])
savefig('../notes/figures/CF4_TScorr'+savename+'.pdf',bbox_inches='tight')

#####################################################################
###################### More drift at CF7! #########################
#####################################################################


for key in sort(list(sal[7].keys()))[:5]:
    plot(date[7][key],sal[7][key],label=key)
legend(loc=(1.05,0.5))



TSplot(6,100)
TSplot(7,100)
spruce_TS()



TSplot(6,300)
# TSplot(8,300)
TSplot(7,250)
spruce_TS()


offstime=180+7.358e5


sal['7corr']={}
tmp['7corr']={}
for rr in arange(0.01,0.11,0.01):
    sal['7corr'][rr]={}
    tmp['7corr'][rr]={}
    sal['7corr'][rr]=sal[7][100].copy()
    tmp['7corr'][rr]=tmp[7][100].copy()
    sal['7corr'][rr][time[7][100]>offstime]=sal[7][100].copy()[time[7][100]>offstime]+rr


figure()
TSplot(6,100)
TSplot(7,100)
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF7_TS_CF6comp'+savename+'.pdf',bbox_inches='tight')

figure()
TSplot(6,100)
# TSplot(6,300)
# TSplot('6corr_tot',150)
# TSplot(7,50)
# TSplot(7,250)
TSplot('7corr',sort(list(sal['7corr'].keys()))[-4])
spruce_TS()
xlim([34,35.5])
savefig('../notes/figures/CF7_TScorr_CF6comp'+savename+'.pdf',bbox_inches='tight')


for rr in arange(0.01,0.11,0.01):
    figure()
    TSplot(7,50)
    TSplot(7,250)
    TSplot('7corr',rr)
    spruce_TS()

figure(figsize=(14,3))
for key in sort(list(sal[7].keys()))[:5]:
    plot(time[7][key],sal[7][key],label=key)
plot(time[7][100],sal['7corr'][sort(list(sal['7corr'].keys()))[-4]],label='100, corrected')
axvline(offstime,color='k')
legend(loc=(1.05,0.25),title='nominal pressure')
xlabel('time (days)')
ylabel('salinity')
title('CF7 salinity correction')
savefig('../notes/figures/CF7_salcorr'+savename+'.pdf',bbox_inches='tight')

figure()
TSplot(6,100)
# TSplot('6corr',150)
TSplot('6corr_tot',150)
TSplot(6,300)
spruce_TS()


TSplot(7,50)
TSplot(7,250)
TSplot(7,100)
spruce_TS()
xlim([34.5,35.2])
savefig('../notes/figures/CF7_TS'+savename+'.pdf',bbox_inches='tight')


TSplot(7,50)
TSplot(7,250)
TSplot('7corr',sort(list(sal['7corr'].keys()))[-4])
spruce_TS()
xlim([34.5,35.2])
savefig('../notes/figures/CF7_TScorr'+savename+'.pdf',bbox_inches='tight')


#####################################################################
# Create corrected sal and tmp dictionaries to reflect changes
#####################################################################

salfinal={}
tmpfinal={}
datefinal={}
for key in range(4,8):
    salfinal[key]={}
    tmpfinal[key]={}
    for key2 in sal[key]:
        salfinal[key][key2]=sal[key][key2]
        tmpfinal[key][key2]=tmp[key][key2]


salfinal[4][100]=sal['4corr'][100]
tmpfinal[4][100]=tmp['4corr'][100]

salfinal[5][250]=sal['5corr'][250]
tmpfinal[5][250]=tmp['5corr'][250]

salfinal[6][150]=sal['6corr_tot'][150]
tmpfinal[6][150]=tmp['6corr_tot'][150]

salfinal[7][100]=sal['7corr'][sort(list(sal['7corr'].keys()))[-4]]
tmpfinal[7][100]=tmp['7corr'][sort(list(sal['7corr'].keys()))[-4]]

len(sal[4][100])
len(salfinal[4][100])

len(sal[5][250])
len(salfinal[5][250])

len(sal[6][150])
len(salfinal[6][150])

len(sal[7][100])
len(salfinal[7][100])

date.keys()

del date['4corr']
del month['4corr']

salfinal.keys()

pickle.dump([date,month,prs,salfinal,tmpfinal],
            open('../pickles/TSinterp/TS_daily_dic_wcorr_nofilter.pickle','wb'))
