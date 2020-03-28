
from aux_funcs import *


########################################################################################
############################ 2005,2006 and 2008 cruises ##################################
########################################################################################

## Penny's bathymetry file
pennybathdat=io.loadmat(datadir+'Shipboard/jr302_osnap_sim_bathfromPenny.mat')['osnapnogap'][0][0]
pennylat=pennybathdat[0].flatten()
pennylon=pennybathdat[1].flatten()
pennybath=pennybathdat[2].flatten()
pennydist=pennybathdat[4].flatten()

plot(pennylon,pennylat);

plot(pennydist,pennybath);

#dtype=[('lat', 'O'), ('lon', 'O'), ('stn', 'O'), ('maxpress', 'O'), ('press', 'O'), ('temp', 'O'), ('potemp', 'O'), ('cond', 'O'), ('sal', 'O'), ('oxy', 'O'), ('trans', 'O'), ('sig0', 'O'), ('dist', 'O')])
datadir

alldat=sort(glob.glob(datadir+'Shipboard/CTD_fromPenny/*/*'))
alldat


dicvec=['2005','2006','2008','2014','2015','2016']
lon={}
lat={}
prs={}
sal={}
tmp={}
date={}
for ii,aa in enumerate(alldat[:3]):
    print(aa)
    dat=io.loadmat(aa)['ctd'][0][0]
    lat[dicvec[ii]]=dat[0].flatten()
    lon[dicvec[ii]]=dat[1].flatten()
    prstmp=sort(unique(dat[4]))
    prstmp=prstmp[~isnan(prstmp)]
    prs[dicvec[ii]]=prstmp
    sal[dicvec[ii]]=dat[8][:len(prstmp),:]
    tmp[dicvec[ii]]=dat[5][:len(prstmp),:]


colvec=['r','b','orange','grey','c','purple']
for ii,key in enumerate(sal):
    figure(1)
    plot(lon[key],lat[key],'o',color=colvec[ii])
    # plot(sal[key],tmp[key],'o',color=colvec[ii])
    print(shape(sal[key]),shape(prs[key]),shape(lat[key]))
    figure()
    contourf(lon[key],prs[key],sal[key])
    colorbar()



########################################################################################
###################### 2014 dedicated hydrography cruise ###############################
########################################################################################
#For JR302 (alldat[3]['ctd']) 2014
#dtype=[('lat', 'O'), ('lon', '1'), ('stn', '2'), ('yr', '3'), ('month', '4'), ('day', '5'), ('timebotHH', 'O'), ('timebotMM', 'O'), ('maxpress', 'O'), ('press', '9'), ('temp', '1O'), ('potemp', 'O'), ('cond', 'O'), ('sal', '13'), ('oxy', 'O'), ('trans', 'O'), ('fluor', 'O'), ('sig0', 'O'), ('dist', 'O')]),
jrnum=3
dat=io.loadmat(alldat[jrnum])['ctd'][0][0]
dat

lat[dicvec[jrnum]]=dat[0].flatten()
lon[dicvec[jrnum]]=dat[1].flatten()
prstmp=sort(unique(dat[9]))
prstmp=prstmp[~isnan(prstmp)]
prs[dicvec[jrnum]]=prstmp
sal[dicvec[jrnum]]=dat[13][:len(prstmp),:]
tmp[dicvec[jrnum]]=dat[10][:len(prstmp),:]

lonind14=(lon['2014']<-30.5) & (lon['2014']>-44)
lat['2014']=lat['2014'][lonind14]
lon['2014']=lon['2014'][lonind14]
prs['2014']=prs['2014']
sal['2014']=sal['2014'][:,lonind14]
tmp['2014']=tmp['2014'][:,lonind14]

########################################################################################
################################## 2014 KN221 section ########################################
########################################################################################

dates={}

ctd_kn221=sort(glob.glob(datadir+'Shipboard/kn221_2014/cnv_files/k3/*.cnv'))

prs['kn221']=arange(2,3172,2) # determined in retrospect after inspecting
lat['kn221']=zeros(len(ctd_kn221))
lon['kn221']=zeros(len(ctd_kn221))
sal['kn221']=zeros((len(prs['kn221']),len(ctd_kn221)))
tmp['kn221']=zeros((len(prs['kn221']),len(ctd_kn221)))
for ii,dd in enumerate(ctd_kn221):
    profile = cnv.fCNV(dd)
    lat['kn221'][ii]=profile.attributes['LATITUDE']
    lon['kn221'][ii]=profile.attributes['LONGITUDE']
    print(profile.attributes['LONGITUDE'],profile.attributes['datetime'])
    salfnc=interpolate.interp1d(profile['PRES'],profile['PSAL'],kind='linear',fill_value='NaN',bounds_error=False)
    tmpfnc=interpolate.interp1d(profile['PRES'],profile['TEMP'],kind='linear',fill_value='NaN',bounds_error=False)
    sal['kn221'][:,ii]=salfnc(prs['kn221'])
    tmp['kn221'][:,ii]=tmpfnc(prs['kn221'])

########################################################################################
################################## 2015 section ########################################
########################################################################################
# #For PE 400
# #dtype=[('JD0', 'O'), ('T', '1'), ('S', '2'), ('C', 'O'), ('time0', 'O'), ('datetime', 'O'), ('P', '6'), ('file', 'O'), ('lon', '8'), ('lat', '9'), ('dlon', '10'), ('dlat', '11')])}
penum=4


dat=io.loadmat(datadir+'Shipboard/CTD_fromPenny/PE400 CTD data/pe400_Ar7E_CTD_section.mat')
# dat=io.loadmat(datadir+'Shipboard/CTD_fromPenny/PE400 CTD data/PE400_1Hz.mat')['ctd_cal']

dat['dd']

datelist_2015=[datetime.datetime(1,1,1)+datetime.timedelta(days=tt-366) for tt in dat['dd'][0]]
# datelist_2015=[datetime.datetime(1,1,1)+datetime.timedelta(days=tt-366) for tt in dat['time0'][0][0].flatten()]

datelist_2015

lat[dicvec[penum]]=dat['LAT'].flatten()
lon[dicvec[penum]]=dat['LON'].flatten()

plot(lon['2015'], datelist_2015,'o')


prs[dicvec[penum]]=dat['P'].flatten()
sal[dicvec[penum]]=dat['S']
tmp[dicvec[penum]]=dat['T']
shape(sal['2015'])
shape(prs['2015'])
shape(lon['2015'])

plot(dat['POTEMP'],'r');
plot(dat['T'],'b');
plot(tmp['2015'],'g');
axhline(2)


contourf(lon['2015'],prs['2015'],tmp['2015'])
colorbar()

# lat[dicvec[penum]]=[NaN]
# lon[dicvec[penum]]=[NaN]
# sal[dicvec[penum]]=[NaN]
# prs[dicvec[penum]]=[NaN]
# tmp[dicvec[penum]]=[NaN]

########################################################################################
################################## 2016 section ########################################
########################################################################################
# name 0 = prDM: Pressure, Digiquartz [db]
# name 1 = depSM: Depth [salt water, m]
# name 2 = t090C: Temperature [ITS-90, deg C]
# name 3 = t190C: Temperature, 2 [ITS-90, deg C]
# name 4 = c0mS/cm: Conductivity [mS/cm]
# name 5 = c1mS/cm: Conductivity, 2 [mS/cm]
# name 6 = sal00: Salinity, Practical [PSU]
# name 7 = sal11: Salinity, Practical, 2 [PSU]
# name 8 = sbeox0ML/L: Oxygen, SBE 43 [ml/l]
# name 9 = sbeox0V: Oxygen raw, SBE 43 [V]
# name 10 = flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]
# name 11 = CStarTr0: Beam Transmission, WET Labs C-Star [%]
# name 12 = v0: Voltage 0
# name 13 = v2: Voltage 2
# name 14 = v4: Voltage 4
# name 15 = v6: Voltage 6
# name 16 = svCM: Sound Velocity [Chen-Millero, m/s]
# name 17 = altM: Altimeter [m]
# name 18 = scan: Scan Count
# name 19 = timeS: Time, Elapsed [seconds]
# name 20 = flag: flag
########################################################################################
ctdlist_16=sort(glob.glob(datadir+'Shipboard/ar07_2016/dcc_files/*'))

arnum=5
prs['2016']=arange(1,3180,2) # determined in retrospect after inspecting
lat[dicvec[arnum]]=[]
lon[dicvec[arnum]]=[]
sal[dicvec[arnum]]=zeros((len(prs['2016']),len(ctdlist_16)))
tmp[dicvec[arnum]]=zeros((len(prs['2016']),len(ctdlist_16)))
for ii,cc in enumerate(ctdlist_16):
    for line in open(cc):
        if 'Latitude:' in line:
            lattmp=float(line[line.find('Latitude:')+9:line.find('Longitude:')])
            lontmp=float(line[line.find('Longitude:')+10:line.find('Date:')])
            print(lontmp,line[line.find('Date:'):])
    lat[dicvec[arnum]]=hstack((lat[dicvec[arnum]],lattmp))
    lon[dicvec[arnum]]=hstack((lon[dicvec[arnum]],lontmp))

    datpd=pd.read_csv(cc,header=2,sep='\s*')
    salfunc=interpolate.interp1d(datpd.iloc[:,0],datpd.iloc[:,3],kind='linear',fill_value='NaN',bounds_error=False)
    sal[dicvec[arnum]][:,ii]=salfunc(prs['2016'])
    tmpfunc=interpolate.interp1d(datpd.iloc[:,0],datpd.iloc[:,1],kind='linear',fill_value='NaN',bounds_error=False)
    tmp[dicvec[arnum]][:,ii]=tmpfunc(prs['2016'])

########################################################################################
################################## Grid all ########################################
########################################################################################

for key in lon:
    sort00=argsort(lon[key])
    lon[key]=lon[key][sort00]
    lat[key]=lat[key][sort00]
    sal[key]=sal[key][:,sort00]
    tmp[key]=tmp[key][:,sort00]


sal.keys()

# Find distance from Jo's zero point
################################################################################################
lon0=-43.0884
lat0=60.0911

def getdist(lonex,latex):
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([lat0,latex[ii]],[lon0,lonex[ii]])[0][0]
    return distout

dist={}
for key in sal:
    dist[key]=getdist(lon[key],lat[key])

for key in sal:
    plot(lon[key],lat[key],'o-',label=key)
legend()

for key in sal:
    plot(dist[key],'o',label=key)
legend()

for key in sal:
    plot(diff(dist[key]),'o',label=key)
legend()
ylim([-10,20])

for key in sal:
    print(shape(sal[key]))

dist['kn221']

# Grid onto same grid as LADCP (10db x 10km)

smprs=range(10, 3125, 10)
distgrid=range(0,int(max(dist['2014']))+10,10)
distlen=len(distgrid)
smprslen=len(smprs)
def gridel(vex,distex,elkey):
    origdist=shape(vex)[1]

    # 1. Smooth each dataset into 10db bins in the vertical (use above and below measurement for each)
    vsm=zeros((smprslen,origdist))
    for jj in range(origdist):
        for ii in range(smprslen):
            if elkey!='2015':
                vsm[ii,jj]=mean(vex[5*ii:5*ii+5,jj])
            else:
                vsm[ii,jj]=mean(vex[10*ii:10*ii+10,jj])

    # 2. Extend each profile downwards
    vext=vsm.copy()
    # for jj in range(origdist):
    #     vext[isnan(vsm[:,jj]),jj]=vsm[~isnan(vsm[:,jj]),jj][-1]

    # 3. Interpolate horizontally onto a uniform 10km dist grid
    vhorz=zeros((smprslen,distlen))
    for ii in range(smprslen):
            distfunc=interpolate.interp1d(distex,vext[ii,:],kind='linear',fill_value='NaN',bounds_error=False)
            vhorz[ii,:]=distfunc(distgrid)

    return vhorz


for key in dist:
    print(shape(sal[key]))

newdic=['2005', '2006','2008', '2014','kn221','2015', '2016']
# newdic=['2005', '2008', '2014','kn221', '2016']


saltot=NaN*ones((smprslen,distlen,7))
tmptot=NaN*ones((smprslen,distlen,7))

for ii,key in enumerate(newdic):
        print(key)
        saltot[:,:,ii]=gridel(sal[key],dist[key],key)
        tmptot[:,:,ii]=gridel(tmp[key],dist[key],key)


pdentot=NaN*ones(shape(saltot))
for kk in range(7):
    print(kk)
    for jj in range(distlen):
        salnan=~isnan(saltot[:,jj,kk])
        if sum(salnan)>2:
            SA=gsw.SA_from_SP(saltot[:,jj,kk][salnan],array(smprs)[salnan],[lon0]*sum(salnan),[lat0]*sum(salnan))
            pdentot[salnan,jj,kk]=gsw.sigma0(SA,gsw.CT_from_t(SA,tmptot[:,jj,kk][salnan],array(smprs)[salnan]))

occvec=['2005', '2006','2008', '2014 (JR302)','2014 (KN221)','2015', '2016']
# occvec=['2005', '2008', '2014 (JR302)','2014 (KN221)', '2016']

#Get into an xarray
grdat=xr.Dataset({'salinity': (['pressure [db]','distance [km]', 'occupation'],  saltot),
                    'temperature [$^\circ$ C]': (['pressure [db]','distance [km]', 'occupation'],  tmptot),
                    'density': (['pressure [db]','distance [km]', 'occupation'],  pdentot),},
                coords={'distance [km]': distgrid,
                        'pressure [db]': smprs,
                        'occupation': occvec})

grdat['density'].sel(occupation='2016').plot()


pickle.dump(grdat,open('../pickles/Shipboard/CTD_xarray_1810_forJo.pickle','wb'))

########################################################################################
################################## Plot! ########################################
########################################################################################
d1=1
d2=30

fullbathdist=bathy['dist'].flatten()
fullbathbath=bathy['bath'].flatten()

g=grdat['temperature [$^\circ$ C]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=0,vmax=8,cmap=cm.RdBu_r,figsize=(15,12))
ylim([3e3,0]);
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    kcont=ax.contour(distgrid[d1:d2],smprs,grdat['temperature [$^\circ$ C]'][:,d1:d2,ii],range(-2,12,1),colors='k')
    clabel(kcont,fmt='%d')
savefig('../figures/CTD_cruises/firstsix_extdpth_tmp.png',bbox_inches='tight')
savefig('../figures/CTD_cruises/firstsix_extdpth_tmp.pdf',bbox_inches='tight')

figure(figsize=(10,6))
contourf(distgrid[d1:d2],smprs,grdat['temperature [$^\circ$ C]'][:,d1:d2,:].mean(dim='occupation'),1001,vmin=0,vmax=8,cmap=cm.RdBu_r)
colorbar(ticks=arange(0,10,2),label='[$^\circ$ C]')
kcont=contour(distgrid[d1:d2],smprs,grdat['temperature [$^\circ$ C]'][:,d1:d2,:].mean(dim='occupation'),range(-2,12,1),colors='k')
clabel(kcont,fmt='%d')
fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
title('Mean temperature from CTD sections',fontsize=20)
xlabel('distance [km]')
ylabel('pressure [db]')
ylim([3e3,0])
xlim([0,300])
savefig('../figures/CTD_cruises/firstsix_mean_extdpth_tmp.png',bbox_inches='tight')
savefig('../figures/CTD_cruises/firstsix_mean_extdpth_tmp.pdf',bbox_inches='tight')
(44,162,95)
colors = [(44,127,184),(35,139,69) ,(237,248,177),(254,196,79),(217,95,14)] # This example uses the 8-bit RGB
sal_cmap = make_cmap(colors,position=[0,0.65,0.85,0.95,1],bit=True)


g=grdat['salinity'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation',vmin=34,vmax=35.1, col_wrap=2,cmap=sal_cmap,figsize=(15,12))
ylim([3e3,0]);
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    kcont=ax.contour(distgrid[d1:d2],smprs,grdat['salinity'][:,d1:d2,ii],[34,34.9,35],colors='k')
    clabel(kcont,fmt='%1.1f')
savefig('../figures/CTD_cruises/firstsix_extdpth_sal.png',bbox_inches='tight')
savefig('../figures/CTD_cruises/firstsix_extdpth_sal.pdf',bbox_inches='tight')


figure(figsize=(10,6))
contourf(distgrid[d1:d2],smprs,grdat['salinity'][:,d1:d2,:].mean(dim='occupation'),linspace(34,35.1,1001),cmap=sal_cmap,extend="both")
colorbar(ticks=arange(34,35.1,0.25))
clim(34,35.1)
kcont=contour(distgrid[d1:d2],smprs,grdat['salinity'][:,d1:d2,:].mean(dim='occupation'),[34,34.9,35],colors='k')
clabel(kcont,fmt='%1.1f')
fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
title('Mean salinity from CTD sections',fontsize=20)
xlabel('distance [km]')
ylabel('pressure [db]')
ylim([3e3,0])
xlim([0,300])
savefig('../figures/CTD_cruises/firstsix_mean_extdpth_sal.png',bbox_inches='tight')
savefig('../figures/CTD_cruises/firstsix_mean_extdpth_sal.pdf',bbox_inches='tight')

colors = [(237,248,177),(127,205,187),(5,112,176),(110,1,107)]
pden_cmap=make_cmap(colors,position=[0,0.666666666666,0.833333333333,1],bit=True)

g=grdat['density'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation',vmin=27,vmax=27.9, col_wrap=2,cmap=pden_cmap,figsize=(15,12))
ylim([3e3,0]);
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    kcont=ax.contour(distgrid[d1:d2],smprs,grdat['density'][:,d1:d2,ii],[27.2,27.6,27.7,27.8,27.9],colors='k')
    clabel(kcont,fmt='%1.1f')
savefig('../figures/CTD_cruises/firstsix_extdpth_pden.png',bbox_inches='tight')
savefig('../figures/CTD_cruises/firstsix_extdpth_pden.pdf',bbox_inches='tight')


figure(figsize=(10,6))
colorcont=contourf(distgrid[d1:d2],smprs,grdat['density'][:,d1:d2,:].mean(dim='occupation'),linspace(27,27.9,1001),cmap=pden_cmap,extend="both")
colorbar(colorcont,ticks=arange(27,28,0.2))
kcont=contour(distgrid[d1:d2],smprs,grdat['density'][:,d1:d2,:].mean(dim='occupation'),[27.2,27.6,27.7,27.8,27.9],colors='k')
clabel(kcont,fmt='%1.1f')
fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
title('Mean density from CTD sections',fontsize=20)
xlabel('distance [km]')
ylabel('pressure [db]')
ylim([3e3,0])
xlim([0,300])
savefig('../figures/CTD_cruises/firstsix_mean_extdpth_pden.png',bbox_inches='tight')
savefig('../figures/CTD_cruises/firstsix_mean_extdpth_pden.pdf',bbox_inches='tight')


########################################################################################
######################### 2014 KN221 sub-sections #######################################
########################################################################################


def colstopanda(keylab):
    field=pd.DataFrame(data=ctd_2014[keylab],columns=['distance','depth',keylab])
    field=field.pivot(columns='distance',index='depth')
    del field.index.name
    return field

sal_14={}
tmp_14={}
pden_14={}
bathdist_14={}
bathbath_14={}
for snum in range(1,4):
    ctd_2014=io.loadmat(datadir+'Shipboard/kn221_2014/mat_files/k'+str(snum)+'.mat')
    sal_14[snum]=colstopanda('k'+str(snum)+'_salt')
    tmp_14[snum]=colstopanda('k'+str(snum)+'_ptmp')
    pden_14[snum]=colstopanda('k'+str(snum)+'_dens')
    bathdist_14[snum]=ctd_2014['k'+str(snum)+'_bath'][:,0]
    bathbath_14[snum]=ctd_2014['k'+str(snum)+'_bath'][:,1]


ctd_2014

def onecont(apanda,tit,vrange,coloor,hlevs,bathdist,bathbath,ylim1,ylim2,xlim1,xlim2,savename):
    # figure()
    ax1=contourf(apanda.columns.levels[1],apanda.index,apanda.values,vrange,cmap=coloor)
    colorbar()
    contour(apanda.columns.levels[1].values,apanda.index,apanda.values,levels=hlevs,colors='k')
    fill_between(bathdist,bathbath,3000*ones(len(bathbath)),color='k',zorder=22)
    gca().invert_yaxis()
    ylim([ylim2,ylim1])
    xlim([xlim1,xlim2])
    xlabel('distance [km]')
    ylabel('pressure [db]')
    title(tit,fontsize=18)
    savefig('../figures/CTD/'+savename+'.png')

    return ax1,vrange


onecont(sal_14[1],'Salinity on k1 section of k221, 2014',linspace(26,35.5,21),cm.YlGnBu_r,[34.0,34.4,34.8],bathdist_14[1],bathbath_14[1],0,400,-2,60,'Sal_2014_1')
onecont(sal_14[2],'Salinity on k2 section of k221, 2014',linspace(28,35.5,21),cm.YlGnBu_r,[34.0,34.4,34.8],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Sal_2014_2')
onecont(sal_14[3],'Salinity on k3 section of k221, 2014',linspace(28,35.5,21),cm.YlGnBu_r,[34.0,34.4,34.8],bathdist_14[3],bathbath_14[3],0,400,-2,60,'Sal_2014_3')

onecont(sal_14[3],'Salinity on k3 section of k221, 2014',linspace(28,35.5,21),cm.YlGnBu_r,[34.0,34.4,34.8],bathdist_14[3],bathbath_14[3],0,2500,-2,160,'Sal_2014_3_zoomout')

onecont(tmp_14[1],'Temperature on k1 section of k221, 2014',linspace(-2,10,13),cm.RdYlBu_r,[3],bathdist_14[1],bathbath_14[1],0,400,-2,60,'Tmp_2014_1')
onecont(tmp_14[2],'Temperature on k2 section of k221, 2014',linspace(-2,10,13),cm.RdYlBu_r,[3],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Tmp_2014_2')
onecont(tmp_14[3],'Temperature on k3 section of k221, 2014',linspace(-2,10,13),cm.RdYlBu_r,[3],bathdist_14[3],bathbath_14[3],0,400,-2,60,'Tmp_2014_3')


onecont(pden_14[1],'Density on k1 section of k221, 2014',linspace(23,28,21),cm.RdPu,[25,26,27],bathdist_14[1],bathbath_14[1],0,400,-2,60,'Pden_2014_1')
onecont(pden_14[2],'Density on k2 section of k221, 2014',linspace(23,28,21),cm.RdPu,[25,26,27],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Pden_2014_2')
onecont(pden_14[3],'Density on k3 section of k221, 2014',linspace(23,28,21),cm.RdPu,[25,26,27],bathdist_14[3],bathbath_14[3],0,400,-2,60,'Pden_2014_3')

onecont(pden_14[2],'Density on k2 section of k221, 2014',linspace(23,28,21),cm.RdPu,[25,26,27],bathdist_14[2],bathbath_14[2],0,2500,-2,160,'Pden_2014_2_zoomout')
onecont(pden_14[3],'Density on k3 section of k221, 2014',linspace(23,28,21),cm.RdPu,[25,26,27],bathdist_14[3],bathbath_14[3],0,2500,-2,160,'Pden_2014_3_zoomout')

linspace(-2,10,49)

figure(figsize=(14,3))
subplot(131)
onecont(tmp_14[2],'Temperature on k2 section of k221, 2014',linspace(-2,10,49),cm.RdYlBu_r,range(-1,9,2),bathdist_14[2],bathbath_14[2],0,400,-2,60,'Tmp_2014_2')
title('Potential Temperature [$^\circ$C]')
subplot(132)
onecont(sal_14[2],'Salinity on k2 section of k221, 2014',linspace(28,35.5,51),cm.YlGnBu_r,[30,32,34.0,34.4,34.8],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Sal_2014_2')
title('Salinity')
ylabel('')
subplot(133)
onecont(pden_14[2],'Density on k2 section of k221, 2014',linspace(23,28,51),cm.RdPu,[25,26,27,27.5],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Pden_2014_2')
ylabel('')
title('Potential Density [kg/m$^3$]')
savefig('../figures/CTD/K2sections_triptich.pdf',bbox_inches='tight')

adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')


def acrosstrackname(adic):
    adic['across track velocity']=adic['v_dt']*sin(theta)+adic['u_dt']*cos(theta)

    return adic

adcp_14=acrosstrackname(adcp_14)

def contadcp(distvec,adic,tit):
    contourf(distvec,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.8,0.05),cmap=cm.RdBu_r);
    colorbar()
    contour(distvec,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2],colors='k');
    fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    title('Velocity along bathymetry [m/s]')
    ylim([400,0])
    xlim([-10,60])
    xlabel('distance [km]')

adcp_dist=array(pd.DataFrame.from_csv(datadir+'Shipboard/adcp_distances.dat').index)

with open(datadir+'Shipboard/adcp_distances.dat', 'r') as f:
    reader = csv.reader(f)
    adcp_dist = list(reader)

adcp_dist=[float(dd[0]) for dd in adcp_dist]

figure(figsize=(10,8))
subplot(221)
onecont(tmp_14[2],'Temperature on k2 section of k221, 2014',linspace(-2,10,49),cm.RdYlBu_r,range(-1,9,2),bathdist_14[2],bathbath_14[2],0,400,-2,60,'Tmp_2014_2')
title('Potential Temperature [$^\circ$C]')
xlabel('')
subplot(222)
onecont(sal_14[2],'Salinity on k2 section of k221, 2014',linspace(28,35.5,51),cm.YlGnBu_r,[30,32,34.0,34.4,34.8],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Sal_2014_2')
title('Salinity')
ylabel('')
xlabel('')
subplot(223)
onecont(pden_14[2],'Density on k2 section of k221, 2014',linspace(23,28,51),cm.RdPu,[25,26,27,27.5],bathdist_14[2],bathbath_14[2],0,400,-2,60,'Pden_2014_2')
title('Potential Density [kg/m$^3$]')
subplot(224)
contadcp(adcp_dist,adcp_14,'2014')
savefig('../figures/CTD/K2sections_4plusvel.pdf',bbox_inches='tight')


def vert_turner(salarray,prsvec,tmparray,pdenarray):

    # get vertical turner angles for k2 section...
    SA=salarray.copy()
    for ii in range(shape(salarray)[1]):
        SA[:,ii]=gsw.SA_from_SP(salarray[:,ii],prsvec,CFlon[3],CFlat[3]).T

    CT=gsw.CT_from_pt(SA,tmparray)

    prsint=prsvec[:-1]+diff(prsvec)

    SAint=SA[:-1,:]+diff(SA,axis=0)/2
    CTint=CT[:-1,:]+diff(CT,axis=0)/2
    alphaint=CTint.copy()
    betaint=SAint.copy()
    for ii in range(shape(salarray)[1]):
        alphaint[:,ii]=gsw.alpha(SAint[:,ii].T,CTint[:,ii].T,prsint).T
        betaint[:,ii]=gsw.beta(SAint[:,ii].T,CTint[:,ii].T,prsint).T

    alphaT=-alphaint*diff(tmparray,axis=0)
    betaS=betaint*diff(salarray,axis=0)

    R=alphaT/betaS

    dendiff=diff(pdenarray,axis=0)

    turner=arctan2((alphaT-betaS)*dendiff/abs(dendiff),(alphaT+betaS)*dendiff/abs(dendiff))*180/pi

    turner_noden=arctan2((alphaT-betaS),(alphaT+betaS))*180/pi

    return alphaint,betaint,alphaT,betaS,dendiff,turner,turner_noden

alphaint,betaint,alphaT,betaS,dendiff,turner,turner_noden=vert_turner(sal_14[2].values,sal_14[2].index,tmp_14[2],pden_14[2])
pint=sal_14[2].index[:-1]+diff(sal_14[2].index)/2

plot(alphaT,pint,'r.');
ylim([400,0])

plot(betaS,pint,'b.');
ylim([400,0])

plot(alphaint,'r');
plot(betaint,'b');

pint
shape(turner)

dd=2
[plot(turner[::dd,ii],pint[::dd],'.') for ii in range(0,21,4)];
[plot(turner[::dd,ii],pint[::dd],'-',alpha=0.7) for ii in range(0,21,4)];
ylim([400,0])
axvline(0,color='k')
text(10,-15,'Temperature dominated',color='red')
text(-85,-15,'Salinity dominated',color='blue')
xticks(arange(-90,100,45))
xlim([-100,100])
xlabel('Turner Angle [$^\circ$]')
ylabel('pressure [db]')
savefig('../figures/shipboard_turb/vert_turner_2014.pdf',bbox_inches='tight')


# plot(turner_noden,,'.');

shape(sal_14[2])

########################################################################################
##################### Look at CTD profile detail for turbulence #######################
########################################################################################

pden={}
for occ in ['2005','2008','2014','2016']:
    pden[occ]=sal[occ].copy()
    for ii in range(shape(sal[occ])[1]):
        SA=gsw.SA_from_SP(sal[occ][:,ii],prs[occ],[lon0]*len(prs[occ]),[lat0]*len(prs[occ]))
        pden[occ][:,ii]=gsw.sigma0(SA,gsw.CT_from_t(SA,tmp[occ][:,ii],prs[occ]))

shape(sal['2005'])
shape(dist['2005'])
dist['2005']<60

for occ in ['2005','2008','2014','2016']:
    figure(figsize=(10,3))
    subplot(121)
    plot(sal[occ][:,dist[occ]<100],prs[occ]);
    ylim([200,0])
    ylabel('pressure [db]')
    xlabel('salinity')
    subplot(122)
    plot(tmp[occ][:,dist[occ]<100],prs[occ]);
    ylim([200,0])
    xlabel('temperature')
    suptitle(occ)
    savefig('../figures/shipboard_turb/TS_profs_'+occ+'.png',bbox_inches='tight')

for occ in ['2005','2008','2014','2016']:
    figure(figsize=(14,3))
    subplot(131)
    plot(sal[occ][:,dist[occ]<100],prs[occ]);
    ylim([200,0])
    ylabel('pressure [db]')
    xlabel('salinity')
    subplot(132)
    plot(tmp[occ][:,dist[occ]<100],prs[occ]);
    ylim([200,0])
    xlabel('potential temperature [$^\circ$ C]')
    subplot(133)
    plot(pden[occ][:,dist[occ]<100],prs[occ]);
    ylim([200,0])
    xlabel('potential density [kg/$m^3$]')
    suptitle(occ)
    savefig('../figures/shipboard_turb/TSD_profs_'+occ+'.png',bbox_inches='tight')

univec['sal'][-2]

occ='2016'
figure(figsize=(14,3))
subplot(131)
contourf(dist[occ][dist[occ]<100],prs[occ],sal[occ][:,dist[occ]<100],101,cmap=univec['sal'][2]);
ylim([400,0])
colorbar()
contour(dist[occ][dist[occ]<100],prs[occ],sal[occ][:,dist[occ]<100],levels=univec['sal'][-2],colors='k');
ylabel('pressure [db]')
xlabel('distance [km]')
subplot(132)
contourf(dist[occ][dist[occ]<100],prs[occ],tmp[occ][:,dist[occ]<100],101,cmap=univec['tmp'][2]);
ylim([400,0])
colorbar()
contour(dist[occ][dist[occ]<100],prs[occ],tmp[occ][:,dist[occ]<100],levels=univec['tmp'][-2],colors='k');
xlabel('distance [km]')
subplot(133)
contourf(dist[occ][dist[occ]<100],prs[occ],pden[occ][:,dist[occ]<100],101,cmap=univec['pden'][2]);
colorbar()
contour(dist[occ][dist[occ]<100],prs[occ],pden[occ][:,dist[occ]<100],levels=univec['pden'][-2],colors='k');
ylim([400,0])
xlabel('distance [km]')
savefig('../figures/CTD/TSDsections_2016.png',bbox_inches='tight')


xxxxxxxxxxxxx
    # subplot(132)
    # plot(tmp[occ][:,dist[occ]<100],prs[occ]);
    # ylim([200,0])
    # xlabel('potential temperature [$^\circ$ C]')
    # subplot(133)
    # plot(pden[occ][:,dist[occ]<100],prs[occ]);
    # ylim([200,0])
    # xlabel('potential density [kg/$m^3$]')
    # suptitle(occ)
    # savefig('../figures/shipboard_turb/TSD_profs_'+occ+'.png',bbox_inches='tight')







########################################################################################
######################### Compare with Nuka and moorings #######################################
########################################################################################

OSNAPgrid=pickle.load(open('../pickles/CF_xarray_gridplot_notid_pdenfix.pickle','rb'))

OSNAPgrid
oshelf=OSNAPgrid.where(OSNAPgrid.distance<45).resample('M',dim='date')

nuka=io.loadmat(datadir+'aux_data/nuka/jamie/osnapsssfromnuka.mat')

shape(nuka['osnapsi'])

purps=['darkgrey','lightgrey','grey']

nukasal=nuka['osnapsi'][:5,4:10]
nukatmp=nuka['osnapti'][:5,4:10]

plot(oshelf.salinity[:,:,ii]);

oshelf

shape(oshelf)


ii=0
plot(oshelf.salinity[:,:,ii],oshelf.temperature[:,:,ii].to_array(),'go');
oshelf

salpd=oshelf.salinity.to_pandas()
tmppd=oshelf.temperature.to_pandas()

plot(salpd[0,:,:],tmppd[0,:,:],'go')

oshelf.salinity[ii,:,:].dropna(dim='distance').dropna(dim='depth')
oshelf.date
sumvec=[0,9,10,11,12,21,22,23]
[oshelf.date[ss] for ss in sumvec]

# [plot(oshelf.salinity[ii,:,:].to_pandas(),oshelf.temperature[ii,:,:].to_pandas(),'go'); for ii in range(24)];

figure()
for ii in range(1,4):
    plot(sal_14[ii],tmp_14[ii],'o',color=purps[ii-1]);
for ii in range(6):
    plot(saltot[:,:,ii],tmptot[:,:,ii],'o',color=colvec[ii],label='');
    plot(saltot[0,0,ii],tmptot[0,0,ii],'o',color=colvec[ii],label=occvec[ii]);
plot(nukasal,nukatmp,'k*',label='')
plot(nukasal[0,0],nukatmp[0,0],'k*',label='nuka may-sep')
legend(loc=(1.05,0.25))
xlabel('salinity')
ylabel('temeperature ($^\circ$ C)')
savefig('../figures/CTDcomp/TS_hydroNuka.pdf',bbox_inches='tight')
[plot(oshelf.salinity[ii,:,:].to_pandas(),oshelf.temperature[ii,:,:].to_pandas(),'go',alpha=0.1) for ii in sumvec];
plot(oshelf.salinity[0,0,0],oshelf.temperature[0,0,0],'go',label='OSNAP may-sep')
legend(loc=(1.05,0.25))
savefig('../figures/CTDcomp/TS_hydroNuka_OSNAP.pdf',bbox_inches='tight')

plot(oshelf.depth[0],oshelf.salinity[0,0,0].to_pandas(),'g',label='OSNAP may-sep')

figure()
for ii in range(1,4):
    plot(sal_14[ii],sal_14[ii].index,'-',color=purps[ii-1]);
for ii in range(6):
    plot(saltot[0,0,ii],grdat['pressure [db]'][0],'-',color=colvec[ii],linewidth=3,label=occvec[ii]);
    plot(saltot[:,:5,ii],grdat['pressure [db]'],'-',color=colvec[ii],linewidth=3,label='');
plot(nukasal,10*ones(shape(nukasal)),'k*',label='')
plot(nukasal[0,0],10,'k*',label='nuka may-sep')
axhline(50,color='k')
legend(loc=(1.05,0.25))
ylim([200,0])
savefig('../figures/CTDcomp/Sprof_hydroNuka.pdf',bbox_inches='tight')
[plot(oshelf.salinity[ii,:,:].to_pandas().T,oshelf.depth,'g',alpha=0.5,label='') for ii in sumvec];
plot(oshelf.salinity[0,0,0].to_pandas(),oshelf.depth[0],'g',label='OSNAP may-sep')
legend(loc=(1.05,0.25))
savefig('../figures/CTDcomp/Sprof_hydroNuka_OSNAP.pdf',bbox_inches='tight')

newdic
figure()
for ii in range(1,4):
    plot(tmp_14[ii],tmp_14[ii].index,'-',color=purps[ii-1]);
for ii in range(6):
    plot(tmptot[0,0,ii],grdat['pressure [db]'][0],'-',color=colvec[ii],linewidth=3,label=occvec[ii]);
    plot(tmptot[:,:5,ii],grdat['pressure [db]'],'-',color=colvec[ii],linewidth=3,label='');
axhline(50,color='k')
plot(nukatmp,10*ones(shape(nukasal)),'k*',label='')
plot(nukatmp[0,0],10,'k*',label='nuka may-sep')
legend(loc=(1.05,0.25))
ylim([200,0])
savefig('../figures/CTDcomp/Tprof_hydroNuka.pdf',bbox_inches='tight')
[plot(oshelf.temperature[ii,:,:].to_pandas().T,oshelf.depth,'g',alpha=0.5,label='') for ii in sumvec];
plot(oshelf.temperature[0,0,0].to_pandas(),oshelf.depth[0],'g',label='OSNAP may-sep')
legend(loc=(1.05,0.25))
savefig('../figures/CTDcomp/Tprof_hydroNuka_OSNAP.pdf',bbox_inches='tight')

########################################################################################
##################### Extract nearby stations for Jo's OOI calib ########################
########################################################################################


## Note, there are likely more CTD stations for k3 than the ones in the above mat files, but not sure its worth extracting those right now, as I don't have corresponding LADCP and I am mostly interested in the shelf properties for my own work.

# Flanking mooring A
#
#                 Deployment 1 : 13-Sep-2014 to 18-Aug-2015 at  59.7668, -39.8425
#
#                 Deployment 2: 21-Aug-2015 to 12-July-2016 at 59.7707, -39.8801
#
#
#
# Flanking mooring B
#
#                 Deployment 1: 18-Sep-2014 to 18-Aug-2015 at 59.7128, -39.3204
#
#                 Deployment 2: 22-Aug-2015 to 13-Jul-2016 at 59.7182, -39.3536

# 2014: KN221 cruise was 5-Aug-2014 to 2-Sep-2014
# 2014: Hydrography only cruise was 9-Jun-2014 to 17-July-2014
# 2015: PE400 cruise was 12-July-2015 to 25-July-2015
# 2016: AR07 cruise 5-Aug-2016 to 3-Sep-2016

lon['FLA1']=-39.8425
lat['FLA1']=59.7668

lon['FLB1']=-39.3204
lat['FLB1']=59.7128

flapos
flbpos
#find closest station to each mooring
flbpos=argmin(abs(lon['2015']-lon['FLB1'])+abs(lat['2015']-lat['FLB1']))
plot(lon['2015'],lat['2015'],'o',label='PE400')
plot(lon['FLA1'],lat['FLA1'],'o',label='FLA')
plot(lon['FLB1'],lat['FLB1'],'o',label='FLB')
plot(lon['2015'][flbpos-1:flbpos+2],lat['2015'][flbpos-1:flbpos+2],'*',label='nearest')
legend()
xlabel('longitude (W)')
ylabel('latitude (N)')
title('PE400 stations nearest to OOI flanking moorings')
savefig('../figures/CTD_cruises/stapos_forJo_OOI.pdf',bbox_inches='tight')
dat=io.loadmat(datadir+'Shipboard/CTD_fromPenny/PE400 CTD data/pe400_Ar7E_CTD_section.mat')

seldat15={}
seldat15['date']=list(dat['dd'][0][flbpos-1:flbpos+2])
seldat15['lat']=list(dat['LAT'].flatten()[flbpos-1:flbpos+2])
seldat15['lon']=list(dat['LON'].flatten()[flbpos-1:flbpos+2])
seldat15['prs']=list(dat['P'].flatten())
seldat15['sal']=dat['S'][:,flbpos-1:flbpos+2]
seldat15['tmp']=dat['T'][:,flbpos-1:flbpos+2]

seldat15

for key in seldat15:
    figure()
    title(key)
    plot(seldat15[key],'o')

io.savemat('../data/Shipboard/forJo/pe400_2015_OOIsubset.mat',seldat15)
testdt=io.loadmat('../data/Shipboard/forJo/pe400_2015_OOIsubset.mat')

datelist_2015[flbpos-1:flbpos+2]
# Nearest stations were on July 18, 21 and 23 -- 2015

for key in testdt:
    if '_' not in key:
        figure()
        title(key)
        plot(testdt[key],'o')




########################################################################################
###################### Historical sections from Joleen ###############################
########################################################################################
matpic=datadir+'Shipboard/hydrography_fromJoleen/inuse/Pickart_CTD.mat'

matseg=datadir+'Shipboard/hydrography_fromJoleen/inuse/SEG_shelf_CTD.mat'
## note: for more data, look at SEG_CTD -- would need to wittle it down to really be usable, though
#dtype=[('station', 'O'), ('lat', 'O'), ('lon', 'O'), ('time', 'O'), ('press', 'O'), ('sal', 'O'), ('sal2', 'O'), ('temp', 'O'), ('temp2', 'O'), ('distance', 'O')])

lon={}
lat={}
prs={}
sal={}
tmp={}
date={}
dicvec=['2001 Pickart','2002 Pickart','2003 Pickart','2004 Pickart']
allpic=io.loadmat(matpic)['Pi_CTD'][0]
for ii,datpic in enumerate(allpic):
    lat[dicvec[ii]]=datpic[1].flatten()
    lon[dicvec[ii]]=datpic[2].flatten()
    timetmp=datpic[3].flatten()
    prstmp=sort(unique(datpic[4]))
    prstmp=prstmp[~isnan(prstmp)]
    prs[dicvec[ii]]=prstmp
    sal[dicvec[ii]]=datpic[5][:len(prstmp),:]
    tmp[dicvec[ii]]=datpic[7][:len(prstmp),:]
    date[dicvec[ii]]=[datetime.datetime(1,1,1)+datetime.timedelta(days=float(tt-366)) for tt in timetmp]

#dtype=[('WOD_cruise', 'O'), ('WOD_station', 'O'), ('WOD_type', 'O'), ('lon', 'O'), ('lat', 'O'), ('bot_depth', 'O'), ('OCL_cruise_number', 'O'), ('originator_cruise', 'O'), ('originator_station', 'O'), ('PI', 'O'), ('time', 'O'), ('depth', 'O'), ('temp', 'O'), ('sal', 'O'), ('num_measure', 'O')])

allseg=io.loadmat(matseg)['SEG_shelf_CTD'][0]

wodicvec=['']
name={}
for ii,datpic in enumerate(allseg[1:]):
    name[ii]=datpic[6][0]
    lat[ii]=datpic[4].flatten()
    lon[ii]=datpic[3].flatten()
    timetmp=datpic[10].flatten()
    prstmp=sort(unique(datpic[11]))
    prstmp=prstmp[~isnan(prstmp)]
    prs[ii]=prstmp
    sal[ii]=datpic[-2][:len(prstmp),:]
    tmp[ii]=datpic[-3][:len(prstmp),:]
    date[ii]=[datetime.datetime(1,1,1)+datetime.timedelta(days=float(tt-366)) for tt in timetmp]

for ii in date:
    print(date[ii][0])

for key in date:
    plot(lon[key],lat[key],'o',label=key)
legend()
lon['2001 Pickart']

for key in date:
    plot(sal[key],tmp[key],'o')


unique(prs[key])

for key in date:
    print(shape(sal[key]),shape(prs[key]),shape(lon[key]))


xxxxxxxxxxxxxxxxxxxx
