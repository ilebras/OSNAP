from aux_funcs import *

alldat=sort(glob.glob(datadir+'Shipboard/CTD_fromPenny/*/*'))

dicvec=['2014','2015','2016']#'2005','2006','2008',
lon={}
lat={}
prs={}
sal={}
tmp={}
date={}
oxy={}
########################################################################################
############################ 2005,2006 and 2008 cruises ##################################
########################################################################################

## Penny's bathymetry file
# pennybathdat=io.loadmat(datadir+'Shipboard/jr302_osnap_sim_bathfromPenny.mat')['osnapnogap'][0][0]
# pennylat=pennybathdat[0].flatten()
# pennylon=pennybathdat[1].flatten()
# pennybath=pennybathdat[2].flatten()
# pennydist=pennybathdat[4].flatten()

# plot(pennylon,pennylat);
#
# plot(pennydist,pennybath);

#dtype=[('lat', 'O'), ('lon', 'O'), ('stn', 'O'), ('maxpress', 'O'), ('press', 'O'), ('temp', 'O'), ('potemp', 'O'), ('cond', 'O'), ('sal', 'O'), ('oxy', 'O'), ('trans', 'O'), ('sig0', 'O'), ('dist', 'O')])

# for ii,aa in enumerate(alldat[:3]):
#     print(aa)
#     dat=io.loadmat(aa)['ctd'][0][0]
#     lat[dicvec[ii]]=dat[0].flatten()
#     lon[dicvec[ii]]=dat[1].flatten()
#     prstmp=sort(unique(dat[4]))
#     prstmp=prstmp[~isnan(prstmp)]
#     prs[dicvec[ii]]=prstmp
#     sal[dicvec[ii]]=dat[8][:len(prstmp),:]
#     tmp[dicvec[ii]]=dat[5][:len(prstmp),:]
#
#
# colvec=['r','b','orange','grey','c','purple']
# for ii,key in enumerate(sal):
#     figure(1)
#     plot(lon[key],lat[key],'o',color=colvec[ii])
#     # plot(sal[key],tmp[key],'o',color=colvec[ii])
#     print(shape(sal[key]),shape(prs[key]),shape(lat[key]))
#     figure()
#     contourf(lon[key],prs[key],sal[key])
#     colorbar()
#




########################################################################################
###################### 2014 dedicated hydrography cruise ###############################
########################################################################################
#For JR302 (alldat[3]['ctd']) 2014
#dtype=[('lat', 'O'), ('lon', '1'), ('stn', '2'), ('yr', '3'), ('month', '4'), ('day', '5'), ('timebotHH', 'O'), ('timebotMM', 'O'), ('maxpress', 'O'), ('press', '9'), ('temp', '1O'), ('potemp', 'O'), ('cond', 'O'), ('sal', '13'), ('oxy', 'O'), ('trans', 'O'), ('fluor', 'O'), ('sig0', 'O'), ('dist', 'O')]),

alldat[3]
dat=io.loadmat(alldat[3])['ctd'][0][0]
jrnum=0
lat[dicvec[jrnum]]=dat[0].flatten()
lon[dicvec[jrnum]]=dat[1].flatten()
prstmp=sort(unique(dat[9]))
prstmp=prstmp[~isnan(prstmp)]

prs[dicvec[jrnum]]=prstmp
sal[dicvec[jrnum]]=dat[13][:len(prstmp),:]
tmp[dicvec[jrnum]]=dat[10][:len(prstmp),:]
oxy[dicvec[jrnum]]=dat[14][:len(prstmp),:]

plot(oxy[dicvec[jrnum]])

lonind14=(lon['2014']<-30.5) & (lon['2014']>-44)
lat['2014']=lat['2014'][lonind14]
lon['2014']=lon['2014'][lonind14]
prs['2014']=prs['2014']
sal['2014']=sal['2014'][:,lonind14]
tmp['2014']=tmp['2014'][:,lonind14]
oxy['2014']=oxy['2014'][:,lonind14]


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
oxy['kn221']=zeros((len(prs['kn221']),len(ctd_kn221)))
for ii,dd in enumerate(ctd_kn221):
    profile = cnv.fCNV(dd)
    lat['kn221'][ii]=profile.attributes['LATITUDE']
    lon['kn221'][ii]=profile.attributes['LONGITUDE']
    print(profile.attributes['LONGITUDE'],profile.attributes['datetime'])
    salfnc=interpolate.interp1d(profile['PRES'],profile['PSAL'],kind='linear',fill_value='NaN',bounds_error=False)
    tmpfnc=interpolate.interp1d(profile['PRES'],profile['TEMP'],kind='linear',fill_value='NaN',bounds_error=False)
    oxyfnc=interpolate.interp1d(profile['PRES'],profile['oxygen_ml_L']*44.661,kind='linear',fill_value='NaN',bounds_error=False)
    sal['kn221'][:,ii]=salfnc(prs['kn221'])
    tmp['kn221'][:,ii]=tmpfnc(prs['kn221'])
    oxy['kn221'][:,ii]=oxyfnc(prs['kn221'])

########################################################################################
################################## 2015 section ########################################
########################################################################################
## NO OXYGEN!!

# #For PE 400
# #dtype=[('JD0', 'O'), ('T', '1'), ('S', '2'), ('C', 'O'), ('time0', 'O'), ('datetime', 'O'), ('P', '6'), ('file', 'O'), ('lon', '8'), ('lat', '9'), ('dlon', '10'), ('dlat', '11')])}
# penum=1
#
# dat=io.loadmat(datadir+'Shipboard/CTD_fromPenny/PE400 CTD data/pe400_Ar7E_CTD_section.mat')
# # dat=io.loadmat(datadir+'Shipboard/CTD_fromPenny/PE400 CTD data/PE400_1Hz.mat')['ctd_cal']
#
# datelist_2015=[datetime.datetime(1,1,1)+datetime.timedelta(days=tt-366) for tt in dat['dd'][0]]
# # datelist_2015=[datetime.datetime(1,1,1)+datetime.timedelta(days=tt-366) for tt in dat['time0'][0][0].flatten()]
#
#
# lat[dicvec[penum]]=dat['LAT'].flatten()
# lon[dicvec[penum]]=dat['LON'].flatten()
#
# plot(lon['2015'], datelist_2015,'o')
#
# dat.keys()
# prs[dicvec[penum]]=dat['P'].flatten()
# sal[dicvec[penum]]=dat['S']
# tmp[dicvec[penum]]=dat['T']
# shape(sal['2015'])
# shape(prs['2015'])
# shape(lon['2015'])
#
# plot(dat['POTEMP'],'r');
# plot(dat['T'],'b');
# plot(tmp['2015'],'g');
# axhline(2)
#
#
# contourf(lon['2015'],prs['2015'],tmp['2015'])
# colorbar()

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
ctdlist_16=sort(glob.glob(datadir+'Shipboard/ar07_2016/cnv_files/OSNAPEast_sec6/*.cnv'))

arnum=2
prs['2016']=arange(1,3180,2) # determined in retrospect after inspecting
lat[dicvec[arnum]]=zeros(len(ctdlist_16))
lon[dicvec[arnum]]=zeros(len(ctdlist_16))
sal[dicvec[arnum]]=zeros((len(prs['2016']),len(ctdlist_16)))
tmp[dicvec[arnum]]=zeros((len(prs['2016']),len(ctdlist_16)))
oxy[dicvec[arnum]]=zeros((len(prs['2016']),len(ctdlist_16)))
for ii,cc in enumerate(ctdlist_16):
    profile = cnv.fCNV(cc)
    lat['2016'][ii]=profile.attributes['LATITUDE']
    lon['2016'][ii]=profile.attributes['LONGITUDE']
    print(profile.attributes['LONGITUDE'],profile.attributes['datetime'])
    salfnc=interpolate.interp1d(profile['PRES'],profile['PSAL'],kind='linear',fill_value='NaN',bounds_error=False)
    tmpfnc=interpolate.interp1d(profile['PRES'],profile['TEMP'],kind='linear',fill_value='NaN',bounds_error=False)
    oxyfnc=interpolate.interp1d(profile['PRES'],profile['oxygen_ml_L']*44.661,kind='linear',fill_value='NaN',bounds_error=False)
    sal['2016'][:,ii]=salfnc(prs['2016'])
    tmp['2016'][:,ii]=tmpfnc(prs['2016'])
    oxy['2016'][:,ii]=oxyfnc(prs['2016'])

########################################################################################
################################## Grid all ########################################
########################################################################################
for key in lon:
    sort00=argsort(lon[key])
    lon[key]=lon[key][sort00]
    lat[key]=lat[key][sort00]
    sal[key]=sal[key][:,sort00]
    tmp[key]=tmp[key][:,sort00]
    oxy[key]=oxy[key][:,sort00]


# Find distance from CF1
lon0=CFlon[0]
lat0=CFlat[0]

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

# Grid onto a new grid!

smprs=range(0, 3125, 10)
distgrid=range(0,int(max(dist['2014']))+10,5)
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
dicvec

newdic=[ '2014','kn221','2016']
# newdic=['2005', '2008', '2014','kn221', '2016']


saltot=NaN*ones((smprslen,distlen,3))
tmptot=NaN*ones((smprslen,distlen,3))
oxytot=NaN*ones((smprslen,distlen,3))

for ii,key in enumerate(newdic):
        print(key)
        saltot[:,:,ii]=gridel(sal[key],dist[key],key)
        tmptot[:,:,ii]=gridel(tmp[key],dist[key],key)
        oxytot[:,:,ii]=gridel(oxy[key],dist[key],key)


pdentot=NaN*ones(shape(saltot))
for kk in range(3):
    print(kk)
    for jj in range(distlen):
        salnan=~isnan(saltot[:,jj,kk])
        if sum(salnan)>2:
            SA=gsw.SA_from_SP(saltot[:,jj,kk][salnan],array(smprs)[salnan],[lon0]*sum(salnan),[lat0]*sum(salnan))
            pdentot[salnan,jj,kk]=gsw.sigma0(SA,gsw.CT_from_t(SA,tmptot[:,jj,kk][salnan],array(smprs)[salnan]))

occvec=['2014 (JR302)','2014 (KN221)', '2016']
# occvec=['2005', '2008', '2014 (JR302)','2014 (KN221)', '2016']

#Get into an xarray
grdat=xr.Dataset({'sal': (['prs','dist', 'occ'],  saltot),
                    'tmp': (['prs','dist', 'occ'],  tmptot),
                    'pden': (['prs','dist', 'occ'],  pdentot),
                    'o2': (['prs','dist', 'occ'],  oxytot),},
                coords={'dist': distgrid,
                        'prs': smprs,
                        'occ':occvec})

grdat.to_netcdf(datadir+'Shipboard/gridded/Shipboard_2014x2_2016_wO2.nc','w',format='netCDF4')

grdat['o2'].sel(occ='2016').plot()


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
