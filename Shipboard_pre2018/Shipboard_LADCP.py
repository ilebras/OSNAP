
#### Load and plot LADCP data, eventually make mean

from aux_funcs import *

############################################################################################################
##MIGHT WANT TO TRY USING PENNY'S BATHYMETRY FOR THIS -- USE CAUTION WITH THIS AS THERE IS ANOMALOUS WHITE SPACE BETWEEN EGC AND SLOPE EMERGING
############################################################################################################
#### Load into same vertical grid where it makes sense to

minprs=5
maxprs=3125
prsvec=range(minprs,maxprs,5)
prslen=len(prsvec)

############################################################################################################
################################## 2005, 2006 and 2008 sections ###########################################
###################################################################################################
#[('name', 'O'), ('lat', 'O'), ('lon', 'O'), ('date', 'O'), ('tim', 'O'), ('z', 'O'), ('u', 'O'), ('v', 'O'), ('u_do', 'O'), ('v_do', 'O'), ('u_up', 'O'), ('v_up', 'O'), ('u_corr', 'O'), ('v_corr', 'O')


mat05=sort(glob.glob(datadir+'Shipboard/LADCP/2005_D298/*.mat'))
mat06=sort(glob.glob(datadir+'Shipboard/LADCP/2006_D309/*.mat'))
mat08=sort(glob.glob(datadir+'Shipboard/LADCP/2008_D332/*.mat'))

def load05_06_08(matname,yy):

    distlen00=len(matname)
    date00=[0]*(distlen00)
    lat00=zeros(distlen00)
    lon00=zeros(distlen00)
    u00=zeros((prslen,distlen00))
    v00=zeros((prslen,distlen00))
    for ii,dd in enumerate(matname):

        dat=io.loadmat(dd)['ladcp'][0][0]
        if '06' in yy:
            date00[ii]=dat[3][0][0];
        else:
            date00[ii]=datetime.datetime(dat[3][0][0],dat[3][0][1],dat[3][0][2],dat[3][0][3],dat[3][0][4]);

        lat00[ii]=dat[1][0][0]
        lon00[ii]=dat[2][0][0]

        if '06' in yy:
            ztmp=dat[4].flatten()
        else:
            ztmp=dat[5].flatten()
        utmp=dat[-2].flatten()
        vtmp=dat[-1].flatten()

        ufunc=interpolate.interp1d(hstack((5,ztmp)),
                                hstack((utmp[0],utmp)),
                                kind='linear',fill_value='NaN',bounds_error=False)

        vfunc=interpolate.interp1d(hstack((5,ztmp)),
                                    hstack((vtmp[0],vtmp)),
                                    kind='linear',fill_value='NaN',bounds_error=False)

        u00[:,ii]=ufunc(prsvec)
        v00[:,ii]=vfunc(prsvec)

    lonind00=(lon00>-44)&(lon00<-20)
    v00=v00[:,lonind00]
    u00=u00[:,lonind00]
    lon00=lon00[lonind00]
    lat00=lat00[lonind00]

    if '05' in yy:
        ind5=22
        lat00=lat00[:ind5]
        lon00=lon00[:ind5]
        date00=date00[:ind5]
        v00=v00[:,:ind5]
        u00=u00[:,:ind5]

    elif '08' in yy:
        lon00=hstack((lon00[8:20],lon00[23:]))
        lat00=hstack((lat00[8:20],lat00[23:]))
        date00=hstack((date00[8:20],date00[23:]))
        v00=hstack((v00[:,8:20],v00[:,23:]))
        u00=hstack((u00[:,8:20],u00[:,23:]))

    sort00=argsort(lon00)
    lon00=lon00[sort00]
    lat00=lat00[sort00]
    date00=[date00[ss] for ss in sort00]
    u00=u00[:,sort00]
    v00=v00[:,sort00]


    return lon00,lat00,u00,v00,date00

lon08,lat08,u08,v08,date08=load05_06_08(mat08,'08')
plot(v08);

lon05,lat05,u05,v05,date05=load05_06_08(mat05,'05')

lon06,lat06,u06,v06,date06=load05_06_08(mat06,'06')

plot(v06);

figure()
plot(lon05,lat05,'o',label='2005')
plot(lon06,lat06,'o',label='2006')
plot(lon08,lat08,'o',label='2008')
legend()


############################################################################################################
################################## 2011 section ####################################################
############################################################################################################
##('name', 'O'), ('date', 'O'), ('lat', 'O'), ('lon', 'O'), ('zbot', 'O'), ('ubot', 'O'), ('vbot', 'O'), ('uerrbot', 'O'), ('z', 'O'), ('u', 'O'), ('v', 'O'), ('nvel', 'O'), ('ubar', 'O'), ('vbar', 'O'), ('tim', 'O'), ('shiplon', 'O'), ('shiplat', 'O'), ('xship', 'O'), ('yship', 'O'), ('uship', 'O'), ('vship', 'O'), ('zctd', 'O'), ('wctd', 'O'), ('uctd', 'O'), ('vctd', 'O'), ('xctd', 'O'), ('yctd', 'O'), ('uerr', 'O'), ('range', 'O'), ('ts', 'O'), ('ts_out', 'O'), ('p', 'O'), ('uctderr', 'O'), ('u_do', 'O'), ('v_do', 'O'), ('u_up', 'O'), ('v_up', 'O'), ('ensemble_vel_err', 'O'), ('u_shear_method', 'O'), ('v_shear_method', 'O'), ('ctd_t', 'O'), ('ctd_s', 'O'), ('ctd_N2', 'O'), ('Kz', 'O')

mat11=sort(glob.glob(datadir+'Shipboard/LADCP/2011_fromFemke/*.mat'))

distlen11=len(mat11)
date11=[0]*(distlen11)
lat11=zeros(distlen11)
lon11=zeros(distlen11)
u11=zeros((prslen,distlen11))
v11=zeros((prslen,distlen11))
for ii,dd in enumerate(mat11):

    dat=io.loadmat(dd)['dr'][0][0]

    date11[ii]=datetime.datetime(dat[1][0][0],dat[1][0][1],dat[1][0][2],dat[1][0][3],dat[1][0][4])

    lat11[ii]=dat[2][0][0]
    lon11[ii]=dat[3][0][0]

    ztmp=dat[8].flatten()
    utmp=dat[9].flatten()
    vtmp=dat[10].flatten()

    ufunc=interpolate.interp1d(hstack((5,ztmp)),
                            hstack((utmp[0],utmp)),
                            kind='linear',fill_value='NaN',bounds_error=False)

    vfunc=interpolate.interp1d(hstack((5,ztmp)),
                                hstack((vtmp[0],vtmp)),
                                kind='linear',fill_value='NaN',bounds_error=False)

    u11[:,ii]=ufunc(prsvec)
    v11[:,ii]=vfunc(prsvec)


lonind11=(lon11<-30.5)
v11=v11[:,lonind11]
lon11=lon11[lonind11]
lat11=lat11[lonind11]

sort11=argsort(lon11)
lon11=lon11[sort11]
lat11=lat11[sort11]
u11=u11[:,sort11]
v11=v11[:,sort11]

#######################################################################################################
################################## 2014 KN221 section #################################################
######################################################################################################
# dtype=[('name', 'O'), ('date', 'O'), ('lat', 'O'), ('lon', 'O'), ('zbot', 'O'), ('ubot', 'O'), ('vbot', 'O'), ('uerrbot', 'O'), ('z_sadcp', 'O'), ('u_sadcp', 'O'), ('v_sadcp', 'O'), ('uerr_sadcp', 'O'),
#('z', 'O'), ('u', 'O'), ('v', 'O'), ('nvel', 'O'), ('ubar', 'O'), ('vbar', 'O'), ('tim', 'O'), ('tim_hour', 'O'), ('zctd', 'O'), ('wctd', 'O'), ('uctd', 'O'), ('vctd', 'O'), ('xctd', 'O'), ('yctd', 'O'),
#('uerr', 'O'), ('range', 'O'), ('range_do', 'O'), ('range_up', 'O'), ('ts', 'O'), ('ts_out', 'O'), ('p', 'O'), ('uctderr', 'O'), ('u_do', 'O'), ('v_do', 'O'), ('u_up', 'O'), ('v_up', 'O'),
#('ensemble_vel_err', 'O'), ('u_shear_method', 'O'), ('v_shear_method', 'O'), ('w_shear_method', 'O')])

matkn221=sort(glob.glob(datadir+'Shipboard/LADCP/2014_KN221_dtorres/*.mat'))

knlen=len(matkn221)
latkn=zeros(knlen)
lonkn=zeros(knlen)
ukn=zeros((prslen,knlen))
vkn=zeros((prslen,knlen))
for ii,dd in enumerate(matkn221):

    dat=io.loadmat(dd)['dr'][0][0]

    latkn[ii]=float(dat[2])
    lonkn[ii]=float(dat[3])


    ztmp=dat[12].flatten()
    utmp=dat[13].flatten()
    vtmp=dat[14].flatten()

    ufunc=interpolate.interp1d(hstack((10,ztmp)),
                            hstack((utmp[0],utmp)),
                            kind='linear',fill_value='NaN',bounds_error=False)

    vfunc=interpolate.interp1d(hstack((10,ztmp)),
                                hstack((vtmp[0],vtmp)),
                                kind='linear',fill_value='NaN',bounds_error=False)

    ukn[:,ii]=ufunc(prsvec)
    vkn[:,ii]=vfunc(prsvec)


plot(lonkn,latkn,'o')

############################################################################################################
#################################### 2015 section ##########################################################
############################################################################################################


#################################### Variable names for LADCP_64PE400_deSteur_2015 #######################################
##       dtype=[('name', 'O'), ('date', 'O'), ('lat', 'O'), ('lon', 'O'), ('zbot', 'O'), ('ubot', 'O'), ('vbot', 'O'), ('uerrbot', 'O'), ('z', 'O'), ('u', 'O'), ('v', 'O'), ('nvel', 'O'), ('ubar', 'O'), ('vbar', 'O'), ('tim', 'O'), ('tim_hour', 'O'), ('shiplon', 'O'), ('shiplat', 'O'), ('xship', 'O'), ('yship', 'O'), ('uship', 'O'), ('vship', 'O'), ('zctd', 'O'), ('wctd', 'O'), ('uctd', 'O'), ('vctd', 'O'), ('xctd', 'O'), ('yctd', 'O'), ('uerr', 'O'), ('range', 'O'), ('range_do', 'O'), ('range_up', 'O'), ('ts', 'O'), ('ts_out', 'O'), ('p', 'O'), ('uctderr', 'O'), ('u_do', 'O'), ('v_do', 'O'), ('u_up', 'O'), ('v_up', 'O'), ('ensemble_vel_err', 'O'), ('u_shear_method', 'O'), ('v_shear_method', 'O'), ('w_shear_method', 'O'), ('ctd_t', 'O'), ('ctd_s', 'O'), ('ctd_ss', 'O'), ('u_detide', 'O'), ('v_detide', 'O'), ('utide', 'O'), ('vtide', 'O')])}

mat15=sort(glob.glob(datadir+'Shipboard/LADCP/2015_64PE400_deSteur/*.mat'))

distlen15=36
date15=[0]*(distlen15)
lat15=zeros(distlen15)
lon15=zeros(distlen15)
u15=zeros((prslen,distlen15))
v15=zeros((prslen,distlen15))
for ii,dd in enumerate(hstack((mat15[:24],mat15[27:39]))):

    dat=io.loadmat(dd)['dr'][0][0]

    date15[ii]=datetime.datetime(dat[1][0][0],dat[1][0][1],dat[1][0][2],dat[1][0][3],dat[1][0][4])

    lat15[ii]=dat[2][0][0]
    lon15[ii]=dat[3][0][0]


    ztmp=dat[8].flatten()
    utmp=dat[-4].flatten()
    vtmp=dat[-3].flatten()

    ufunc=interpolate.interp1d(hstack((5,ztmp)),
                            hstack((utmp[0],utmp)),
                            kind='linear',fill_value='NaN',bounds_error=False)

    vfunc=interpolate.interp1d(hstack((5,ztmp)),
                                hstack((vtmp[0],vtmp)),
                                kind='linear',fill_value='NaN',bounds_error=False)

    u15[:,ii]=ufunc(prsvec)
    v15[:,ii]=vfunc(prsvec)



sort15=argsort(lon15)
lon15=lon15[sort15]
lat15=lat15[sort15]
u15=u15[:,sort15]
v15=v15[:,sort15]




############################################################################################################
################################## 2014 & 2016 sections from Penny ########################################
############################################################################################################

mat14=sort(glob.glob(datadir+'Shipboard/LADCP/2014_2016_gridded_fromPenny/2014*'))[0]
dat14=io.loadmat(mat14)

mat16=sort(glob.glob(datadir+'Shipboard/LADCP/2014_2016_gridded_fromPenny/2016*'))[0]
dat16=io.loadmat(mat16)

u14=dat14['corr_lad_u']
v14=dat14['corr_lad_v']
z14=dat14['lad_z']
lon14=dat14['lon'][0]
lat14=dat14['lat'][0]

u16=dat16['corr_lad_u']
v16=dat16['corr_lad_v']
z16=dat16['lad_z']
lon16=dat16['lon'][0]
lat16=dat16['lat'][0]



# For now, focus on data between 44W and 30.5W, and discard anomalous values

lonind14=(lon14<-30.5) & (lon14>-44)
u14=u14[:,lonind14]
v14=v14[:,lonind14]
z14=z14[:,lonind14]
lon14=lon14[lonind14]
lat14=lat14[lonind14]

lonind16=(lon16<-30.5) & (lon16>-44)
u16=u16[:,lonind16]
v16=v16[:,lonind16]
z16=z16[:,lonind16]
lon16=lon16[lonind16]
lat16=lat16[lonind16]

z14vec=sort(unique(z14))
z14vec=z14vec[~isnan(z14vec)]

z16vec=sort(unique(z16))
z16vec=z16vec[~isnan(z16vec)]

u14=u14[:len(z14vec),:]
u16=u16[:len(z16vec),:]
v14=v14[:len(z14vec),:]
v16=v16[:len(z16vec),:]


#################################################################################################
################################## Plot all ####################################################
################################################################################################

figure()
plot(lon16,lat16,'o',label='2016')
plot(lon15,lat15,'o',label='2015')
plot(lon14,lat14,'o',label='2014 (jr302)')
plot(lonkn,latkn,'o',label='2014 (kn221)')
plot(lon11,lat11,'o',label='2011')
legend()
ylabel('Latitude ($^\circ$ N)')
xlabel('Longitude ($^\circ$ W)')
savefig('../figures/LADCP/measpos.png',bbox_inches='tight')

figure()
plot(lon14,lat14,'o',label='2014 (jr302)')
plot(lon14,lat14,'o',label='2014 (kn221)')
plot(lon15,lat15,'o',label='2015')
plot(lon16,lat16,'o',label='2016')
plot(lon05,lat05,'o',label='2005')
plot(lon06,lat06,'o',label='2006')
plot(lon08,lat08,'o',label='2008')
legend()
ylabel('Latitude ($^\circ$ N)')
xlabel('Longitude ($^\circ$ W)')
savefig('../figures/Shipboard_forPenny/measpos_all.png',bbox_inches='tight')




# def addcontplot(savelab):
#     colorbar(ticks=arange(-0.4,0.5,0.2),label='[m/s]')
#     title('Meridional velocity, '+savelab,fontsize=18)
#     xlabel('longitude [$^\circ$ W]')
#     ylabel('pressure [db]')
#     ylim([3000,0])
#     xlim([-43,-30.5])
#     savefig('../figures/Shipboard_forPenny/v'+savelab+'_full.png')
#     xlim([-43,-38])
#     savefig('../figures/Shipboard_forPenny/v'+savelab+'_zoom.png')
#
#
# figure(figsize=(12,6))
# contourf(lon11,prsvec,v11,21,cmap=cm.RdBu_r,vmin=-0.5,vmax=0.5)
# addcontplot('2011')
#
# figure(figsize=(12,6))
# contourf(lon15,prsvec,v15,21,cmap=cm.RdBu_r,vmin=-0.5,vmax=0.5)
# addcontplot('2015')
#
#
# figure(figsize=(12,6))
# contourf(lon14,z14vec,v14[:len(z14vec),:],21,cmap=cm.RdBu_r,vmin=-0.5,vmax=0.5)
# addcontplot('2014')
#
# figure(figsize=(12,6))
# contourf(lon16,z16vec,v16[:len(z16vec),:],21,cmap=cm.RdBu_r,vmin=-0.5,vmax=0.5)
# addcontplot('2016')


#################################################################################################
################################## Grid all ####################################################
################################################################################################

# Find distance from Jo's zero point
################################################################################################
lon0=-43.0884
lat0=60.0911




def getdist(lonex,latex):
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([lat0,latex[ii]],[lon0,lonex[ii]])[0][0]
    return distout

dist05=getdist(lon05,lat05)
dist06=getdist(lon06,lat06)
dist08=getdist(lon08,lat08)

dist14=getdist(lon14,lat14)
dist11=getdist(lon11,lat11)
dist15=getdist(lon15,lat15)
dist16=getdist(lon16,lat16)

dist06

distkn=getdist(lonkn,latkn)

# # Use true, old across track definition here!
xdist=sw.dist([CFlat[-2],CFlat[-2]],[CFlon[0],CFlon[-2]])[0][0]
ydist=sw.dist([CFlat[0],CFlat[-2]],[CFlon[-2],CFlon[-2]])[0][0]
theta=abs(arctan((xdist)/(ydist)))
#
# Angle for 2005,2006
xdist2=sw.dist([lat05[5],lat05[5]],[lon05[0],lon05[10]])[0][0]
ydist2=sw.dist([lat05[0],lat05[10]],[lon05[5],lon05[5]])[0][0]
xdist2
ydist2
theta2=abs(arctan((xdist2)/(ydist2)))

90-theta*180/pi
pi/2

figure()
plot(lon05[dist05<300],lat05[dist05<300],'o',label='2005')
plot(lon06[dist06<300],lat06[dist06<300],'o',label='2006')
plot(lon08[dist08<300],lat08[dist08<300],'o',label='2008')
plot(lon14[dist14<300],lat14[dist14<300],'o',label='2014 (jr302)')
plot(lonkn[distkn<300],latkn[distkn<300],'o',label='2014 (kn221)')
plot(lon15[dist15<300],lat15[dist15<300],'o',label='2015')
plot(lon16[dist16<300],lat16[dist16<300],'o',label='2016')
plot(CFlon[2:],CFlat[2:],'o',label='CF3-M1')
legend()
ylabel('Latitude ($^\circ$ N)')
xlabel('Longitude ($^\circ$ W)')
savefig('../figures/Shipboard_forPenny/measpos_dwreg.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/measpos_dwreg.pdf',bbox_inches='tight')

### Get across track velocity:
################################################################################################

def getacrosstrack(u00,v00,theta00):
    uac00=u00*cos(theta00)+v00*sin(theta00)

    return uac00


##For now, get mer and zon vels to share w Jo
theta=pi/2

uac05=getacrosstrack(u05,v05,theta)
uac06=getacrosstrack(u06,v06,theta)
uac08=getacrosstrack(u08,v08,theta)
uac14=getacrosstrack(u14,v14,theta)
uackn=getacrosstrack(ukn,vkn,theta)
uac15=getacrosstrack(u15,v15,theta)
uac16=getacrosstrack(u16,v16,theta)

################################################################################################
### Smooth/interpolate each dataset onto same grid --> make an xarray!
################################################################################################

distgrid=range(0,int(max(dist14))+10,10)
distlen=len(distgrid)
smprslen=int(prslen/2)
smprs=prsvec[1::2]
def gridv(vex,distex):
    # 1. Smooth each dataset into 10db bins in the vertical (use above and below measurement for each)
    vsm=(vex[::2,:][:-1,:]+vex[1::2,:][:-1,:]+vex[2::2])/3

    # 2. Extend each profile downwards
    vext=vsm.copy()
    # for jj in range(shape(vsm)[1]):
    #     vext[isnan(vsm[:,jj]),jj]=vsm[~isnan(vsm[:,jj]),jj][-1]

    # 3. Interpolate horizontally onto a uniform 10km dist grid
    vhorz=zeros((smprslen,distlen))
    for ii in range(shape(vsm)[0]):
            distfunc=interpolate.interp1d(distex,vext[ii,:],kind='linear',fill_value='NaN',bounds_error=False)
            vhorz[ii,:]=distfunc(distgrid)

    return vhorz

v14grid=gridv(uac14,dist14)
v15grid=gridv(uac15,dist15)
# v11grid=gridv(uac11,dist11)
v16grid=gridv(uac16,dist16)
v05grid=gridv(uac05,dist05)
v06grid=gridv(uac06,dist06)
v08grid=gridv(uac08,dist08)

vkngrid=gridv(uackn,distkn)

vtot=zeros((smprslen,distlen,7))
vtot[:,:,0]=v05grid
vtot[:,:,1]=v06grid
vtot[:,:,2]=v08grid
# vtot[:,:,3]=v11grid
vtot[:,:,3]=v14grid
vtot[:,:,4]=vkngrid
vtot[:,:,5]=v15grid
vtot[:,:,6]=v16grid

u14grid=gridv(u14,dist14)
u15grid=gridv(u15,dist15)
u16grid=gridv(u16,dist16)
u05grid=gridv(u05,dist05)
u06grid=gridv(u06,dist06)
u08grid=gridv(u08,dist08)
ukngrid=gridv(ukn,distkn)

utot=zeros((smprslen,distlen,7))
utot[:,:,0]=u05grid
utot[:,:,1]=u05grid
utot[:,:,2]=u08grid
utot[:,:,3]=u14grid
utot[:,:,4]=ukngrid
utot[:,:,5]=u15grid
utot[:,:,6]=u16grid

## Make into an xarray
grdat=xr.Dataset({'zonal velocity, u [m/s]': (['pressure [db]','distance [km]', 'occupation'],  utot),
                    'meridional velocity, v [m/s]': (['pressure [db]','distance [km]', 'occupation'],  vtot)},
                coords={'distance [km]': distgrid,
                        'pressure [db]': smprs,
                        'occupation': ['2005','2006','2008','2014 (JR302)','2014 (KN221)','2015','2016']})

## Make into an xarray
# grdat=xr.Dataset({'across track velocity [m/s]': (['pressure [db]','distance [km]', 'occupation'],  vtot)},
#                 coords={'distance [km]': distgrid,
#                         'pressure [db]': smprs,
#                         'occupation': ['2005','2008','2014 (JR302)','2014 (KN221)','2015','2016']})



d1=0
d2=30

g=grdat['zonal velocity, u [m/s]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-0.4,vmax=0.4,cmap=cm.RdBu_r,figsize=(15,12))
ylim([3e3,0]);
xlim([0,250])
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    # ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    ax.contour(distgrid[d1:d2],smprs,grdat['zonal velocity, u [m/s]'][:,d1:d2,ii],arange(-0.4,0.5,0.1),colors='k')
savefig('../figures/Shipboard_forPenny/LADCP_uzonal_panels.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/LADCP_uzonal_panels.pdf',bbox_inches='tight')

g=grdat['meridional velocity, v [m/s]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-0.4,vmax=0.4,cmap=cm.RdBu_r,figsize=(15,12))
ylim([3e3,0]);
xlim([0,250])
for ii, ax in enumerate(g.axes.flat):
    ax.set_title(grdat.occupation.values[ii],fontsize=20)
    # ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
    ax.contour(distgrid[d1:d2],smprs,grdat['meridional velocity, v [m/s]'][:,d1:d2,ii],arange(-0.4,0.5,0.1),colors='k')
savefig('../figures/Shipboard_forPenny/LADCP_umer_panels.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/LADCP_umer_panels.pdf',bbox_inches='tight')

# g=grdat['across track velocity [m/s]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-0.4,vmax=0.4,cmap=cm.RdBu_r,figsize=(15,12))
# ylim([3e3,0]);
# xlim([0,250])
# for ii, ax in enumerate(g.axes.flat):
#     ax.set_title(grdat.occupation.values[ii],fontsize=20)
#     ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
#     ax.contour(distgrid[d1:d2],smprs,grdat['across track velocity [m/s]'][:,d1:d2,ii],arange(-0.4,0.5,0.1),colors='k')
# savefig('../figures/Shipboard_forPenny/LADCP_uacross_panels.png',bbox_inches='tight')
# savefig('../figures/Shipboard_forPenny/LADCP_uacross_panels.pdf',bbox_inches='tight')
#
# d1=0
# d2=10
# g=grdat['across track velocity [m/s]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-0.6,vmax=0.6,cmap=cm.RdBu_r,figsize=(15,12))
# ylim([400,0]);
# for ii, ax in enumerate(g.axes.flat):
#     ax.set_title(grdat.occupation.values[ii],fontsize=20)
#     ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
#     ax.contour(distgrid[d1:d2],smprs,grdat['across track velocity [m/s]'][:,d1:d2,ii],arange(-0.4,0.5,0.1),colors='k')
# savefig('../figures/Shipboard_forPenny/firstsix_extdpth_uacross_shelfslope.png',bbox_inches='tight')
#
# d1=2
# d2=-1
#
# g=grdat['across track velocity [m/s]'][:,d1:d2,:].plot(x='distance [km]', y='pressure [db]', col='occupation', col_wrap=2,vmin=-0.6,vmax=0.6,cmap=cm.RdBu_r,figsize=(15,12))
# ylim([3e3,0])
# for ii, ax in enumerate(g.axes.flat):
#     ax.set_title(grdat.occupation.values[ii],fontsize=20)
#     ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
#     ax.contour(distgrid[d1:d2],smprs,grdat['across track velocity [m/s]'][:,d1:d2,ii],arange(-0.4,0.5,0.1),colors='k')
# savefig('../figures/Shipboard_forPenny/firstsix_extdpth_uacross_fullbasin.png',bbox_inches='tight')



d1=0
d2=30

figure(figsize=(10,6))
contourf(distgrid[d1:d2],smprs,grdat['meridional velocity, v [m/s]'][:,d1:d2,:].mean(dim='occupation'),21,vmin=-0.3,vmax=0.3,cmap=cm.RdBu_r,extend='both')
colorbar(ticks=arange(-0.4,0.5,0.1),label='[m/s]')
contour(distgrid[d1:d2],smprs,grdat['meridional velocity, v [m/s]'][:,d1:d2,:].mean(dim='occupation'),arange(-0.4,0.5,0.1),colors='k')
title('Mean meridional velocity from LADCP',fontsize=20)
xlabel('distance [km]')
ylabel('pressure [db]')
ylim([3e3,0])
xlim([0,250])
savefig('../figures/Shipboard_forPenny/LADCP_umer_mean.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/LADCP_umer_mean.pdf',bbox_inches='tight')

figure(figsize=(10,6))
contourf(distgrid[d1:d2],smprs,grdat['zonal velocity, u [m/s]'][:,d1:d2,:].mean(dim='occupation'),21,vmin=-0.3,vmax=0.3,cmap=cm.RdBu_r,extend='both')
colorbar(ticks=arange(-0.4,0.5,0.1),label='[m/s]')
contour(distgrid[d1:d2],smprs,grdat['zonal velocity, u [m/s]'][:,d1:d2,:].mean(dim='occupation'),arange(-0.4,0.5,0.1),colors='k')
title('Mean zonal velocity from LADCP',fontsize=20)
xlabel('distance [km]')
ylabel('pressure [db]')
ylim([3e3,0])
xlim([0,250])
savefig('../figures/Shipboard_forPenny/LADCP_uzon_mean.png',bbox_inches='tight')
savefig('../figures/Shipboard_forPenny/LADCP_uzon_mean.pdf',bbox_inches='tight')

grdat['distance [km]']

grdat.occupation.values

CTD=pickle.load(open('../pickles/Shipboard/CTD_xarray_1810_forJo.pickle','rb'))



datdic={}
for oo in grdat.occupation.values:
    if 'KN221' in oo:
        datdic['u_2014_KN221']=grdat['zonal velocity, u [m/s]'].sel(occupation=oo).values
        datdic['v_2014_KN221']=grdat['meridional velocity, v [m/s]'].sel(occupation=oo).values
        datdic['sig_2014_KN221']=CTD.density.sel(occupation=oo).values
    elif 'JR302' in oo:
        datdic['u_2014_JR302']=grdat['zonal velocity, u [m/s]'].sel(occupation=oo).values
        datdic['v_2014_JR302']=grdat['meridional velocity, v [m/s]'].sel(occupation=oo).values
        datdic['sig_2014_JR302']=CTD.density.sel(occupation=oo).values
    else:
        datdic['u_'+oo]=grdat['zonal velocity, u [m/s]'].sel(occupation=oo).values
        datdic['v_'+oo]=grdat['meridional velocity, v [m/s]'].sel(occupation=oo).values
        datdic['sig_'+oo]=CTD.density.sel(occupation=oo).values
datdic['pressure']=grdat['pressure [db]'].values
datdic['distance']=grdat['distance [km]'].values

datdic.keys()

io.savemat(open('../data/Shipboard/forJo/LADCP_gridded_1810.mat','wb'),datdic)


# datdic={}
# for oo in grdat.occupation.values:
#     datdic[oo]=grdat['zonal velocity, u [m/s]'].sel(occupation=oo).values
# datdic['pressure']=grdat['pressure [db]'].values
# datdic['distance']=grdat['distance [km]'].values
#
# io.savemat(open('../data/Shipboard/forJo/LADCP_u_gridded_180427.mat','wb'),datdic)
#
#
# datdic={}
#Ou for oo in grdat.occupation.values:
#     datdic[oo]=grdat['meridional velocity, v [m/s]'].sel(occupation=oo).values
# datdic['pressure']=grdat['pressure [db]'].values
# datdic['distance']=grdat['distance [km]'].values
#
# io.savemat(open('../data/Shipboard/forJo/LADCP_v_gridded_180427.mat','wb'),datdic)
#
#
# datdic={}
# for oo in grdat.occupation.values:
#     datdic[oo]=CTD.density.sel(occupation=oo).values
# datdic['pressure']=grdat['pressure [db]'].values
# datdic['distance']=grdat['distance [km]'].values
#
# io.savemat(open('../data/Shipboard/forJo/LADCP_sig0_gridded_180427.mat','wb'),datdic)

xxxxxxxxxxxxxxxxxxxxxxx
# figure(figsize=(10,6))
# contourf(distgrid[d1:d2],smprs,grdat['across track velocity [m/s]'][:,d1:d2,:].mean(dim='occupation'),21,vmin=-0.3,vmax=0.3,cmap=cm.RdBu_r,extend='both')
# colorbar(ticks=arange(-0.4,0.5,0.1),label='[m/s]')
# contour(distgrid[d1:d2],smprs,grdat['across track velocity [m/s]'][:,d1:d2,:].mean(dim='occupation'),arange(-0.4,0.5,0.1),colors='k')
# # fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
# title('Mean across track velocity from LADCP',fontsize=20)
# ylim([1e3,0])
# xlim([0,100])
# savefig('../figures/Shipboard_forPenny/firstsix_mean_extdpth_uacross_shelfslope.png',bbox_inches='tight')
# xlabel('distance [km]')
# ylabel('pressure [db]')
# ylim([3e3,0])
# xlim([0,250])
# savefig('../figures/Shipboard_forPenny/LADCP_uacross_mean.png',bbox_inches='tight')
# savefig('../figures/Shipboard_forPenny/LADCP_uacross_mean.pdf',bbox_inches='tight')

########################################################################################################
## Choose a distance range and plot vertical profiles of across-track velocity for each year (as well as mean)
# In mean -- looks like 80-200 gets everything
# Individual sections -- that should def do it
########################################################################################################
# grdat
# grdat['pressure [db]'].values
#
# distgrid[12:18]
# colvec=['b','orange','purple','lightgreen','red','cyan']
#
# for rr in range(6):
#     figure(1000)
#     plot(nanmean(vtot[:,12:18,rr],axis=1),smprs,color=colvec[rr],label=(grdat['occupation'].values[rr]));
#     figure()
#     plot(vtot[:,12:18,rr],smprs,color=colvec[rr]);
#     title(grdat['occupation'].values[rr])
#     ylim([2600,1600])
#     axvline(0,color='k')
#     axvline(-0.1,color='k')
#     xlim([-0.3,0.1])
#     xlabel('across track velocity [m/s]')
#     ylabel('pressure [db]')
#     savefig('../figures/Shipboard_forPenny/LADCP_profiles/LADCP_'+grdat['occupation'].values[rr]+'_120-170.pdf',bbox_inches='tight')
# figure(1000)
# ylim([3e3,0])
# legend(loc=(1.05,0.5))
# axvline(0,color='k')
# axvline(-0.1,color='k')
# xlabel('across track velocity [m/s]')
# ylabel('pressure [db]')
# title('Average profile between 120km and 170km')
# savefig('../figures/Shipboard_forPenny/LADCP_profiles/LADCP_meanprofcomp_120-170.pdf',bbox_inches='tight')
#
#
#
# pickle.dump(grdat,open('../pickles/Shipboard/LADCP_xarray.pickle','wb'))
