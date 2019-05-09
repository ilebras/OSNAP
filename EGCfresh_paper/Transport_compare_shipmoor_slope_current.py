################################################################################
################################################################################
##### Compare sections and transports between shipboard and moored measurements
################################################################################
################################################################################

from aux_funcs import *

#################################################################################
#### Load Mooring dailya and set up transport calc framework
#################################################################################

newgrid=pickle.load(open('../../pickles/xarray/CF_xarray_gridplot_notid_1807tide.pickle','rb'))

savename='tide'
#
# bathf=interpolate.interp1d(bathdist,bathbath)
# bathonmygrid=bathf(newgrid.distance)

daily=newgrid.copy()
# for vv in newgrid:
#     if vv[0]!='d':
#         print(vv)
#         for dd,adist in enumerate(daily.distance):
#             daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=bathonmygrid[dd])
#
# bathdist=hstack((-20,bathdist))
# bathbath=hstack((bathbath[0],bathbath))

# some experimentation with current position
# note, further snippets that I used to generate plots in transport_errors doc are at the end of this
# mindxr=daily.date.copy()
# mindxr.values=daily.distance[minind]
#
# maxind=zeros(len(daily.date))
# for tt,na in enumerate(daily.date):
#     maxind[tt]=where(min(daily['across track velocity'][1:12,0,tt])==daily['across track velocity'][1:12,0,tt])[0][0]+1
#
# maxind=[int(mm) for mm in maxind]
# maxind[maxind==11]=1
# maxdxr=daily.date.copy()
# maxdxr.values=daily.distance[maxind]

#
# mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
# mid_dist[1]=1.5 ## being strict -- only on one side of mooring
# mid_dist[-1]=2.5 ## being strict -- only on one side of mooring
# middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
# depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))


mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))

mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))



daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3

sref=34.9

egic={}
egic['trans']=daily.where(daily['potential density']<27.8).xport[6:,:,:].sum('depth').sum('distance')
egic['fresh no sum']=(daily.where(daily['potential density']<27.8)['xport'][6:,:-1,:]*1e3*(daily.where(daily['potential density']<27.8).salinity[6:,:-1,:]-sref)/sref)
egic['fresh']=(daily.where(daily['potential density']<27.8)['xport'][6:,:-1,:]*1e3*(daily.where(daily['potential density']<27.8).salinity[6:,:-1,:]-sref)/sref).sum('depth').sum('distance')

egic['trans'].std()


################################################################################
##### Load LADCP and calc trans
################################################################################
LADCP=pickle.load(open('../../pickles/Shipboard/LADCP_xarray.pickle','rb'))

CTD=pickle.load(open('../../pickles/Shipboard/CTD_xarray_1803.pickle','rb'))

grdat=xr.merge([LADCP,CTD])

grdat

#Going from 10km to 100km in basic accordance with moorings
d1=1
d2=11
mid_dist_ladcp=hstack((diff(grdat['distance [km]'][d1:d2])[0],(diff(grdat['distance [km]'][d1:d2][:-1])+diff(grdat['distance [km]'][d1:d2])[1:])/2,diff(grdat['distance [km]'][d1:d2])[-1]))
middistmat_ladcp=transpose((tile(mid_dist_ladcp,[len(grdat['pressure [db]'])-1,len(grdat['occupation']),1])),(0,2,1))
depthdiffmat_ladcp=transpose((tile(diff(grdat['pressure [db]']),[len(grdat['distance [km]'][d1:d2]),len(grdat['occupation']),1])),(2,0,1))
mid_dist_ladcp
trans_ladcp=(grdat.where(grdat['density']<27.8)['across track velocity [m/s]'][:-1,d1:d2,:]*depthdiffmat_ladcp*middistmat_ladcp/1e3).sum('pressure [db]').sum('distance [km]')
trans_ladcp_noden=(grdat['across track velocity [m/s]'][:-1,d1:d2,:]*depthdiffmat_ladcp*middistmat_ladcp/1e3).sum('pressure [db]').sum('distance [km]')
fresh_ladcp_nosum=(grdat.where(grdat['density']<27.8)['across track velocity [m/s]'][:-1,d1:d2,:]*(grdat.where(grdat['density']<27.8)['salinity'][:-1,d1:d2,:]-sref)/sref*depthdiffmat_ladcp*middistmat_ladcp)
fresh_ladcp=(grdat.where(grdat['density']<27.8)['across track velocity [m/s]'][:-1,d1:d2,:]*(grdat.where(grdat['density']<27.8)['salinity'][:-1,d1:d2,:]-sref)/sref*depthdiffmat_ladcp*middistmat_ladcp).sum('pressure [db]').sum('distance [km]')

fresh_ladcp_nosum.sel(occupation='2014 (KN221)').plot()

fresh_ladcp.sel(occupation='2014 (KN221)').values

egic['fresh'][0].values
mean(egic['fresh'][:30]).values


fresh_ladcp.sel(occupation='2016')
egic['fresh'][-1].values
mean(egic['fresh'][-30:]).values



daily

daily.where(daily['potential density']<27.8).salinity[6,:,:30].mean(axis=1)

def compsalfreshprof(occ,datechoose):
        figure(figsize=(12,4))
        subplot(121)
        plot(daily.where(daily['potential density']<27.8).salinity[6,:,:].isel(date=datechoose).T,daily.depth,'k',label='Moorings, daily: '+str('{:5.2f}'.format(egic['fresh'][datechoose].values))+' mSv')
        plot(daily.where(daily['potential density']<27.8).salinity[6:,:,:].isel(date=datechoose).T,daily.depth,'k',label='')
        if datechoose==0:
            plot(daily.where(daily['potential density']<27.8).salinity[6,:,:30].mean(axis=1).T,daily.depth,'b',label='Moorings, monthly: '+str('{:5.2f}'.format(mean(egic['fresh'][:30]).values))+' mSv')
            plot(daily.where(daily['potential density']<27.8).salinity[6:,:,:30].mean(axis=2).T,daily.depth,'b',label='')
        if datechoose==-1:
            plot(daily.where(daily['potential density']<27.8).salinity[6,:,-30:].mean(axis=1).T,daily.depth,'b',label='Moorings, monthly: '+str('{:5.2f}'.format(mean(egic['fresh'][-30:]).values))+' mSv')
            plot(daily.where(daily['potential density']<27.8).salinity[6:,:,-30:].mean(axis=2).T,daily.depth,'b',label='')
        plot(grdat.where(grdat['density']<27.8)['salinity'][:,d1,:].sel(occupation=occ),grdat['pressure [db]'],color='orange',label='Shipboard: '+str('{:5.2f}'.format(fresh_ladcp.sel(occupation=occ).values))+' mSv')
        plot(grdat.where(grdat['density']<27.8)['salinity'][:,d1:d2,:].sel(occupation=occ),grdat['pressure [db]'],color='orange',label='')
        ylim([1600,0])
        axvline(34.9,color='r',linewidth=2)
        title('profiles of salinity')
        legend()
        subplot(122)
        plot(egic['fresh no sum'].isel(date=datechoose).T,daily.depth[:-1],'k',label='')
        if datechoose==0:
            plot(egic['fresh no sum'][:,:,:30].mean(axis=2).T,daily.depth[:-1],'b',label='')
        if datechoose==-1:
            plot(egic['fresh no sum'][:,:,-30:].mean(axis=2).T,daily.depth[:-1],'b',label='')
        plot(fresh_ladcp_nosum.sel(occupation=occ)/5,grdat['pressure [db]'][:-1],color='orange',label='')
        ylim([1600,0])
        gca().set_yticklabels('')
        title('profiles of freshwater flux (ref to 34.9)')
        suptitle(occ[:4])
        savefig('../../figures/compare_shipmoor/salfreshprofcomp_slope'+occ[:4]+savename+'.pdf')

compsalfreshprof('2014 (KN221)',0)
compsalfreshprof('2016',-1)

grdat

################################################################################
##### Salinity comparison paper figure
################################################################################


def salcomp_paper():
        f=figure(figsize=(8,5))
        for ii in range(1,3):
            subplot(1,2,ii)
            plot(daily.salinity[6:,:,:].isel(date=0).mean(axis=0).T,daily.depth,'k',label='Moorings, first day: '+str('{:5.2f}'.format(egic['fresh'][0].values))+' mSv')
            plot(daily.salinity[6:,:,:].isel(date=-1).mean(axis=0).T,daily.depth,'k--',label='Moorings, last day: '+str('{:5.2f}'.format(egic['fresh'][-1].values))+' mSv')
            plot(daily.salinity[6:,:,:30].mean(axis=0).mean(axis=1).T,daily.depth,'b',label='Moorings, first month: '+str('{:5.2f}'.format(mean(egic['fresh'][:30]).values))+' mSv')
            plot(daily.salinity[6:,:,-30:].mean(axis=0).mean(axis=1).T,daily.depth,'b--',label='Moorings, last month: '+str('{:5.2f}'.format(mean(egic['fresh'][-30:]).values))+' mSv')
            plot(grdat.salinity[:,d1:d2,:].sel(occupation='2014 (KN221)').mean(axis=1),grdat['pressure [db]'],'r',label='Shipboard 2014: '+str('{:5.2f}'.format(fresh_ladcp.sel(occupation='2014 (KN221)').values))+' mSv')
            plot(grdat.salinity[:,d1:d2,:].sel(occupation='2016').mean(axis=1),grdat['pressure [db]'],'r--',label='Shipboard 2016: '+str('{:5.2f}'.format(fresh_ladcp.sel(occupation='2016').values))+' mSv')
            ylim([1600,0])
        subplot(121)
        xlim([34,34.9])
        legend()
        ylabel('depth [m]',fontsize=14)
        subplot(122)
        gca().set_xticks(arange(34.9,35.1,0.05))
        gca().set_yticklabels('')
        xlim([34.9,35.05])
        tight_layout()
        f.text(0.5, -0.03, 'mean slope current salinity', ha='center',fontsize=14)
        savefig('../../figures/paperfigs/salcomp_slope.pdf')
salcomp_paper()


################################################################################
##### Load LADCP pre-smoothing
################################################################################

matkn221=sort(glob.glob(datadir+'Shipboard/LADCP/2014_KN221_dtorres/*.mat'))

prsvec=range(2,3120,2)
prslen=len(prsvec)

knlen=len(matkn221)
lat14_ladcp=zeros(knlen)
lon14_ladcp=zeros(knlen)
ukn=zeros((prslen,knlen))
vkn=zeros((prslen,knlen))
for ii,dd in enumerate(matkn221):

    dat=io.loadmat(dd)['dr'][0][0]

    lat14_ladcp[ii]=float(dat[2])
    lon14_ladcp[ii]=float(dat[3])

    ztmp=dat[12].flatten()
    utmp=dat[13].flatten()
    vtmp=dat[14].flatten()

    ufunc=interpolate.interp1d(hstack((2,ztmp)),
                            hstack((utmp[0],utmp)),
                            kind='linear',fill_value='NaN',bounds_error=False)

    vfunc=interpolate.interp1d(hstack((2,ztmp)),
                                hstack((vtmp[0],vtmp)),
                                kind='linear',fill_value='NaN',bounds_error=False)

    ukn[:,ii]=ufunc(prsvec)
    vkn[:,ii]=vfunc(prsvec)

def getdist(lonex,latex):
    # get distance from longitude of CF1, but projected onto the nearest point in latitude for line
    # this dist is only really useful for the lines close to CF moorings at this point
    distout=zeros(len(lonex))
    lat0=latex[argmin(lonex-CFlon[0])]
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([lat0,latex[ii]],[CFlon[0],lonex[ii]])[0][0]
        if lonex[ii]<CFlon[0]:
            distout[ii]=-1*distout[ii]
    return distout


dist14=getdist(lon14_ladcp,lat14_ladcp)

def getacrosstrack(u00,v00,theta00):
    uac00=u00*cos(theta00)+v00*sin(theta00)
    return uac00

v14_ladcp=getacrosstrack(ukn,vkn,theta)

dsortind14=argsort(dist14)
lon14_ladcp=lon14_ladcp[dsortind14]
lat14_ladcp=lat14_ladcp[dsortind14]
v14_ladcp=v14_ladcp[:,dsortind14]

dist14=sort(dist14)
ladcp14=xr.Dataset({'across track velocity [m/s]': (['pressure [db]','distance [km]'],  v14_ladcp)},
                coords={'distance [km]': dist14,
                        'pressure [db]': prsvec})

d1=0
d2=18
mid_dist_ladcp14=hstack((diff(ladcp14['distance [km]'][d1:d2])[0],(diff(ladcp14['distance [km]'][d1:d2][:-1])+diff(ladcp14['distance [km]'][d1:d2])[1:])/2,diff(ladcp14['distance [km]'][d1:d2])[-1]))
middistmat_ladcp14=(tile(mid_dist_ladcp14,[prslen-1,1]))
depthdiffmat_ladcp14=(tile(diff(prsvec),[len(ladcp14['distance [km]'][d1:d2]),1])).T


trans_ladcp14=(ladcp14['across track velocity [m/s]'][:-1,d1:d2]*depthdiffmat_ladcp14*middistmat_ladcp14/1e3).sum('pressure [db]').sum('distance [km]')
trans_ladcp14

## Haven't looked at these ones
# mat16=sort(glob.glob(datadir+'Shipboard/LADCP/2016*/*'))
# dat16=io.loadmat(mat16[0])

#This version from Penny's gridded product (after Shipboard_LADCP)
# mat16=  ...2014_2016_gridded_fromPenny
# dat16=io.loadmat(mat16)
#
# u16=dat16['corr_lad_u']
# v16=dat16['corr_lad_v']
# z16=dat16['lad_z']
# lon16=dat16['lon'][0]
# lat16=dat16['lat'][0]
#
# lonind16=(lon16<-30.5) & (lon16>-44)
# u16=u16[:,lonind16]
# v16=v16[:,lonind16]
# z16=z16[:,lonind16]
# lon16=lon16[lonind16]
# lat16=lat16[lonind16]


#################################################################################
# Do the same for geostrophic velocity from CTD
#################################################################################

[prs_ctd,dist,middist,lon,lat,sal,tmp,pden,geoshear,geovel]=pickle.load(open('../../pickles/Shipboard/CTD_1416geo.pickle','rb'))

# geostrophic transport calcs
#going to focus on '2014_2' and '2016_6' here

key14='2014_2'
key16='2016_6'

mid_dist_geo={}
middistmat_geo={}
minind_geo={}
for key in [key14,key16]:
    mid_dist_geo[key]=diff(dist[key])
    middistmat_geo[key]=tile(mid_dist_geo[key],[len(prs_ctd),1])
    minind_geo[key]=where(geovel[key][1,:]>0)[0][0]

depthdiff_geo=2

minind_geo

transgeo={}
for key in minind_geo:
    transgeo[key]=nansum(geovel[key][:,:minind_geo[key]]*middistmat_geo[key][:,:minind_geo[key]]*depthdiff_geo)/1e3

transgeo

moorind_geo={}
for key in transgeo:
    moorind_geo[key]=zeros(len(middist[key]))
    for dd in range(len(distvec)):
        moorind_geo[key][dd]=int(argmin(abs(middist[key]-distvec[dd])))
    moorind_geo[key]=[int(kk) for kk in moorind_geo[key]]

transgeo_moorloc={}
for key in transgeo:
    transgeo_moorloc[key]=nansum(geovel[key][:,moorind_geo[key][0]:moorind_geo[key][1]]*middistmat_geo[key][:,moorind_geo[key][0]:moorind_geo[key][1]]*depthdiff_geo)/1e3

transgeo_moorloc

#calculate transport from CF1 and CF2 positions, using nominal width of 12km each

geovel_sub={}
transgeo_sub={}
for key in transgeo:
    geovel_sub[key]=array([geovel[key][:,dd] for dd in moorind_geo[key][:2]]).T
    transgeo_sub[key]=nansum(geovel_sub[key]*12*depthdiff_geo)/1e3


#################################################################################
#### Plots comparing shipboard and mooring
#################################################################################

grdat['distance [km]']

# minimal framework for adcp contouring
def contladcp(axspec,occ,field,uniname,salinst=0):
    axspec.contourf(grdat['distance [km]'],grdat['pressure [db]'],grdat[field].sel(occupation=occ),univec[uniname][1],cmap=univec[uniname][2],extend='both');
    axspec.contour(grdat['distance [km]'],grdat['pressure [db]'],grdat[field].sel(occupation=occ),univec[uniname][3],colors='k');
    axspec.contour(grdat['distance [km]'],grdat['pressure [db]'],grdat['density'].sel(occupation=occ),27.8,colors='purple',linewidth=3)
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    if salinst==1:
        plotinstpos(axspec,'sal')

#same for mooring snap
def contmoor_pnt(axspec,datechoice,field,uniname,salinst=0):
    im=axspec.contourf(daily.distance,daily.depth,daily[field][:,:,datechoice].T,univec[uniname][1],cmap=univec[uniname][2],extend='both');
    axspec.contour(daily.distance,daily.depth,daily[field][:,:,datechoice].T,univec[uniname][3],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axspec.contour(daily.distance,daily.depth,daily['potential density'][:,:,datechoice].T,27.8,colors='purple',linewidth=3)
    if salinst==1:
        plotinstpos(axspec,'sal')
    return im


#mooring version for longer mean!

def contmoor_ave(axspec,date1,date2,field,uniname,salinst=0):
    im=axspec.contourf(daily.distance,daily.depth,daily[field].sel(date=slice(date1,date2)).mean(dim='date').T,univec[uniname][1],cmap=univec[uniname][2],extend='both');
    axspec.contour(daily.distance,daily.depth,daily[field].sel(date=slice(date1,date2)).mean(dim='date').T,univec[uniname][3],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axspec.contour(daily.distance,daily.depth,daily['potential density'].sel(date=slice(date1,date2)).mean(dim='date').T,27.8,colors='purple',linewidth=3)
    if salinst==1:
        plotinstpos(axspec,'sal')
    return im

def contctd(axspec,key,uniname,ctdfld,salinst=0):
    axspec.contourf(dist[key],prs_ctd,ctdfld[key],univec[uniname][1],cmap=univec[uniname][2],extend='both');
    axspec.contour(dist[key],prs_ctd,ctdfld[key],univec[uniname][3],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    if salinst==1:
        plotinstpos(axspec,'sal')

# #adapted for geostrophic vel
# def contgeo(axspec,key):
#     axspec.contourf(dist[key][:-1],prs_ctd,geovel[key],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r,extend='both');
#     axspec.contour(dist[key][:-1],prs_ctd,geovel[key],[-0.4,-0.2,0],colors='k');
#     axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
#     [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
#
# key14='2014_3'
# key16='2016_6'
#
# def plt_moorgeo():
#     fs=11
#     fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8,6))
#     contgeo(ax1,key14)
#     ax1.set_xlim([-15,110])
#     ax1.set_ylim([1800,0])
#     im=contmoor_ave(ax2,'2014-08-07','2014-09-07')
#     contgeo(ax3,key16)
#     im=contmoor_ave(ax4,'2016-06-28','2016-07-28')
#     ax1.set_title('Geostrophic velocity',fontsize=fs+2)
#     ax2.set_title('Monthly average from moorings',fontsize=fs+2)
#     ax1.set_ylabel('August, 2014 \n \n depth [m]',fontsize=fs)
#     ax3.set_ylabel('July, 2016 \n \n depth [m]',fontsize=fs)
#     ax3.set_xlabel('distance [km]',fontsize=fs)
#     ax4.set_xlabel('distance [km]',fontsize=fs)
#     fig.subplots_adjust(hspace=0.1,wspace=0.1)
#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
#     fig.colorbar(im, cax=cbar_ax,label='across-track velocity [m/s]')
#
#     savefig('../../figures/compare_shipmoor/Slopecurr_geomoorcomp'+savename+'.pdf')
#     savefig('../../figures/compare_shipmoor/Slopecurr_geomoorcomp'+savename+'.png')
#
# plt_moorgeo()

univec['sal'][1]=array([32,33,33.6 , 34.  , 34.4 , 34.8 , 34.9 , 34.92, 34.94, 34.96, 34.98,
       35. ,35.02,35.04,35.06 ])

univec['sal'][3]=array([33.6 , 34.  , 34.4 , 34.8 , 34.9 , 34.92, 34.94, 34.96, 34.98,
       35.,35.02,35.04,35.06  ])

def plt_moorladcp(field1,field2,instlab,uniname,ctdfld=sal,smallsec=0,salinst=0):
    fs=11
    fig, ((ax1, ax2,ax3), (ax4,ax5, ax6)) = plt.subplots(2, 3, sharex=True, sharey=True,figsize=(12,6))
    if smallsec==0:
        contladcp(ax1,'2014 (KN221)',field1,uniname,salinst)
    else:
        contctd(ax1,key14,uniname,ctdfld,salinst)
    ax1.set_xlim([-15,110])
    ax1.set_ylim([1800,0])
    im=contmoor_pnt(ax2,0,field2,uniname,salinst)
    contmoor_ave(ax3,'2014-08-07','2014-09-07',field2,uniname,salinst)
    if smallsec==0:
        contladcp(ax4,'2016',field1,uniname,salinst)
    else:
        contctd(ax4,key16,uniname,ctdfld,salinst)
    im=contmoor_pnt(ax5,-1,field2,uniname,salinst)
    contmoor_ave(ax6,'2016-06-28','2016-07-28',field2,uniname,salinst)
    ax1.set_title(instlab,fontsize=fs+2)
    ax2.set_title('Daily snap from moorings',fontsize=fs+2)
    ax3.set_title('Monthly mean from moorings',fontsize=fs+2)
    ax1.set_ylabel('August, 2014 \n \n depth [m]',fontsize=fs)
    ax4.set_ylabel('July, 2016 \n \n depth [m]',fontsize=fs)
    if 'velocity' in field1:
        ax1.text(-5,1500,str('{:5.2f}'.format(trans_ladcp.sel(occupation='2014 (KN221)').values))+' Sv',color='white',zorder=1e3,fontsize=fs+2)
        ax2.text(-5,1500,str('{:5.2f}'.format(egic['trans'][0].values))+' Sv',color='white',zorder=1e3,fontsize=fs+2)
        ax3.text(-5,1500,str('{:5.2f}'.format(mean(egic['trans'][:30]).values))+' Sv',color='white',zorder=1e3,fontsize=fs+2)
        ax4.text(-5,1500,str('{:5.2f}'.format(trans_ladcp.sel(occupation='2016').values))+' Sv',color='white',zorder=1e3,fontsize=fs+2)
        ax5.text(-5,1500,str('{:5.2f}'.format(egic['trans'][-1].values))+' Sv',color='white',zorder=1e3,fontsize=fs+2)
        ax6.text(-5,1500,str('{:5.2f}'.format(mean(egic['trans'][-30:]).values))+' Sv',color='white',zorder=1e3,fontsize=fs+2)
    ax4.set_xlabel('distance [km]',fontsize=fs)
    ax5.set_xlabel('distance [km]',fontsize=fs)
    ax6.set_xlabel('distance [km]',fontsize=fs)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=field1)

    if smallsec==0:
        savefig('../../figures/compare_shipmoor/Slopecurr_'+uniname+savename+'.pdf')
        savefig('../../figures/compare_shipmoor/Slopecurr_'+uniname+savename+'.png')
    else:
        savefig('../../figures/compare_shipmoor/Slopecurr_'+uniname+'_smallscaleCTD'+savename+'.pdf')
        savefig('../../figures/compare_shipmoor/Slopecurr_'+uniname+'_smallscaleCTD'+savename+'.png')


plt_moorladcp('across track velocity [m/s]','across track velocity', 'LADCP','uacross')

plt_moorladcp('salinity','salinity', 'CTD','sal',smallsec=1,salinst=1)


plt_moorladcp('temperature [$^\circ$ C]','temperature', 'CTD','tmp',smallsec=1,salinst=1)

plt_moorladcp('temperature [$^\circ$ C]','temperature', 'CTD','tmp',salinst=1)


plt_moorladcp('salinity','salinity', 'CTD','sal',salinst=1)

XXXXXXXXXXXXXXXXXXX

#################################################################################
# Compare vertically integrated velocity profiles as a function of distance from the coast
#################################################################################

## This was not done properly! If to redo, make sure averages and integrals are pressure-grid aware!!

# def compdistvel(afunc,tit,ylab):
#     figure()
#     intdpth=300
#     plot(adcp_dist,afunc(adcp_14['across track velocity'][adcp_14['depth'][:,0]<intdpth,:],axis=0),label='August 2014')
#
#     plot(adcp_dist[:-3],afunc(adcp_16['across track velocity'][adcp_16['depth'][:,0]<intdpth,:],axis=0),label='August 2016')
#     plot(0,0,'k',label='Vessel mounted ADCP')
#     plot(0,0,'k--',label='Moored array measurements')
#     plot(daily.distance,afunc(daily['across track velocity'][:,daily.depth<intdpth,0],axis=1),'--',color='blue')
#     # plotinstpos('v')
#     plot(daily.distance,afunc(daily['across track velocity'][:,daily.depth<intdpth,-1],axis=1),'--',color='orange')
#     xlim([-15,110])
#     axhline(0,color='k')
#     xlabel('distance (km)')
#     ylabel(ylab)
#     legend(loc=(1.05,0.2))
#     title(tit)
#
# help(plotinstpos)
# daily.salinity.min()
#
# compdistvel(nanmean,'Average currents (to 300m)','velocity [m/s]')
# savefig('../../figures/shipboard/ADCPmoorcomp_avevel.png',bbox_inches='tight')
#
# # compdistvel(nansum,'Integrated currents (to 300m)','integrated velocity [m/s]')
# # savefig('../../figures/shipboard/ADCPmoorcomp_intvel.png',bbox_inches='tight')


#################################################################################
### Extra figures
#################################################################################

#
# figure(figsize=(10,3))
# psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
# psf(cc['trans plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# savefig('../../figures/xport/cc_trans_fixedpos_salcomp.pdf',bbox_inches='tight')
#
# figure(figsize=(10,3))
# psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
# psf(cc['trans min vel plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# savefig('../../figures/xport/cc_trans_minvel_salcomp.pdf',bbox_inches='tight')
#
#
# figure(figsize=(10,9))
# subplot(311)
# psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
# psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
# psf(cc['trans sal def'],'orange',-1.5,0,'cc_trans',labit='sal def')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# gca().set_xticklabels('')
# ylabel('Transport [Sv]')
# subplot(312)
# psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
# psf(cc['trans min vel plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# gca().set_xticklabels('')
# subplot(313)
# psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
# psf(cc['trans plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# savefig('../../figures/xport/cc_trans_allcomp.pdf',bbox_inches='tight')
#
# figure(figsize=(10,3))
# psf(cc['trans min vel plus sal'],'r',-1.5,0,'cc_trans',labit='min vel plus sal')
# psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
# psf(cc['trans plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# savefig('../../figures/xport/cc_trans_mvelsalpcomp.pdf',bbox_inches='tight')
#
#
# figure(figsize=(10,3))
# psf(mindxr,'purple',0,40,'position of surface velocity minimum')
# ylabel('')
# title('position of surface velocity minimum between currents [km]')
# savefig('../../figures/xport/minpos_smooth.pdf',bbox_inches='tight')
#
#
# figure(figsize=(10,3))
# psf(maxdxr,'grey',-5,30,'position of coastal current core (max out to 30km)')
# ylabel('')
# title('position of coastal current core (max out to 30km) [km]')
# savefig('../../figures/xport/maxpos_smooth.pdf',bbox_inches='tight')
