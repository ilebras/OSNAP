## Save properties in density space as well as gridded density fields for use in other plotting scripts (so I don't have to re-run all this...)
from aux_funcs import *

#mooring resolution gridded CF data
dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1808lpfilt.pickle','rb'))

# dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))
#OOI profiling mooring
ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb'))
##OOI merged density profiles (sfc mooring and flanking mooring A)
oom=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_denmerged_xray.pickle','rb'))

# univec['PV']=['PV']

denvec=arange(27.5,27.85,0.005)
middenvec=denvec[:-1]+diff(denvec)/2
varvec=['sal','tmp','dpth','thick']#,'PV'
denmats={}
for var in varvec:
    denmats[var]=zeros((len(dat.distance),len(middenvec),len(dat.date)))
    for moornum in range(8):
        for dd in range(len(middenvec)):
            if var=='dpth':
                denmats[var][moornum,dd,:]=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).min(dim='depth');
                minnan=dat['potential density'][moornum,:,:].min(dim='depth')>denvec[dd] #
                denmats[var][moornum,dd,minnan]=0
            elif var=='thick':
                alldpths=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1])
                denmats[var][moornum,dd,:]=alldpths.max(dim='depth')-alldpths.min(dim='depth');
            else:
                denmats[var][moornum,dd,:]=dat[univec[var][0]][moornum,:,:].where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');

ooidenmat={}
for var in varvec:
    ooidenmat[var]=zeros((len(middenvec),len(ooi.date))) # note: OOI data is every 12 hours instead of 24
    for dd in range(len(middenvec)):
        if var=='dpth':
            ooidenmat[var][dd,:]=ooi['depth'].where(ooi['potential density']>denvec[dd]).min(dim='prs');
            minnan=ooi['potential density'].min(dim='prs')>denvec[dd]
            ooidenmat[var][dd,minnan]=0
        elif var=='thick':
            alldpths=ooi.depth.where(ooi['potential density']>denvec[dd]).where(ooi['potential density']<=denvec[dd+1])
            ooidenmat[var][dd,:]=alldpths.max(dim='prs')-alldpths.min(dim='prs');
        elif var=='PV':
            ooidenmat[var][dd,:]=ooi[univec[var][0]].where(ooi['pden_mid']>denvec[dd]).where(ooi['pden_mid']<=denvec[dd+1]).mean(dim='prs_mid');
        else:
            ooidenmat[var][dd,:]=ooi[univec[var][0]].where(ooi['potential density']>denvec[dd]).where(ooi['potential density']<=denvec[dd+1]).mean(dim='prs');

oomden={}
for var in ['thick','dpth','sal','ptmp']:
    oomden[var]=zeros((len(middenvec),len(oom.date)))
    for dd in range(len(middenvec)):
        if var=='thick':
            alldpths=oom.depth.where(oom['pden']>denvec[dd]).where(oom['pden']<=denvec[dd+1])
            oomden['thick'][dd,:]=alldpths.max(dim='prs')-alldpths.min(dim='prs');
        elif var=='dpth':
            oomden[var][dd,:]=oom['depth'].where(oom['pden']>denvec[dd]).min(dim='prs');
            minnan=oom['pden'].min(dim='prs')>denvec[dd]
            oomden[var][dd,minnan]=0
        else:
            oomden[var][dd,:]=oom[var].where(oom['pden']>denvec[dd]).where(oom['pden']<=denvec[dd+1]).mean(dim='prs');


oom_dendat=xr.Dataset({'thickness': (['den','date'],  oomden['thick']),'depth': (['den','date'],  oomden['dpth']),'sal': (['den','date'],  oomden['sal']),'tmp': (['den','date'],  oomden['ptmp'])},coords={'den': middenvec,'date': oom.date.values,'distance': oom_dist})


def makedendat():
    dendat=xr.Dataset({'depth': (['distance','den_bnd','date'],  denmats['dpth']),
                       'thickness': (['distance','den','date'],  denmats['thick']),
                       # 'PV': (['distance','den','date'],  denmats['PV']),
                       'sal': (['distance','den','date'],  denmats['sal']),
                       'tmp': (['distance','den','date'],  denmats['tmp'])},
                       coords={'den': middenvec,
                               'den_bnd': denvec[:-1],
                               'distance': dat.distance.values,
                               'date': dat.date.values})
# 'PV': (['distance','den','date'],  denmats['PV']),

    return dendat

dendat=makedendat()


def make_ooi_dendat():
    dendat=xr.Dataset({'depth': (['den_bnd','date'],  ooidenmat['dpth']),
                       'thickness': (['den','date'],  ooidenmat['thick']),
                       # 'PV': (['den','date'],  ooidenmat['PV']),
                       'sal': (['den','date'],  ooidenmat['sal']),
                       'tmp': (['den','date'],  ooidenmat['tmp'])},
                       coords={'den': middenvec,
                               'den_bnd': denvec[:-1],
                               'date': ooi.date.values,
                               'distance': ooi_dist})
# 'PV': (['distance','den','date'],  denmats['PV']),

    return dendat

ooi_dendat=make_ooi_dendat()

pickle.dump([dendat,ooi_dendat,oom_dendat],open(datadir+'OSNAP2016recovery/pickles/convection_dengrid/dendats_cf-ooi-oom.pickle','wb'),protocol=2)


##############################################################################################################################################################################################################################################################
################################## Combine densities into one xarray for manip   #################################################
###############################################################################################################################
###############################################################################################################################


ooi_den_daily=ooi['potential density'].resample(date='1D').mean(dim='date')
oom_den_daily=oom['pden']

dmin=5
dmax=1975
grid_den={}
grid_dpth=range(dmin,dmax,25)
grid_den=zeros((7,len(grid_dpth),len(ooi_den_daily.date)))
grid_date=ooi_den_daily.date.values
for ii,dd in enumerate(grid_date):
    ## cf+m1
    jj=0
    for mm in range(2,8):
        if mm!=6:
            mden=dat['potential density'][mm,:,:].sel(date=dd)
            if sum(~isnan(mden))>1:
                nini=~isnan(mden).values
                nanind=argwhere(nini)[0][0]
                fm=interp1d(dat['depth'][nini],mden[nini],bounds_error=False)
                grid_den[jj,:,ii]=fm(grid_dpth)
            else:
                grid_den[jj,:,ii]=NaN*ones(len(grid_dpth))
            jj+=1


    ## oom
    oomden=oom_den_daily.sel(date=dd)
    if sum(~isnan(oomden))>1:
        nini=~isnan(oomden).values
        foom=interp1d(oom.depth[nini],oomden[nini],bounds_error=False)
        grid_den[-2,:,ii]=foom(grid_dpth)
    else:
        grid_den[-2,:,ii]=NaN*ones(len(grid_dpth))

    ## ooi
    ooiden=ooi_den_daily.sel(date=dd)
    if sum(~isnan(ooiden))>1:
        nini=~isnan(ooiden).values
        nanind=argwhere(nini)[0][0]
        fooi=interp1d(ooi.depth[nini],ooiden[nini],bounds_error=False)
        grid_den[-1,:,ii]=fooi(grid_dpth)
    else:
        grid_den[-1,:,ii]=NaN*ones(len(grid_dpth))

datelen=len(grid_date)
## now go through every depth level and interpolate horizontally
for jj,zz in enumerate(grid_dpth):
    for mm in range(7):
        ##ooi
        nani=~isnan(grid_den[mm,jj,:])
        if sum(nani)!=0:
            f_time=interp1d(arange(datelen)[nani],grid_den[mm,jj,nani],bounds_error=False)
            grid_den[mm,jj,:]=f_time(arange(datelen))


#make it into an xarray!
grden=xr.Dataset({'pden': (['distance','depth','date'],grid_den),'moor': (['distance'],['cf3','cf4','cf5','cf6','m1','ooi_merged','ooi_prf'])}, coords={'depth': grid_dpth,'date': ooi_den_daily.date.values,'distance':hstack(([distvec[2:6],distvec[-1],oom_dist,ooi_dist]))})

pickle.dump(grden,open(datadir+'OSNAP2016recovery/pickles/convection_dengrid/grden_cf3-ooi.pickle','wb'),protocol=2)



# Leave the geostrophic tangent here since all the necessary pieces are above and its a dead end for now so not worth organizing properly.
# Basically, the changing geostrophic velocities are not a big part of the story. Layer thickness changes are much more important.
# Also had a hard time getting a good reference, but its kind of a moot point.
# Plots are in MixedLayer/geotrans
# ###############################################################################################################################
# ###############################################################################################################################
# ############################################## Geostrophic velocity calc   #################################################
# ###############################################################################################################################
# ###############################################################################################################################
#
# fcor=gsw.f(60)
# diffden=squeeze(diff(grid_den[:-1,:,:],axis=0))
# geoshear=[diffden[ii,::-1,:]*9.8/fcor/1028/diff(grden.distance[:-1])[ii]/1e3 for ii in range(5)]
#
# geovel=nancumsum(geoshear,axis=1)*diff(grid_dpth[::-1])[0]
# # geovel[isnan(geoshear)]=NaN
# grid_dpth
#
# geodist=grden.distance[:-2]+diff(grden.distance[:-1])/2
#
# gvel=xr.Dataset({'u': (['distance','depth','date'],geovel),'ones': (['depth','date'],ones((len(grid_dpth),len(ooi_den_daily.date)))),'moor': (['distance'],['cf3-4','cf4-5','cf5-6','cf6-m1','m1-ooi_merged'])},coords={'depth': grid_dpth[::-1],'date': ooi_den_daily.date.values,'distance':geodist})
#
# ooi_date1=str(ooi_den_daily.date[0].values)[:10]
# ooi_date2=str(ooi_den_daily.date[-1].values)[:10]
#
#
# mvel=dat['across track velocity'].sel(date=slice(ooi_date1,ooi_date2))
# mveli=mvel[2:-1,:,:]+diff(mvel[2:,:,:],axis=0)/2
#
# cf6m1_vel=(mvel.mean(dim='date')[5,:]+mvel.mean(dim='date')[7,:])/2
#
# def comp_geoveli_m1():
#     plot(geovel_m1,gvel.depth)
#
#
# geovel={}
# for ii,mm in enumerate([4,5,7]):
#     geovel[mm]=gvel.u[ii+1,:,:]+(gvel.u[ii+2,:,:]-gvel.u[ii+1,:,:])*(mvel.distance[mm]-gvel.distance[ii+1])/(gvel.distance[ii+2]-gvel.distance[ii+1])
#
#
# def comp_geoveli_cf6():
#     plot(mean(geovel[5],axis=1),gvel.depth,label='geo vel at cf6')
#     plot(mean(gvel.u[3,:,:],axis=1),gvel.depth, label='geo vel CF6-M1')
#     plot(mean(gvel.u[2,:,:],axis=1),gvel.depth,label='geo vel CF5-CF6')
#     plot(mean(mvel[5,:,:]+0.2,axis=1),dat.depth,label='direct vel at CF6 + 0.2')
#     legend()
#     ylim(2000,0)
#
# comp_geoveli_cf6()
#
# #
# # def comp_geoveli_m1():
# #     plot(mean(geovel_m1,axis=1),gvel.depth,label='geo vel at M1')
# #     plot(mean(gvel.u[3,:,:],axis=1),gvel.depth, label='geo vel CF6-M1')
# #     plot(mean(gvel.u[4,:,:],axis=1),gvel.depth,label='geo vel M1-OOI')
# #     plot(mean(mvel[-1,:,:]+0.1,axis=1),dat.depth,label='direct vel at M1 + 0.1')
# #     legend()
# #     ylim(2000,0)
# #
# # comp_geoveli_m1()
#
#
#
# def comp_mveli_geovel():
#     for ii in range(4):
#         figure()
#         if ii==3:
#             plot(mveli[ii,:,:].mean(dim='date'),dat.depth,label='CF6-CF7')
#             plot(mveli[ii+1,:,:].mean(dim='date'),dat.depth,label='CF7-M1')
#             plot(mveli[ii:ii+2,:,:].mean(dim='distance').mean(dim='date'),dat.depth,label='CF6-CF7-M1')
#             plot(cf6m1_vel,dat.depth,label='CF6-M1 mean')
#             plot(nanmean(geovel[ii,:,:],axis=1)-0.175,grid_dpth[::-1],label='geovel between CF6-M1 minus 0.175')
#             gca().invert_yaxis()
#             legend()
#             figure()
#             plot(nanmean(geovel[ii,:,:],axis=1)-0.175,grid_dpth[::-1],label='geovel between CF6-M1 minus 0.175')
#             plot(mvel.mean(dim='date')[5,:],dat.depth,label='CF6 mean vel')
#             plot(mvel.mean(dim='date')[7,:],dat.depth,label='M1 mean vel')
#
#         else:
#             plot(mveli[ii,:,:].mean(dim='date')+nanmin(abs(mveli[ii,:,:].mean(dim='date'))),dat.depth,label='direct vel between moorings')
#             plot(nanmean(geovel[ii,:,:],axis=1),grid_dpth[::-1],label='geovel')
#             title(str(nanmin(abs(mveli[ii,:,:].mean(dim='date')))))
#         gca().invert_yaxis()
#         xlim(-0.4,0)
#         xlabel('speed [m/s]')
#         ylabel('depth [m]')
#         legend()
#
#
# comp_mveli_geovel()
#
# def comp_cont():
#     figure()
#     contourf(dat.distance,dat.depth,dat['across track velocity'].mean(dim='date').T,31,cmap=cm.YlGnBu_r,vmin=-0.5,vmax=0)
#     ylim(1500,0)
#     colorbar()
#     xlim(25,150)
#     figure()
#     contourf(geodist,grid_dpth[::-1],mean(geovel,axis=2).T-0.175,31,cmap=cm.YlGnBu_r,vmin=-0.5,vmax=0)
#     ylim(1500,0)
#     colorbar()
#     xlim(25,150)
#
#
# comp_cont()
#
# barovec=[-0.1,-0.15,-0.2,-0.175,-0.05]
#
# def comp_Atdepth():
#     title('~surface')
#     plot(gvel.distance,mean(gvel.u[:,-1,:],axis=1)+barovec,'o-',linewidth=3,label='geostrophic estimate');
#     plot(dat.distance[2:],dat['across track velocity'][2:,3,:].mean(dim='date'),'o-',linewidth=3,label='mooring vel');
#     legend()
#     savefig(figdir+'MixedLayer/geotrans/Velcomp_timemean_surface.png',bbox_inches='tight')
#     figure()
#     title('~500m')
#     plot(gvel.distance,mean(gvel.u[:,-21,:],axis=1)+barovec,'o-',linewidth=3,label='geostrophic estimate');
#     plot(dat.distance[2:],dat['across track velocity'][2:,254,:].mean(dim='date'),'o-',linewidth=3,label='mooring vel');
#     axhline(0)
#     legend()
#     savefig(figdir+'MixedLayer/geotrans/Velcomp_timemean_500m.png',bbox_inches='tight')
#     figure()
#     title('~1000m')
#     plot(gvel.distance,mean(gvel.u[:,38,:],axis=1)+barovec,'o-',linewidth=3,label='geostrophic estimate');
#     plot(dat.distance[2:],dat['across track velocity'][2:,508,:].mean(dim='date'),'o-',linewidth=3,label='mooring vel');
#     axhline(0)
#     legend()
#     savefig(figdir+'MixedLayer/geotrans/Velcomp_timemean_1000m.png',bbox_inches='tight')
#     figure()
#     title('~1500m')
#     plot(gvel.distance,mean(gvel.u[:,18,:],axis=1)+barovec,'o-',linewidth=3,label='geostrophic estimate');
#     plot(dat.distance[2:],dat['across track velocity'][2:,761,:].mean(dim='date'),'o-',linewidth=3,label='mooring vel');
#     axhline(0)
#     savefig(figdir+'MixedLayer/geotrans/Velcomp_timemean_1500m.png',bbox_inches='tight')
#
# comp_Atdepth()
#
# ###############################################################################
# #####First, calculate the geostrophic transport of u and d IIW between CF5 and M1
# ###############################################################################
# btwndpths={}
# btwndpths['cf5-6']={}
# btwndpths['cf6-m1']={}
# btwndpths['m1-ooi']={}
# for ii,dd in enumerate([p1,p2,p3]):
#     btwndpths['cf5-6'][ii]=(dendat.depth[4,dd,:]+dendat.depth[5,dd,:])/2
#     btwndpths['cf6-m1'][ii]=(dendat.depth[5,dd,:]+dendat.depth[7,dd,:])/2
#     btwndpths['m1-ooi'][ii]=(dendat.depth[7,dd,:]+oom_dendat.depth[dd,:])/2
#
# #need to make corrections to deal with outcrops
# #step 1: if outcrops at offshore point, make the between depth 0. I will then subtract
# for ii,dd in enumerate([p1,p2,p3]):
#     btwndpths['cf5-6'][ii][dendat.depth[5,dd,:]==0]=0
#     btwndpths['cf6-m1'][ii][dendat.depth[7,dd,:]==0]=0
#     btwndpths['m1-ooi'][ii][oom_dendat.depth[dd,:]==0]=0
#
# btwndpths['cf5-6'][0].plot()
# btwndpths['cf5-6'][0].plot()
# plot(xout['cf5-6']['upper'],'o')
#
#
# surfden={}
# for var in ['cf5','cf6','m1','ooi_merged']:
#     surfden[var]=grden.where(grden.moor==var).sel(depth=30).dropna(dim='distance').pden.values.flatten()
#
#
# ## get horizontal position of outcrop between moorings for the outcrop correction transport
#
# xout={}
# xout['cf5-6']={}
# xout['cf5-6']['upper']=(distvec[5]-distvec[4])*(d1-surfden['cf5'])/(surfden['cf6']-surfden['cf5'])
# xout['cf5-6']['upper'][surfden['cf6']<d1]=0
# # plot(xout['cf5-6']['upper'],'o')
# xout['cf5-6']['deep']=0*surfden['cf5']
# xout['cf6-m1']={}
# xout['cf6-m1']['upper']=(distvec[7]-distvec[5])*(d1-surfden['cf6'])/(surfden['m1']-surfden['cf6'])
# xout['cf6-m1']['upper'][surfden['m1']<d1]=0
# xout['cf6-m1']['upper'][surfden['cf6']>d1]=0
# # plot(xout['cf6-m1']['upper'],'o')
# xout['cf6-m1']['deep']=(distvec[7]-distvec[5])*(d2-surfden['cf6'])/(surfden['m1']-surfden['cf6'])
# xout['cf6-m1']['deep'][surfden['m1']<d2]=0
# # plot(xout['cf6-m1']['deep'],'o')
# xout['m1-ooi']={}
# xout['m1-ooi']['upper']=(oom_dist-distvec[7])*(d1-surfden['m1'])/(surfden['ooi_merged']-surfden['m1'])
# xout['m1-ooi']['upper'][surfden['ooi_merged']<d1]=0
# xout['m1-ooi']['upper'][surfden['m1']>d1]=0
# # plot(xout['m1-ooi']['upper'],'o')
# xout['m1-ooi']['deep']=(oom_dist-distvec[7])*(d2-surfden['m1'])/(surfden['ooi_merged']-surfden['m1'])
# xout['m1-ooi']['deep'][surfden['ooi_merged']<d2]=0
# xout['m1-ooi']['deep'][surfden['m1']>d2]=0
# # plot(xout['m1-ooi']['deep'],'o')
#
#
# geodistvec=hstack((distvec,oom_dist))
#
# gdd=abs(hstack((diff(gvel.depth),25/2+5)))
#
#
# [uIIW,dIIW,IIW,egic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_xport/IIWtrans_direct.pickle','rb'))
#
# m_inner=[4,5,7]
# m_outer=[5,7,8]
#
# ulayer={}
# ulayer_baro={}
# outcropcorr={}
# outcropcorr_baro={}
# for mm,mchoo in enumerate(['cf5-6','cf6-m1','m1-ooi']):
#     bvel=barovec[mm+2]
#     for ii,lev in enumerate(['upper','deep']):
#         if lev=='upper':
#             p_up=p1
#             p_dn=p2
#         elif lev=='deep':
#             p_up=p2
#             p_dn=p3
#
#         ulayer[mchoo+'_'+lev]=(gvel.u[mm+2]).where(gvel.depth>btwndpths[mchoo][ii]).where(gvel.depth<btwndpths[mchoo][ii+1]).T
#         ulayer_baro[mchoo+'_'+lev]=bvel*(gvel.ones).where(gvel.depth>btwndpths[mchoo][ii]).where(gvel.depth<btwndpths[mchoo][ii+1]).T
#
#         outcrop_inner=((gvel.u[mm+2]).where(gvel.depth<dendat.depth[m_inner[mm],p_up,:]).T*gdd).sum(dim='depth')*(xout[mchoo][lev])/2e3
#         outcrop_inner_baro=(bvel*(gvel.ones).where(gvel.depth<dendat.depth[m_inner[mm],p_up,:]).T*gdd).sum(dim='depth')*(xout[mchoo][lev])/2e3
#         if lev=='upper':
#             dout=(geodistvec[m_outer[mm]]-geodistvec[m_inner[mm]])-xout[mchoo]['deep']
#             dout[xout[mchoo]['deep']==0]=0
#             if 'ooi' in mchoo:
#                 outcrop_outer=((gvel.u[mm+2]).where(gvel.depth<oom_dendat.depth[p_dn,:]).T*gdd).sum(dim='depth')*dout/1e3
#                 outcrop_outer_baro=(bvel*(gvel.ones).where(gvel.depth<oom_dendat.depth[p_dn,:]).T*gdd).sum(dim='depth')*dout/1e3
#             else:
#                 outcrop_outer=((gvel.u[mm+2]).where(gvel.depth<dendat.depth[m_outer[mm],p_dn,:]).T*gdd).sum(dim='depth')*dout/1e3
#                 outcrop_outer_baro=(bvel*(gvel.ones).where(gvel.depth<dendat.depth[m_outer[mm],p_dn,:]).T*gdd).sum(dim='depth')*dout/1e3
#             outcropcorr[mchoo+'_upper']=outcrop_inner+outcrop_outer
#             outcropcorr_baro[mchoo+'_upper']=outcrop_inner_baro+outcrop_outer_baro
#         else:
#             outcropcorr[mchoo+'_deep']=outcrop_inner
#             outcropcorr_baro[mchoo+'_deep']=outcrop_inner_baro
#
# xport_geo={}
# xport_baro={}
# for mm,mchoo in enumerate(['cf5-6','cf6-m1','m1-ooi']):
#     for ii,lev in enumerate(['upper','deep']):
#         # print(mchoo,geodistvec[m_outer[mm]],geodistvec[m_inner[mm]])
#         xport_geo[mchoo+'_'+lev]=-(ulayer[mchoo+'_'+lev]*gdd*(geodistvec[m_outer[mm]]-geodistvec[m_inner[mm]])).sum(dim='depth')/1e3+outcropcorr[mchoo+'_'+lev]
#         xport_baro[mchoo+'_'+lev]=-(ulayer_baro[mchoo+'_'+lev]*gdd*(geodistvec[m_outer[mm]]-geodistvec[m_inner[mm]])).sum(dim='depth')/1e3+outcropcorr_baro[mchoo+'_'+lev]
#
#
# uIIW['geotrans']=xport_geo['cf5-6_upper']+xport_geo['cf6-m1_upper']+xport_baro['cf5-6_upper']+xport_baro['cf6-m1_upper']
# dIIW['geotrans']=xport_geo['cf5-6_deep']+xport_geo['cf6-m1_deep']+xport_baro['cf5-6_deep']+xport_baro['cf6-m1_deep']
#
# uIIW['tottrans']=uIIW['geotrans']+xport_geo['m1-ooi_upper']+xport_baro['m1-ooi_upper']
# dIIW['tottrans']=dIIW['geotrans']+xport_geo['m1-ooi_deep']+xport_baro['m1-ooi_deep']
#
# def savefilt(field):
#     klist=list(field.keys())
#     for key in klist:
#         field[key+' filt']=sig.filtfilt(B,A,field[key])
#
#     return field
#
# egic=savefilt(egic)
# uIIW=savefilt(uIIW)
# dIIW=savefilt(dIIW)
# IIW=savefilt(IIW)
# xport_baro=savefilt(xport_baro)
# xport_geo=savefilt(xport_geo)
#
# def compgeo(var,ud):
#     figure(figsize=(12,3))
#     var['geotrans'].plot(alpha=0.5)
#     var['trans'].plot(alpha=0.5)
#     plot(gvel.date,var['geotrans filt'],color='C0',linewidth=3,label='geostrophic approx')
#     plot(gvel.date,outcropcorr['cf5-6_'+ud]+outcropcorr['cf6-m1_'+ud]+outcropcorr_baro['cf5-6_'+ud]+outcropcorr_baro['cf6-m1_'+ud],label='outcropcorr_tot')
#     # plot(gvel.date,xport_geo['cf5-6_'+ud]+xport_geo['cf6-m1_'+ud],label='geo component')
#     # plot(gvel.date,xport_baro['cf5-6_'+ud]+xport_baro['cf6-m1_'+ud],label='baro component')
#     # plot(gvel.date,xport_geo['cf5-6_'+ud]+xport_baro['cf5-6_'+ud],label='cf5-6 component')
#     # plot(gvel.date,xport_geo['cf6-m1_'+ud]+xport_baro['cf6-m1_'+ud],label='cf6-m1 component')
#     axhline(0)
#     # plot(gvel.date,outcropcorr_baro['cf6-m1_'+ud],label='outcropcorr_baro_cf6-m1')
#     plot(dat.date,var['trans filt'],color='C1',linewidth=3,label='direct calc')
#     legend(loc=(1.05,0))
#     # ylim(0,15)
#     figure(figsize=(12,3))
#     (var['geotrans']-var['trans']).plot(color='k')
#     axhline(0,color='grey')
#
#
# compgeo(uIIW,'upper')
#
# compgeo(dIIW,'deep')
#
#
# def tot_trans_decomp(var,ud):
#     figure(figsize=(12,3))
#     var['tottrans'].plot(alpha=0.5)
#     plot(gvel.date,xport_geo['cf5-6_'+ud]+xport_baro['cf5-6_'+ud],alpha=0.5,label='')
#     plot(gvel.date,xport_geo['cf6-m1_'+ud]+xport_baro['cf6-m1_'+ud],alpha=0.5,label='')
#     plot(gvel.date,xport_geo['m1-ooi_'+ud]+xport_baro['m1-ooi_'+ud],alpha=0.5,label='')
#     plot(gvel.date,var['tottrans filt'],color='C0',linewidth=3,label='total trans')
#     plot(gvel.date,xport_geo['cf5-6_'+ud+' filt']+xport_baro['cf5-6_'+ud+' filt'],linewidth=3,label='cf5-6 trans',color='C1')
#     plot(gvel.date,xport_geo['cf6-m1_'+ud+' filt']+xport_baro['cf6-m1_'+ud+' filt'],linewidth=3,label='cf6-m1 trans',color='C2')
#     plot(gvel.date,xport_geo['m1-ooi_'+ud+' filt']+xport_baro['m1-ooi_'+ud+' filt'],linewidth=3,label='m1-ooi trans',color='C3')
#     legend(loc=(1.05,0.2))
#     ylabel('[Sv]')
#     title(ud+' IIW transport: breakdown by mooring',fontsize=16)
#     axhline(0,color='grey')
#     savefig(figdir+'MixedLayer/geotrans/IIW_trans_bymoor_'+ud+'.png',bbox_inches='tight')
#     figure(figsize=(12,3))
#     var['tottrans'].plot(alpha=0.5)
#     plot(gvel.date,xport_geo['cf5-6_'+ud]+xport_geo['cf6-m1_'+ud]+xport_geo['m1-ooi_'+ud],alpha=0.5,label='')
#     plot(gvel.date,xport_baro['cf5-6_'+ud]+xport_baro['cf6-m1_'+ud]+xport_baro['m1-ooi_'+ud],alpha=0.5,label='')
#     plot(gvel.date,var['tottrans filt'],color='C0',linewidth=3,label='total trans')
#     plot(gvel.date,xport_geo['cf5-6_'+ud+' filt']+xport_geo['cf6-m1_'+ud+' filt']+xport_geo['m1-ooi_'+ud+' filt'],color='C1',linewidth=3,label='geostrophic, ref to 0 at bottom')
#     plot(gvel.date,xport_baro['cf5-6_'+ud+' filt']+xport_baro['cf6-m1_'+ud+' filt']+xport_baro['m1-ooi_'+ud+' filt'],color='C2',linewidth=3,label='barotropic correction')
#     legend(loc=(1.05,0.2))
#     ylabel('[Sv]')
#     title(ud+' IIW transport: breakdown into geostrophic components',fontsize=16)
#     axhline(0,color='grey')
#     savefig(figdir+'MixedLayer/geotrans/IIW_trans_bygeo_'+ud+'.png',bbox_inches='tight')
#     # plot(gvel.date,var['tottrans filt'],color='C0',linewidth=3,label='total trans')
#
# tot_trans_decomp(uIIW,'upper')
# tot_trans_decomp(dIIW,'deep')

#################################################################################
######## Some plotting beautification help for ref
####################################################################################
#
# def eachpanel(field,colo,axit,labit='',xr=daily,pnofilt=1,letlab='',ls='-'):
#     if pnofilt==1:
#         axit.plot(xr.date,field,alpha=0.5,color=colo,label='',linewidth=0.75)
#     axit.plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit,linestyle=ls)
#     axit.set_xlabel('')
#     axit.set_ylabel('')
#     axit.text(0.01, 0.85,letlab,transform = axit.transAxes,fontsize=15)
#
# years=matplotlib.dates.YearLocator()
# months=matplotlib.dates.MonthLocator()
# threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
# monthFMT=matplotlib.dates.DateFormatter('%B')
# yearFMT=matplotlib.dates.DateFormatter('\n %Y')
#
#
# def plottrans():
#     f,ax3=subplots(1,1,figsize=(8,2.5),sharex=True)
#     # eachpanel(egic['trans'],egicol,ax1,labit='Slope current transport')
#     # eachpanel(IIW['trans'],'k',ax1,labit='Total IIW transport')
#     eachpanel(uIIW['trans'],'k',ax3,labit='upper IIW transport')
#     eachpanel(dIIW['trans'],'C7',ax3,labit='deep IIW transport')
#     # ax1.legend(loc=2)
#     ax3.legend(loc=2)
#     ax3.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
#     ax3.xaxis.set_major_locator(years)
#     ax3.xaxis.set_minor_locator(threemonth)
#     ax3.xaxis.set_minor_formatter(monthFMT)
#     ax3.xaxis.set_major_formatter(yearFMT)
#     ax3.set_yticks(arange(0,15,4))
#     # ax1.set_ylabel('[Sv]')
#     ax3.set_ylabel('[Sv]')
#     # ax1.set_ylim(5,35)
#     savefig(figdir+'MixedLayer/paperfigs/IIW_trans_moorings.pdf',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/paperfigs/IIW_trans_moorings.png',bbox_inches='tight')
#
#
# plottrans()
#
# def comptrans_postcf5():
#     f,ax3=subplots(1,1,figsize=(8,2.5),sharex=True)
#     eachpanel(uIIW['trans'],'k',ax3,labit='all')
#     eachpanel(uIIW['trans cf5+'],'r',ax3,labit='from cf5')
#     title('upper IIW')
#     legend()
#     f,ax3=subplots(1,1,figsize=(8,2.5),sharex=True)
#     eachpanel(dIIW['trans'],'k',ax3,labit='all')
#     eachpanel(dIIW['trans cf5+'],'r',ax3,labit='from cf5')
#     legend()
#     title('deep IIW')
#
#
# comptrans_postcf5()

# would have to load IIW trans and OSNAP_MOC_context to plot this:
#
# figure(figsize=(12,3))
# (uIIW['tottrans']-uIIW['geotrans']).plot(alpha=0.5)
# plot(gvel.date,uIIW['tottrans filt']-uIIW['geotrans filt'],color='C0',linewidth=3,label='uIIW')
# (dIIW['tottrans']-dIIW['geotrans']).plot(alpha=0.5)
# plot(gvel.date,dIIW['tottrans filt']-dIIW['geotrans filt'],color='C1',linewidth=3,label='dIIW')
# legend(loc=(1.05,0.5))
# title('Between CF and gyre interior')
#
# figure(figsize=(12,3))
# uIIW['tottrans'].plot(alpha=0.5)
# plot(gvel.date,uIIW['tottrans filt'],color='C0',linewidth=3,label='uIIW (geostrophy)')
# dIIW['tottrans'].plot(alpha=0.5)
# plot(gvel.date,dIIW['tottrans filt'],color='C1',linewidth=3,label='dIIW (geostrophy)')
# uIIW['osnap fulltrans'].plot(marker='o',color='C0',label='OSNAP gridded')
# dIIW['osnap fulltrans'].plot(marker='o',color='C1',label='OSNAP gridded')
# legend(loc=(1.05,0.1))
# title('IIW transport up to OOI FLA mooring')
# ylabel('Transport [Sv]')
# xlabel('date')
# savefig(figdir+'MixedLayer/geotrans/IIW_trans_comp_toOOI.png',bbox_inches='tight')
