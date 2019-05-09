################################################################################
################################################################################
##### Compare sections and transports between shipboard and moored measurements
################################################################################
################################################################################

from aux_funcs import *

#################################################################################
# Load and calc trans from geostrophic velocity from CTD
#################################################################################
#
[prs_ctd,dist,middist,lon,lat,sal,tmp,pden,geoshear,geovel]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/Shipboard/CTD_1416geo.pickle','rb'))
#
# # geostrophic transport calcs
# #going to focus on '2014_2' and '2016_6' here
#
key14='2014_2'
key16='2016_6'
#
# mid_dist_geo={}
# middistmat_geo={}
# minind_geo={}
# for key in [key14,key16]:
#     mid_dist_geo[key]=diff(dist[key])
#     middistmat_geo[key]=tile(mid_dist_geo[key],[len(prs_ctd),1])
#     minind_geo[key]=where(geovel[key][1,:]>0)[0][0]
#
# depthdiff_geo=2
#
# minind_geo
#
# transgeo={}
# for key in minind_geo:
#     transgeo[key]=nansum(geovel[key][:,:minind_geo[key]]*middistmat_geo[key][:,:minind_geo[key]]*depthdiff_geo)/1e3
#
# transgeo
#
# moorind_geo={}
# for key in transgeo:
#     moorind_geo[key]=zeros(len(middist[key]))
#     for dd in range(len(distvec)):
#         moorind_geo[key][dd]=int(argmin(abs(middist[key]-distvec[dd])))
#     moorind_geo[key]=[int(kk) for kk in moorind_geo[key]]
#
# transgeo_moorloc={}
# for key in transgeo:
#     transgeo_moorloc[key]=nansum(geovel[key][:,moorind_geo[key][0]:moorind_geo[key][1]]*middistmat_geo[key][:,moorind_geo[key][0]:moorind_geo[key][1]]*depthdiff_geo)/1e3
#
# transgeo_moorloc
#
# #calculate transport from CF1 and CF2 positions, using nominal width of 12km each
#
# geovel_sub={}
# transgeo_sub={}
# for key in transgeo:
#     geovel_sub[key]=array([geovel[key][:,dd] for dd in moorind_geo[key][:2]]).T
#     transgeo_sub[key]=nansum(geovel_sub[key]*12*depthdiff_geo)/1e3
#
# transgeo_sub

################################################################################
##### Load ADCP and calc trans
################################################################################

adcp_dist=array(pd.DataFrame.from_csv(datadir+'Shipboard/adcp_distances.dat').index)



with open(datadir+'Shipboard/adcp_distances.dat', 'r') as f:
    reader = csv.reader(f)
    adcp_dist = list(reader)

adcp_dist=[float(dd[0]) for dd in adcp_dist]

adcp_dist_16=array(adcp_dist[:-3])

# adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')
adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/vmadcp/k2_2014_vmadcp_ilebras.mat')
plot(adcp_14['lon'],adcp_14['lat'],'o');
plot(adcp_14['lon'][0][-1],adcp_14['lat'][0][-1],'w*');


adcp_16=io.loadmat(datadir+'Shipboard/ar07_2016/ar07_2016_vm_adcp_ilebras.mat')


# plot(adcp_14['lon'],adcp_14['lat'],'ko');
# plot(adcp_16['lon'],adcp_16['lat'],'ro');
# plot(360+CFlon,CFlat,'yo');
# ylim([60,60.2])
# xlim([316.75,318])

def acrosstrackname(adic):
    adic['across track velocity']=-(adic['v_dt']*sin(theta)+adic['u_dt']*cos(theta))

    return adic

adcp_14=acrosstrackname(adcp_14)
adcp_16=acrosstrackname(adcp_16)

def getdist(lonex,latex):
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([latex[-1],latex[ii]],[lonex[-1],lonex[ii]])[0][0]
    return distout

adcp_dist_14=getdist(adcp_14['lon'][0],adcp_14['lat'][0])+min(dist[key14]) #to conform with CTD dist.

del adcp_dist

# adcp transport calcs
mid_dist_adcp14=hstack((diff(adcp_dist_14)[0]/2,(diff(adcp_dist_14)[:-1]+diff(adcp_dist_14)[1:])/2,diff(adcp_dist_14)[-1]/2))
mid_dist_adcp16=hstack((diff(adcp_dist_16)[0]/2,(diff(adcp_dist_16)[:-1]+diff(adcp_dist_16)[1:])/2,diff(adcp_dist_16)[-1]/2))

middistmat_a14=tile(mid_dist_adcp14,[len(adcp_14['depth'][:,0]),1])
shape(middistmat_a14)
middistmat_a16=tile(mid_dist_adcp16,[len(adcp_16['depth'][:,0]),1])
shape(middistmat_a16)
depthdiffmat_a14=tile(diff(hstack((0,adcp_14['depth'][:,0]))),[len(adcp_dist_14),1]).T
shape(depthdiffmat_a14)
depthdiffmat_a16=tile(diff(hstack((0,adcp_16['depth'][:,0]))),[len(adcp_dist_16),1]).T
shape(depthdiffmat_a16)

# minind14=where(adcp_14['across track velocity'][0,-20:]==max(adcp_14['across track velocity'][0,-20:]))[0][0]+len(adcp_14['across track velocity'][0,:])-20
# minind16=where(adcp_16['across track velocity'][0,-10:-7]==max(adcp_16['across track velocity'][0,-10:-7]))[0][0]+len(adcp_16['across track velocity'][0,:])-10
#
# trans14adcp_less=nansum(adcp_14['across track velocity'][:,minind14:-3]*depthdiffmat_a14[:,minind14:-3]*middistmat_a14[:,minind14:-3])/1e3
# trans14adcp=nansum(adcp_14['across track velocity'][:,minind14:]*depthdiffmat_a14[:,minind14:]*middistmat_a14[:,minind14:])/1e3
#
# trans16adcp=nansum(adcp_16['across track velocity'][:,minind16:]*depthdiffmat_a16[:,minind16:]*middistmat_a16[:,minind16:])/1e3
## transport from subsampling at mooring locations
# adcp_dist_sub=[adcp_dist[dd] for dd in moorind]
# mid_dist_adcp_sub=hstack((diff(adcp_dist_sub)[0]/2,(diff(adcp_dist_sub)[:-1]+diff(adcp_dist_sub)[1:])/2,diff(adcp_dist_sub)[-1]/2))
#
#
# # subsample adcp dailya at mooring locations
# adcpvel_14_sub=[adcp_14['across track velocity'][:,dd] for dd in moorind14[-2:]]
# adcpvel_16_sub=[adcp_16['across track velocity'][:,dd] for dd in moorind16[-2:]]
# #calculate transport from CF1 and CF2 positions, using nominal width of 12km each
# depthdiffmat_a14_sub=tile(diff(hstack((0,adcp_14['depth'][:,0]))),[2,1])
# depthdiffmat_a16_sub=tile(diff(hstack((0,adcp_16['depth'][:,0]))),[2,1])
#
# trans14adcp_sub=nansum(adcpvel_14_sub*depthdiffmat_a14_sub*12)/1e3
# trans16adcp_sub=nansum(adcpvel_16_sub*depthdiffmat_a16_sub*12)/1e3


# def getmoorind(adcp_dist_func):
#     moorind=zeros(len(distvec))
#     for dd in range(len(distvec)):
#         moorind[dd]=argmin(abs(adcp_dist_func-distvec[dd]))
#
#     moorind=[int(dd) for dd in moorind][::-1]
#     return moorind
# moorind14=getmoorind(adcp_dist_14)
# moorind16=getmoorind(adcp_dist_16)

min(adcp_dist_14)
min(adcp_dist_16)

adcpind_14=adcp_dist_14<15
adcp_dist_14[adcpind_14]

adcpind_16=adcp_dist_16<25
adcp_dist_16[adcpind_16]

## transport between mooring locations (including as much shelf as possible/fair comparison of range)
trans14adcp=nansum(adcp_14['across track velocity'][:,adcpind_14]*depthdiffmat_a14[:,adcpind_14]*middistmat_a14[:,adcpind_14])/1e3
trans16adcp=nansum(adcp_16['across track velocity'][:,adcpind_16]*depthdiffmat_a16[:,adcpind_16]*middistmat_a16[:,adcpind_16])/1e3

## read out the values
trans14adcp

trans16adcp

def coastal_cont(field,distfield,depthfield):
    contourf(distfield,depthfield,field,vmin=-0.8,vmax=0,cmap=cm.Blues_r)
    ylim(180,0)
    colorbar()

coastal_cont(adcp_14['across track velocity'][:,adcpind_14],adcp_dist_14[adcpind_14],adcp_14['depth'][:,-1])

coastal_cont(adcp_16['across track velocity'][:,adcpind_16],adcp_dist_16[adcpind_16],adcp_16['depth'][:,-1])



#################################################################################
# Combine shipboard ADCP and CTD to get freshwater transport
#################################################################################

## Note, for 2014 using offset salinity/velprs because it goes farther onto shelf
figure()
plot(adcp_14['lon']-360,adcp_14['lat'],'bo');
plot(adcp_14['lon'][0]-360,adcp_14['lat'][0],'bo',label='adcp');
plot(lon[key14],lat[key14],'ro',label='ctd');
plot(CFlon,CFlat,'g*',label='moorings');
legend()
## For 2016, going to repeat shelfmost measurement of vel to get all that freshwater
figure()
title('2016')
plot(lon[key16],lat[key16],'ro',label='ctd');
plot(adcp_16['lon'][0]-360,adcp_16['lat'][0],'bo',label='adcp');
plot(adcp_16['lon']-360,adcp_16['lat'],'bo');
plot(CFlon,CFlat,'g*',label='moorings');
legend()

adcp_dist_14
dist[key14]
adcp_dist_16
dist[key16]

#interpolate adcp onto salinity distances and depths

def interptosal(velfield,veldistvec,velprs,salfield,newdist):
    [prslen,distlen]=shape(velfield)
    # newdist=saldist#[saldist>=min(veldistvec)]
    newprs=prs_ctd[prs_ctd<=max(velprs)]
    salvel1=zeros((prslen,len(newdist[newdist<max(veldistvec)])))
    salvel2=zeros((len(newprs),len(newdist[newdist<max(veldistvec)])))
    #interpolate to distances
    for pp in range(prslen):
        if min(newdist)<min(veldistvec):
            distfunc=interpolate.interp1d(hstack((min(newdist),veldistvec)),hstack((0,velfield[pp,:])))
        else:
            distfunc=interpolate.interp1d(veldistvec,velfield[pp,:])
        salvel1[pp,:]=distfunc(newdist[newdist<max(veldistvec)])
    #interpolate to depth
    for dd in range(len(newdist[newdist<max(veldistvec)])):
        prsfunc=interpolate.interp1d(hstack((0,velprs[~isnan(salvel1[:,dd])],180)),hstack((salvel1[0,dd],salvel1[~isnan(salvel1[:,dd]),dd],salvel1[~isnan(salvel1[:,dd]),dd][-1])))
        salvel2[:,dd]=prsfunc(newprs)

    return salvel2,newprs,newdist[newdist<max(veldistvec)]


salvel={}
salprs={}
salvel[key14],salprs[key14],dist[key14]=interptosal(adcp_14['across track velocity'][adcp_14['depth'][:,0]<180,:],adcp_dist_14,adcp_14['depth'][:,0][adcp_14['depth'][:,0]<180],sal[key14],dist[key14])
salvel[key16],salprs[key16],dist[key16]=interptosal(adcp_16['across track velocity'][adcp_16['depth'][:,0]<180,:],adcp_dist_16,adcp_16['depth'][:,0][adcp_16['depth'][:,0]<180],sal[key16],dist[key16])


figure()
plot(adcp_dist_14,adcp_14['across track velocity'][0,:].T,label='adcp grid',color='b')
xlim([-20,50])
plot(dist[key14][:17],salvel[key14][0,:17].T,label='ctd grid',color='r')
legend()

figure()
plot(adcp_dist_16,adcp_16['across track velocity'][0,:].T,label='adcp grid',color='b')
xlim([-20,50])
plot(dist[key16][:17],salvel[key16][0,:17].T,label='ctd grid',color='r')
legend()



ctdind={}
ctdind[key14]=(dist[key14]<15)&(dist[key14]>-10)
ctdind[key16]=(dist[key16]<15)&(dist[key16]>-10)

for kk in ctdind:
    print(kk)
    print(dist[kk][ctdind[kk]])
    sal[kk]=sal[kk][:,:len(dist[kk])]

sref=34.9

# freshwater transport calcs
middist_fresh={}
middistmat_fresh={}
depthdiffmat_fresh={}
shipfresh={}
transtest={}
# minind_fresh={}
# moordistind={}
# shipfresh_moorloc={}
for kk in [key14,key16]:
    middist_fresh[kk]=hstack((diff(dist[kk])[0]/2,(diff(dist[kk])[:-1]+diff(dist[kk])[1:])/2,diff(dist[kk])[-1]/2))
    middistmat_fresh[kk]=tile(middist_fresh[kk],[len(salprs[kk]),1])
    depthdiffmat_fresh[kk]=tile(diff(hstack((0,salprs[kk]))),[len(dist[kk]),1]).T
    # minind_fresh[kk]=where(salvel[kk][0,2:17]==max(salvel[kk][0,2:17]))[0][0]+3
    transtest[kk]=nansum(salvel[kk][:,ctdind[kk]]*depthdiffmat_fresh[kk][:,ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]])/1e3
    shipfresh[kk]=nansum(salvel[kk][:,ctdind[kk]]*(sal[kk][:len(salprs[kk]),ctdind[kk]]-sref)/sref*depthdiffmat_fresh[kk][:,ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]])
    # moordistind[kk]=where((dist[kk]>=0)&(dist[kk]<=distvec[1]))[0]
    # shipfresh_moorloc[kk]=nansum(salvel[kk][:,moordistind[kk]]*(sal[kk][:len(salprs[kk]),moordistind[kk]]-sref)/sref*depthdiffmat_fresh[kk][:,moordistind[kk]]*middistmat_fresh[kk][:,moordistind[kk]])


shipfresh

transtest

versname='1810JHcal'

daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_'+versname+'.pickle','rb'))
daily['across track velocity']=-1*daily['across track velocity']
# versname='1808_CF3'
[cc,egic,eg,ic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/transdic_'+versname+'.pickle','rb'))


cc['trans'].mean()
cc['freshb'].mean()


cc['trans'][0].values
mean(cc['trans'][:30]).values
min(cc['trans'][:30]).values
max(cc['trans'][:30]).values

cc['trans'][-1].values
mean(cc['trans'][-30:]).values
min(cc['trans'][-30:]).values
max(cc['trans'][-30:]).values

cc['freshb'][0].values
mean(cc['freshb'][:30]).values
min(cc['freshb'][:30]).values
max(cc['freshb'][:30]).values

cc['freshb'][-1].values
mean(cc['freshb'][-30:]).values
min(cc['freshb'][-30:]).values
max(cc['freshb'][-30:]).values

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth),len(daily.date),1])),(2,0,1))

sep=9
def salvelcomp_paper():
        figure(figsize=(6.5,4))
        ls=3
        dind=86
        ccmin1=cc['trans']==min(cc['trans'][:30])
        ccmax1=cc['trans']==max(cc['trans'][:30])
        ccmin2=cc['trans']==min(cc['trans'][-30:])
        ccmax2=cc['trans']==max(cc['trans'][-30:])
        subplot(121)
        kk=key14
        salmn={}
        salmn['ship14']=sum(sal[kk][:len(salprs[kk]),ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]],axis=1)/sum(middistmat_fresh[kk][:,ctdind[kk]],axis=1)
        plot(salmn['ship14'],salprs[kk],color='orange',label='Shipboard, Aug 2014',linewidth=ls)
        kk=key16
        salmn['ship16']=sum(sal[kk][:len(salprs[kk]),ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]],axis=1)/sum(middistmat_fresh[kk][:,ctdind[kk]],axis=1)
        plot(salmn['ship16'],salprs[kk],color='red',linewidth=ls,label='Shipboard, July 2016')
        salmn['moor14']=(daily.salinity[:sep,:dind,:30].mean(dim='date')*middistmat_plus[:sep,:dind,0]).sum(dim='distance')/sum(middistmat_plus[:sep,:dind,0],axis=0).T
        plot(salmn['moor14'],daily.depth[:dind],ccol,label='Mooring, Aug 2014',linewidth=ls)
        salmn['moor16']=(daily.salinity[:sep,:dind,-30:].mean(dim='date')*middistmat_plus[:sep,:dind,0]).sum(dim='distance')/sum(middistmat_plus[:sep,:dind,0],axis=0).T
        plot(salmn['moor16'],daily.depth[:dind],'darkblue',label='Mooring, July 2014',linewidth=ls)
        ylabel('depth [m]')
        ylim([180,0])
        xlabel('mean salinity')
        legend()
        yticks([0,50,100,150])
        subplot(122)
        velmn={}
        velmn['moor14']=(daily['across track velocity'][:sep,:dind,:30].mean(dim='date')*middistmat_plus[:sep,:dind,0]).sum(dim='distance')/sum(middistmat_plus[:sep,:dind,0],axis=0).T
        plot(velmn['moor14'],daily.depth[:dind],ccol,linewidth=ls)
        velmn['moor16']=(daily['across track velocity'][:sep,:dind,-30:].mean(dim='date')*middistmat_plus[:sep,:dind,0]).sum(dim='distance')/sum(middistmat_plus[:sep,:dind,0],axis=0).T
        plot(velmn['moor16'],daily.depth[:dind],'darkblue',linewidth=ls)
        kk=key14
        velmn['ship14']=sum(salvel[kk][:len(salprs[kk]),ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]],axis=1)/sum(middistmat_fresh[kk][:,ctdind[kk]],axis=1)
        plot(velmn['ship14'],salprs[kk],color='orange',linewidth=ls)
        kk=key16
        velmn['ship16']=sum(salvel[kk][:len(salprs[kk]),ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]],axis=1)/sum(middistmat_fresh[kk][:,ctdind[kk]],axis=1)
        plot(velmn['ship16'],salprs[kk],color='red',linewidth=ls)
        gca().set_yticklabels('')
        ylim(180,0)
        xlabel('mean velocity [m/s]')
        tight_layout()
        savefig('../figures/paperfigs/salvelcomp_coastal.pdf',bbox_inches='tight')

        return salmn,velmn

salmn,velmn=salvelcomp_paper()


nancumsum(velmn['moor14']*(salmn['moor14']-sref)/sref*2*25,axis=0)[-1]
nancumsum(velmn['ship14']*(salmn['ship14']-sref)/sref*2*25,axis=0)[-1]
nancumsum(velmn['moor14']*(salmn['ship14']-sref)/sref*2*25,axis=0)[-1]
nancumsum(velmn['ship14']*(salmn['moor14']-sref)/sref*2*25,axis=0)[-1]

nancumsum(velmn['moor16']*(salmn['moor16']-sref)/sref*2*25,axis=0)[-1]
nancumsum(velmn['ship16']*(salmn['ship16']-sref)/sref*2*25,axis=0)[-1]
nancumsum(velmn['moor16']*(salmn['ship16']-sref)/sref*2*25,axis=0)[-1]
nancumsum(velmn['ship16']*(salmn['moor16']-sref)/sref*2*25,axis=0)[-1]

def plotfreshcum():
    dind=88
    ls=2
    freshfirst=(daily['across track velocity'].where(daily.distance<=distvec[1])[:,:dind,:30]*(daily.salinity.where(daily.distance<=distvec[1])[:,:dind,:30]-sref)/sref).mean(dim='date')
    freshlast=(daily['across track velocity'].where(daily.distance<=distvec[1])[:,:dind,-30:]*(daily.salinity.where(daily.distance<=distvec[1])[:,:dind,-30:]-sref)/sref).mean(dim='date')

    freshfirst=(daily['across track velocity'].where(daily.distance<=distvec[1])[:,:dind,:30].mean(dim='date')*(daily.salinity.where(daily.distance<=distvec[1])[:,:dind,:30].mean(dim='date')-sref)/sref)
    freshlast=(daily['across track velocity'].where(daily.distance<=distvec[1])[:,:dind,-30:].mean(dim='date')*(daily.salinity.where(daily.distance<=distvec[1])[:,:dind,-30:].mean(dim='date')-sref)/sref)

    fprof1=(freshfirst*middistmat_plus[:,:dind,0]).sum(dim='distance')
    fprof2=(freshlast*middistmat_plus[:,:dind,0]).sum(dim='distance')
    plot((fprof1[::-1]*2).cumsum(dim='depth'),daily.depth[:dind][::-1],ccol,linewidth=ls)
    plot((fprof2[::-1]*2).cumsum(dim='depth'),daily.depth[:dind][::-1],'darkblue',linewidth=ls)
    print('2014 mooring: ',(fprof1[::-1]*2).cumsum(dim='depth')[-1].values)
    print('2016 mooring: ',(fprof2[::-1]*2).cumsum(dim='depth')[-1].values)
    kk=key14
    fprofship1=sum(salvel[kk][:len(salprs[kk]),ctdind[kk]]*(sal[kk][:len(salprs[kk]),ctdind[kk]]-sref)/sref*middistmat_fresh[kk][:,ctdind[kk]],axis=1)
    plot(nancumsum(fprofship1[::-1]*2),salprs[kk][::-1],color='orange',linewidth=ls)
    print('2014 shipboard: ',nancumsum(fprofship1[::-1]*2)[-1])
    kk=key16
    fprofship2=sum(salvel[kk][:len(salprs[kk]),ctdind[kk]]*(sal[kk][:len(salprs[kk]),ctdind[kk]]-sref)/sref*middistmat_fresh[kk][:,ctdind[kk]],axis=1)
    plot(nancumsum(fprofship2[::-1]*2),salprs[kk][::-1],color='red',linewidth=ls)
    print('2016 shipboard: ',nancumsum(fprofship2[::-1]*2)[-1])
    # kk=key16
    # plot(sum(salvel[kk][:len(salprs[kk]),ctdind[kk]]*middistmat_fresh[kk][:,ctdind[kk]],axis=1)/sum(middistmat_fresh[kk][:,ctdind[kk]],axis=1),salprs[kk],color='red',linewidth=ls)
    # gca().set_yticklabels('')
    ylim([180,0])
    xlabel('cumulative freshwater transport [mSv]')

plotfreshcum()

#################################################################################
#### Plots comparing shipboard and mooring
#################################################################################

# # minimal framework for adcp contouring
# def contadcp(axspec,distveca,adic):
#     axspec.contourf(distveca,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
#     axspec.contour(distveca,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2,0],colors='k');
#     axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
#     [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]


#for the version interpolated to ctd
def contadcp_ctd(axspec,kk):
    axspec.contourf(dist[kk],salprs[kk],salvel[kk],arange(0,0.805,0.1),cmap=cm.Blues);
    axspec.contour(dist[kk],salprs[kk],salvel[kk],arange(0,0.805,0.1),colors='k');
    # axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]

# #same for mooring snap
# def contmoor_pnt(axspec,datechoice):
#     im=axspec.contourf(daily.distance,daily.depth,daily['across track velocity'][:,:,datechoice].T,arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
#     axspec.contour(daily.distance,daily.depth,daily['across track velocity'][:,:,datechoice].T,[-0.4,-0.2,0],colors='k');
#     axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
#     [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
#
#     return im


#mooring version for longer mean!

def contmoor_ave(axspec,date1,date2):
    im=axspec.contourf(daily.distance,daily.depth,daily['across track velocity'].sel(date=slice(date1,date2)).mean(dim='date').T,arange(0,0.805,0.1),cmap=cm.Blues);
    axspec.contour(daily.distance,daily.depth,daily['across track velocity'].sel(date=slice(date1,date2)).mean(dim='date').T,arange(0,0.805,0.1),colors='k');
    # axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    return im

#adapted for geostrophic vel
# def contgeo(axspec,key):
#     axspec.contourf(dist[key][:-1],prs_ctd,geovel[key],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r,extend='both');
#     axspec.contour(dist[key][:-1],prs_ctd,geovel[key],[-0.4,-0.2,0],colors='k');
#     axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
#     [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]


def mkbox(axc,x1,x2,y1,y2):
    axc.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],color='purple',linewidth=3)

def plt_mooradcp():
    md2=13
    fs=11
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8,6))
    contadcp_ctd(ax1,key14)
    ax1.set_xlim([-15,40])
    ax1.set_ylim([200,-5])
    # im=contmoor_pnt(ax2,0)
    im=contmoor_ave(ax2,'2014-8-17','2014-9-17')
    contadcp_ctd(ax3,key16)
    # im=contmoor_pnt(ax4,-1)
    im=contmoor_ave(ax4,'2014-6-28','2016-7-28')
    ax1.set_title('Vessel mounted ADCP',fontsize=fs+2)
    ax2.set_title('Monthly mean from moorings',fontsize=fs+2)
    ax1.set_ylabel('August, 2014 \n \n depth [m]',fontsize=fs)
    ax3.set_ylabel('July, 2016 \n \n depth [m]',fontsize=fs)
    mkbox(ax1,dist[key14][ctdind[key14]][0],dist[key14][ctdind[key14]][-1],180,0)
    mkbox(ax2,-10,md2,180,0)
    mkbox(ax3,dist[key16][ctdind[key16]][0],dist[key16][ctdind[key16]][-1],180,0)
    mkbox(ax4,-10,md2,180,0)
    ax3.set_xlabel('distance [km]',fontsize=fs)
    ax4.set_xlabel('distance [km]',fontsize=fs)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax,label='along-stream velocity [m/s]')

    savefig('../figures/paperfigs/Coastalcurr_adcpmoorcomp.pdf')
    savefig('../figures/compare_shipmoor/Coastalcurr_adcpmoorcomp.png')

    savefig('home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/Coastalcurr_adcpmoorcomp.pdf',bbox_inches='tight')

plt_mooradcp()


XXXXXXXXXXX

# def salvelcomp_paper():
#             figure(figsize=(9,5))
#             subplot(121)
#             kk=key14
#             plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'-',color='grey',label='Deployment, 2014')
#             plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'--',color='grey',label='Recovery, 2016')
#             plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'r',label='Full Shipboard')
#             kk=key16
#             plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'r--')
#             kk=key14
#             plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk]]).mean(axis=1),salprs[kk],color='orange',label='Shipboard within mooring range')
#             kk=key16
#             plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk]]).mean(axis=1),salprs[kk],color='orange',linestyle='--',label='')
#             kk=key14
#             plot(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k',label='Mooring daily average')
#             kk=key16
#             plot(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k--',label='')
#             plot(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').mean(dim='distance').T,daily.depth,'b',label='Mooring monthly average')
#             plot(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').mean(dim='distance').T,daily.depth,'b--',label='')
#             ylabel('depth [m]')
#             ylim([190,0])
#             xlabel('mean coastal current salinity')
#             legend()
#             subplot(122)
#             for kk in [key14,key16]:
#                 if kk==key14:
#                     kline='-'
#                 elif kk==key16:
#                     kline='--'
#                 plot(mean(salvel[kk][:len(salprs[kk]),:minind_fresh[kk]],axis=1),salprs[kk],'r',linestyle=kline)
#                 plot(mean(salvel[kk][:len(salprs[kk]),moordistind[kk]],axis=1),salprs[kk],color='orange',linestyle=kline)
#                 plot(daily['across track velocity'].isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k',linestyle=kline,)
#                 if kk==key14:
#                     plot(daily['across track velocity'].where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline)
#                 if kk==key16:
#                     plot(daily['across track velocity'].where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline)
#                 gca().set_yticklabels('')
#                 ylim([190,0])
#                 tight_layout()
#                 xlabel('mean coastal current velocity [m/s]')
#             savefig('../figures/paperfigs/salvelcomp_coastal.pdf',bbox_inches='tight')
#
#
# salvelcomp_paper()
# def plotsalvelprofs():
#     nn=1
#     figure(figsize=(8,8))
#     for kk in [key14,key16]:
#             subplot(2,2,nn)
#             plot(squeeze(sal[kk][:len(salprs[kk]),0]),salprs[kk],color='orange',label='Shipboard',linewidth=2)#: '+str('{:5.2f}'.format(shipfresh[kk]))+' mSv')
#             plot(squeeze(sal[kk][:len(salprs[kk]),ctdind[kk]]),salprs[kk],color='orange',linewidth=2)
#             if kk==key14:
#                 plot(daily.salinity[0,:,:30].mean(dim='date').T,daily.depth,'b',
#                 label='Mooring')#: '+str('{:5.2f}'.format(mean(cc['freshb'][:30].values)))+' mSv')
#                 plot(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').T,daily.depth,'b',label='',linewidth=2)
#             if kk==key16:
#                 plot(daily.salinity[0,:,-30:].mean(dim='date').T,daily.depth,'b',
#                 label='Mooring')#: '+str('{:5.2f}'.format(mean(cc['freshb'][-30:].values)))+' mSv')
#                 plot(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').T,daily.depth,'b',label='',linewidth=2)
#             # ylabel('depth [m]')
#             ylim([190,0])
#             xlim(25,35)
#
#             subplot(2,2,nn+1)
#             plot(squeeze(salvel[kk][:len(salprs[kk]),ctdind[kk]]),salprs[kk],color='orange',linewidth=2)
#             if kk==key14:
#                 moormonthly=daily['across track velocity'].where(daily.distance<distvec[1])[:,:,:30].mean(dim='date')
#             if kk==key16:
#                 moormonthly=daily['across track velocity'].where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date')
#             plot(moormonthly.T,daily.depth,'b',label='',linewidth=2)
#             gca().set_yticklabels('')
#             ylim([190,0])
#             xlim(-0.7,0.1)
#             nn+=2
#
#     subplot(221)
#     ylabel('August, 2014 \n \n depth [m]')
#     legend()
#     subplot(222)
#     subplot(223)
#     ylabel('July, 2016 \n \n depth [m]')
#     xlabel('salinity')
#     subplot(224)
#     xlabel('velocity [m/s]')
#
#
#     savefig('../figures/paperfigs/salvelprofcomp_coastal'+savename+'.pdf')
#
#
# plotsalvelprofs()
# for kk in [key14,key16]:
#         figure(figsize=(12,4))
#         subplot(121)
#         plot(squeeze(sal[kk][:len(salprs[kk]),0]),salprs[kk],color='orange',label='Shipboard: '+str('{:5.2f}'.format(shipfresh[kk]))+' mSv')
#         plot(squeeze(sal[kk][:len(salprs[kk]),ctdind[kk]]),salprs[kk],color='orange')
#         if kk==key14:
#             plot(daily.salinity[0,:,:30].mean(dim='date').T,daily.depth,'b',
#             label='Mooring: '+str('{:5.2f}'.format(mean(cc['freshb'][:30].values)))+' mSv')
#             plot(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').T,daily.depth,'b',label='')
#         if kk==key16:
#             plot(daily.salinity[0,:,-30:].mean(dim='date').T,daily.depth,'b',
#             label='Mooring: '+str('{:5.2f}'.format(mean(cc['freshb'][-30:].values)))+' mSv')
#             plot(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').T,daily.depth,'b',label='')
#         suptitle(kk[:4])
#         title('salinity profiles')
#         ylabel('depth [m]')
#         ylim([190,0])
#         legend()
#
#         subplot(122)
#         plot(squeeze(salvel[kk][:len(salprs[kk]),ctdind[kk]]*(sal[kk][:len(salprs[kk]),ctdind[kk]]-sref)/sref),salprs[kk],color='orange')
#         if kk==key14:
#             moormonthly=daily['across track velocity'].where(daily.distance<distvec[1])[:,:,:30].mean(dim='date')*(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date')-sref)/sref
#         if kk==key16:
#             moormonthly=daily['across track velocity'].where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date')*(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date')-sref)/sref
#         plot(moormonthly.T,daily.depth,'b',label='')
#         gca().set_yticklabels('')
#         title('freshwater transport profiles')
#         ylim([190,0])
#         savefig('../figures/compare_shipmoor/salfreshprofcomp_coastal'+str(kk[:4])+savename+'.pdf')




# def plt_moorgeo():
#     fs=11
#     fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8,6))
#     contgeo(ax1,key14)
#     ax1.set_xlim([-15,110])
#     ax1.set_ylim([800,0])
#     im=contmoor_ave(ax2,'2014-08-07','2014-09-07')
#     contgeo(ax3,key16)
#     im=contmoor_ave(ax4,'2016-06-28','2016-07-28')
#     ax1.set_title('Geostrophic velocity',fontsize=fs+2)
#     ax2.set_title('Monthly average from moorings',fontsize=fs+2)
#     ax1.set_ylabel('August, 2014 \n \n depth [m]',fontsize=fs)
#     ax3.set_ylabel('July, 2016 \n \n depth [m]',fontsize=fs)
#     mkbox(ax1,-12,middist[key14][minind_geo[key14]-1],150,5)
#     mkbox(ax2,-5,daily.distance[int(mean(minind[:30]))],150,5)
#     mkbox(ax3,-12,middist[key16][minind_geo[key16]-1],150,5)
#     mkbox(ax4,-5,daily.distance[int(mean(minind[-30:]))],150,5)
#     ax3.set_xlabel('distance [km]',fontsize=fs)
#     ax4.set_xlabel('distance [km]',fontsize=fs)
#     fig.subplots_adjust(hspace=0.1,wspace=0.1)
#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
#     fig.colorbar(im, cax=cbar_ax,label='across-track velocity [m/s]')
#
#     savefig('../figures/compare_shipmoor/Coastalcurr_geomoorcomp.pdf')
#     savefig('../figures/compare_shipmoor/Coastalcurr_geomoorcomp.png')
#
# plt_moorgeo()

XXXXXXXXXXX

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
# savefig('../figures/shipboard/ADCPmoorcomp_avevel.png',bbox_inches='tight')
#
# # compdistvel(nansum,'Integrated currents (to 300m)','integrated velocity [m/s]')
# # savefig('../figures/shipboard/ADCPmoorcomp_intvel.png',bbox_inches='tight')


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
# savefig('../figures/xport/cc_trans_fixedpos_salcomp.pdf',bbox_inches='tight')
#
# figure(figsize=(10,3))
# psf(cc['trans min vel'],'g',-1.5,0,'cc_trans',labit='min vel')
# psf(cc['trans min vel plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# savefig('../figures/xport/cc_trans_minvel_salcomp.pdf',bbox_inches='tight')
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
# savefig('../figures/xport/cc_trans_allcomp.pdf',bbox_inches='tight')
#
# figure(figsize=(10,3))
# psf(cc['trans min vel plus sal'],'r',-1.5,0,'cc_trans',labit='min vel plus sal')
# psf(cc['trans'],'b',-1.5,0,'cc_trans',labit='fixed pos')
# psf(cc['trans plus sal'],'k',-1.5,0,'cc_trans',labit='plus sal')
# legend()
# gca().set_yticks(arange(-1.5,0.3,0.5))
# ylabel('Transport [Sv]')
# savefig('../figures/xport/cc_trans_mvelsalpcomp.pdf',bbox_inches='tight')
#
#
# figure(figsize=(10,3))
# psf(mindxr,'purple',0,40,'position of surface velocity minimum')
# ylabel('')
# title('position of surface velocity minimum between currents [km]')
# savefig('../figures/xport/minpos_smooth.pdf',bbox_inches='tight')
#
#
# figure(figsize=(10,3))
# psf(maxdxr,'grey',-5,30,'position of coastal current core (max out to 30km)')
# ylabel('')
# title('position of coastal current core (max out to 30km) [km]')
# savefig('../figures/xport/maxpos_smooth.pdf',bbox_inches='tight')
