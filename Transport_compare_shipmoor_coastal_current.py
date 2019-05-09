################################################################################
################################################################################
##### Compare sections and transports between shipboard and moored measurements
################################################################################
################################################################################

from aux_funcs import *

#################################################################################
#### Load Mooring daily and set up transport calc framework
#################################################################################

newgrid=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1803extrap.pickle','rb'))

savename='_extrap'

bathf=interpolate.interp1d(bathdist,bathbath)
bathonmygrid=bathf(newgrid.distance)

daily=newgrid.copy()
for vv in newgrid:
    if vv[0]!='d':
        print(vv)
        for dd,adist in enumerate(daily.distance):
            daily[vv][dd,:,:]=daily[vv][dd,:,:].where(daily[vv][dd,:,:].depth<=bathonmygrid[dd])

bathdist=hstack((-20,bathdist))
bathbath=hstack((bathbath[0],bathbath))

import scipy.signal as sig
# Design the Buterworth filter
N  = 2    # Filter order
Wn = 0.02 # Cutoff frequency
B, A = sig.butter(N, Wn, output='ba')




def pwf(field,colo,nofilt=0,labit='',xr=daily):
    if nofilt==0:
        field.plot(alpha=0.5,color=colo,label='')
    plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit)


def psf(field,colo,ylim1,ylim2,tit,nofilt=0,xr=daily,labit=''):
    pwf(field,colo,nofilt,xr=xr,labit=labit)
    ylim([ylim1,ylim2])
    xlabel('')



# some experimentation with current position
minind=zeros(len(daily.date))
for tt,na in enumerate(daily.date):
    minind[tt]=where(max(daily['across track velocity'][2:15,0,tt])==daily['across track velocity'][2:15,0,tt])[0][0]+2

minind=[int(mm) for mm in minind]


ccvel=daily['across track velocity'].copy()
for tt,mm in enumerate(minind):
        ccvel[mm:,:,tt]=NaN

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

mid_dist=hstack((0,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,0))
mid_dist[1]=1.5 ## being strict -- only on one side of mooring
mid_dist[-1]=2.5 ## being strict -- only on one side of mooring
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

daily['xport']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat/1e3

##adding 10km width to CF1
mid_dist_plus=mid_dist.copy()
mid_dist_plus[0]=10
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))

daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

sref=34.9

srefb=34.8

srefc=35

cc={}
cc['trans']=daily.xport[:6,:,:].sum('depth').sum('distance')
cc['trans plus']=daily['xport plus'][:6,:,:].sum('depth').sum('distance')
# cc['trans sal def']=(daily['across track velocity'].where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
# cc['trans plus sal']=(daily['across track velocity'].where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3)[:6,:,:].sum('depth').sum('distance')
cc['trans min vel']=(ccvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
# cc['trans min vel plus sal']=(ccvel.where(daily.salinity<34)[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cc['fresh']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat[:6,:]*(daily.salinity[:6,:-1,:]-sref)/sref).sum('depth').sum('distance')
cc['fresh plus']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat_plus[:6,:]*(daily.salinity[:6,:-1,:]-sref)/sref).sum('depth').sum('distance')

cc['fresh 34.8']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat[:6,:]*(daily.salinity[:6,:-1,:]-srefb)/srefb).sum('depth').sum('distance')
cc['fresh plus 34.8']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat_plus[:6,:]*(daily.salinity[:6,:-1,:]-srefb)/srefb).sum('depth').sum('distance')

cc['fresh 35.0']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat[:6,:]*(daily.salinity[:6,:-1,:]-srefc)/srefc).sum('depth').sum('distance')
cc['fresh plus 35.0']=(daily['across track velocity'][:6,:-1,:]*depthdiffmat[:6,:]*middistmat_plus[:6,:]*(daily.salinity[:6,:-1,:]-srefc)/srefc).sum('depth').sum('distance')

cc['trans plus'].std()

cc['fresh'].mean()
cc['fresh plus'].mean()

date1=datetime.datetime(2014,8,15)
daily.date[0].values

sig.filtfilt(B,A, field)

date1=datetime.datetime(2014,8,15)

## Print max and min (freshwater) transports
middate=datetime.datetime(2015,5,1,0)
for key in cc:
    print(key)
    print(max(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(date1,middate))))))
    print(daily.date[argmax(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(date1,middate)))))].values)


lastdate=datetime.datetime(2016,7,30)
for key in cc:
    print(key)
    print(max(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(middate,lastdate))))))
    print(daily.date.sel(date=slice(middate,lastdate))[argmax(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(middate,lastdate)))))].values)

mdate1=datetime.datetime(2015,10,1)
mdate2=datetime.datetime(2016,9,1)
for key in cc:
    print(key)
    print(min(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(mdate1,mdate2))))))
    print(daily.date.sel(date=slice(mdate1,mdate2))[argmin(abs(sig.filtfilt(B,A,cc[key].sel(date=slice(mdate1,mdate2)))))].values)

for key in cc:
    print(key)
    print(std(cc[key]).values)


## Print (freshwater) transport for first two weeks of August

date2=datetime.datetime(2014,9,1)
for key in cc:
        print(key)
        print(cc[key].sel(date=slice(date1,date2)).mean().values)


figure(figsize=(10,3))
psf(cc['trans'],'b',-2.5,0.5,'cc_trans',labit='moorings only (~0.5 Sv)')
psf(cc['trans plus'],'g',-2.5,0.5,'cc_trans',labit='adding 10km width (~1 Sv)')
legend()
title('Coastal current transport (plus unmeasured shelf)')
savefig('../figures/compare_shipmoor/Coastalcurr_addcoast.pdf',bbox_inches='tight')

#################################################################################
# Load and calc trans from geostrophic velocity from CTD
#################################################################################

[prs_ctd,dist,middist,lon,lat,sal,tmp,pden,geoshear,geovel]=pickle.load(open('../pickles/Shipboard/CTD_1416geo.pickle','rb'))

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

transgeo_sub

################################################################################
##### Load ADCP and calc trans
################################################################################

adcp_dist=array(pd.DataFrame.from_csv(datadir+'Shipboard/adcp_distances.dat').index)

with open(datadir+'Shipboard/adcp_distances.dat', 'r') as f:
    reader = csv.reader(f)
    adcp_dist = list(reader)

adcp_dist=[float(dd[0]) for dd in adcp_dist]

adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')
adcp_16=io.loadmat(datadir+'Shipboard/ar07_2016/ar07_2016_vm_adcp_ilebras.mat')


# plot(adcp_14['lon'],adcp_14['lat'],'ko');
# plot(adcp_16['lon'],adcp_16['lat'],'ro');
# plot(360+CFlon,CFlat,'yo');
# ylim([60,60.2])
# xlim([316.75,318])

def acrosstrackname(adic):
    adic['across track velocity']=adic['v_dt']*sin(theta)+adic['u_dt']*cos(theta)

    return adic

adcp_14=acrosstrackname(adcp_14)
adcp_16=acrosstrackname(adcp_16)

# adcp transport calcs
mid_dist_adcp=hstack((diff(adcp_dist)[0]/2,(diff(adcp_dist)[:-1]+diff(adcp_dist)[1:])/2,diff(adcp_dist)[-1]/2))
middistmat_a14=tile(mid_dist_adcp,[len(adcp_14['depth'][:,0]),1])
middistmat_a16=tile(mid_dist_adcp[:-3],[len(adcp_16['depth'][:,0]),1])
shape(middistmat_a16)
depthdiffmat_a14=tile(diff(hstack((0,adcp_14['depth'][:,0]))),[len(adcp_dist),1]).T
depthdiffmat_a16=tile(diff(hstack((0,adcp_16['depth'][:,0]))),[len(adcp_dist[:-3]),1]).T
shape(depthdiffmat_a16)

minind14=where(adcp_14['across track velocity'][0,-20:]==max(adcp_14['across track velocity'][0,-20:]))[0][0]+len(adcp_14['across track velocity'][0,:])-20
minind16=where(adcp_16['across track velocity'][0,-10:-7]==max(adcp_16['across track velocity'][0,-10:-7]))[0][0]+len(adcp_16['across track velocity'][0,:])-10

trans14adcp_less=nansum(adcp_14['across track velocity'][:,minind14:-3]*depthdiffmat_a14[:,minind14:-3]*middistmat_a14[:,minind14:-3])/1e3
trans14adcp=nansum(adcp_14['across track velocity'][:,minind14:]*depthdiffmat_a14[:,minind14:]*middistmat_a14[:,minind14:])/1e3
trans16adcp=nansum(adcp_16['across track velocity'][:,minind16:]*depthdiffmat_a16[:,minind16:]*middistmat_a16[:,minind16:])/1e3

moorind=zeros(len(distvec))
for dd in range(len(distvec)):
    moorind[dd]=argmin(abs(adcp_dist-distvec[dd]))

moorind=[int(dd) for dd in moorind][::-1]

## transport between mooring locations
trans14adcp_moorloc=nansum(adcp_14['across track velocity'][:,moorind[-2]:moorind[-1]]*depthdiffmat_a14[:,moorind[-2]:moorind[-1]]*middistmat_a14[:,moorind[-2]:moorind[-1]])/1e3
trans16adcp_moorloc=nansum(adcp_16['across track velocity'][:,moorind[-2]:moorind[-1]]*depthdiffmat_a16[:,moorind[-2]:moorind[-1]]*middistmat_a16[:,moorind[-2]:moorind[-1]])/1e3

## transport from subsampling at mooring locations
adcp_dist_sub=[adcp_dist[dd] for dd in moorind]
mid_dist_adcp_sub=hstack((diff(adcp_dist_sub)[0]/2,(diff(adcp_dist_sub)[:-1]+diff(adcp_dist_sub)[1:])/2,diff(adcp_dist_sub)[-1]/2))

# subsample adcp dailya at mooring locations
adcpvel_14_sub=[adcp_14['across track velocity'][:,dd] for dd in moorind[-2:]]
adcpvel_16_sub=[adcp_16['across track velocity'][:,dd] for dd in moorind[-2:]]
#calculate transport from CF1 and CF2 positions, using nominal width of 12km each
depthdiffmat_a14_sub=tile(diff(hstack((0,adcp_14['depth'][:,0]))),[2,1])
depthdiffmat_a16_sub=tile(diff(hstack((0,adcp_16['depth'][:,0]))),[2,1])

trans14adcp_sub=nansum(adcpvel_14_sub*depthdiffmat_a14_sub*12)/1e3
trans16adcp_sub=nansum(adcpvel_16_sub*depthdiffmat_a16_sub*12)/1e3

## read out the values
trans14adcp
trans14adcp_less
trans14adcp_moorloc
trans14adcp_sub

trans16adcp
trans16adcp_moorloc
trans16adcp_sub


cc['trans'][0].values
cc['trans plus'][0].values
cc['trans min vel'][0].values


mean(cc['trans'][:30]).values
mean(cc['trans plus'][:30]).values
mean(cc['trans min vel'][:30]).values

cc['trans'][-1].values
cc['trans plus'][-1].values
cc['trans min vel'][-1].values

mean(cc['trans'][-30:]).values
mean(cc['trans plus'][-30:]).values
mean(cc['trans min vel'][-30:]).values


#################################################################################
# Combine shipboard ADCP and CTD to get freshwater transport
#################################################################################

## Note, for 2014 using offset salinity because it goes farther onto shelf
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


adcp_14['depth'][:,0]

#interpolate adcp onto salinity distances and depths
def interptosal(velfield,veldistvec,velprs,salfield,newdist):
    [prslen,distlen]=shape(velfield)
    # newdist=saldist#[saldist>=min(veldistvec)]
    newprs=prs_ctd[prs_ctd<=max(velprs)]
    salvel1=zeros((prslen,len(newdist)))
    salvel2=zeros((len(newprs),len(newdist)))
    #interpolate to distances
    for pp in range(prslen):
        distfunc=interpolate.interp1d(hstack((min(newdist),veldistvec)),hstack((velfield[pp,0],velfield[pp,:])))
        salvel1[pp,:]=distfunc(newdist)
    #interpolate to depth
    for dd in range(len(newdist)):
        prsfunc=interpolate.interp1d(hstack((0,velprs[~isnan(salvel1[:,dd])],180)),hstack((salvel1[0,dd],salvel1[~isnan(salvel1[:,dd]),dd],salvel1[~isnan(salvel1[:,dd]),dd][-1])))
        salvel2[:,dd]=prsfunc(newprs)

    return salvel2,newprs

salvel={}
salprs={}
salvel[key14],salprs[key14]=interptosal(adcp_14['across track velocity'][adcp_14['depth'][:,0]<180,:],adcp_dist,adcp_14['depth'][:,0][adcp_14['depth'][:,0]<180],sal[key14],dist[key14])
salvel[key16],salprs[key16]=interptosal(adcp_16['across track velocity'][adcp_16['depth'][:,0]<180,:],adcp_dist[:-3],adcp_16['depth'][:,0][adcp_16['depth'][:,0]<180],sal[key16],dist[key16])

plot(adcp_dist,adcp_14['across track velocity'][0,:])
xlim([-20,50])
plot(dist[key14][:17],salvel[key14][0,:17])

contourf(dist[key14],-salprs[key14],salvel[key14])

dist[key14][8]

dist[key16][8]

sref=34.9

# freshwater transport calcs
middist_fresh={}
middistmat_fresh={}
depthdiffmat_fresh={}
minind_fresh={}
shipfresh={}
shipfresh_moorloc={}
transtest={}
moordistind={}
for kk in [key14,key16]:
    middist_fresh[kk]=hstack((diff(dist[kk])[0]/2,(diff(dist[kk])[:-1]+diff(dist[kk])[1:])/2,diff(dist[kk])[-1]/2))
    middistmat_fresh[kk]=tile(middist_fresh[kk],[len(salprs[kk]),1])
    depthdiffmat_fresh[kk]=tile(diff(hstack((0,salprs[kk]))),[len(dist[kk]),1]).T
    minind_fresh[kk]=where(salvel[kk][0,2:17]==max(salvel[kk][0,2:17]))[0][0]+3
    transtest[kk]=nansum(salvel[kk][:,:minind_fresh[kk]]*depthdiffmat_fresh[kk][:,:minind_fresh[kk]]*middistmat_fresh[kk][:,:minind_fresh[kk]])/1e3
    shipfresh[kk]=nansum(salvel[kk][:,:minind_fresh[kk]]*(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]-sref)/sref*depthdiffmat_fresh[kk][:,:minind_fresh[kk]]*middistmat_fresh[kk][:,:minind_fresh[kk]])
    moordistind[kk]=where((dist[kk]>=0)&(dist[kk]<=distvec[1]))[0]
    shipfresh_moorloc[kk]=nansum(salvel[kk][:,moordistind[kk]]*(sal[kk][:len(salprs[kk]),moordistind[kk]]-sref)/sref*depthdiffmat_fresh[kk][:,moordistind[kk]]*middistmat_fresh[kk][:,moordistind[kk]])

datechoose={}
datechoose[key14]=0
datechoose[key16]=-1


for kk in [key14,key16]:
        figure(figsize=(12,4))
        subplot(121)
        plot(squeeze(sal[kk][:len(salprs[kk]),0]),salprs[kk],color='orange',label='CTD within ADCP coastal current: '+str('{:5.2f}'.format(shipfresh[kk]))+' mSv')
        plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk][0]]),salprs[kk],'r--',label='CTD within mooring range: '+str('{:5.2f}'.format(shipfresh_moorloc[kk]))+' mSv')
        plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]),salprs[kk],color='orange')
        plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk]]),salprs[kk],'r--')
        plot(daily.salinity.isel(date=datechoose[kk],distance=0).T,daily.depth,'k',
        label='Moorings, daily: '+str('{:5.2f}'.format(cc['fresh'][datechoose[kk]].values))+' mSv ('+str('{:5.2f}'.format(cc['fresh plus'][datechoose[kk]].values))+' mSv)')
        plot(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1]).T,daily.depth,'k',label='')
        if kk==key14:
            plot(daily.salinity[0,:,:30].mean(dim='date').T,daily.depth,'b',
            label='Moorings, monthly: '+str('{:5.2f}'.format(mean(cc['fresh'][:30].values)))+' mSv ('+str('{:5.2f}'.format(mean(cc['fresh plus'][:30].values)))+' mSv)')
            plot(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').T,daily.depth,'b',label='')
        if kk==key16:
            plot(daily.salinity[0,:,-30:].mean(dim='date').T,daily.depth,'b',
            label='Moorings, monthly: '+str('{:5.2f}'.format(mean(cc['fresh'][-30:].values)))+' mSv ('+str('{:5.2f}'.format(mean(cc['fresh plus'][:30].values)))+' mSv)')
            plot(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').T,daily.depth,'b',label='')
        suptitle(kk[:4])
        title('salinity profiles')
        ylabel('depth [m]')
        ylim([190,0])
        legend()
        subplot(122)
        plot(squeeze(salvel[kk][:len(salprs[kk]),:minind_fresh[kk]]*(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]-sref)/sref),salprs[kk],color='orange')
        plot(squeeze(salvel[kk][:len(salprs[kk]),moordistind[kk]]*(sal[kk][:len(salprs[kk]),moordistind[kk]]-sref)/sref),salprs[kk],'r--')
        moordaily=daily['across track velocity'].isel(date=datechoose[kk]).where(daily.distance<distvec[1])*(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1])-sref)/sref
        plot(moordaily.T,daily.depth,'k',label='')
        if kk==key14:
            moormonthly=daily['across track velocity'].where(daily.distance<distvec[1])[:,:,:30].mean(dim='date')*(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date')-sref)/sref
        if kk==key16:
            moormonthly=daily['across track velocity'].where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date')*(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date')-sref)/sref
        plot(moormonthly.T,daily.depth,'b',label='')
        gca().set_yticklabels('')
        title('freshwater transport profiles')
        ylim([190,0])
        savefig('../figures/compare_shipmoor/salfreshprofcomp_coastal'+str(kk[:4])+savename+'.pdf')

def salvelcomp_paper():
    figure(figsize=(10,4))
    for kk in [key14,key16]:
            subplot(121)
            if kk==key14:
                kline='-'
            elif kk==key16:
                kline='--'
            plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'r',linestyle=kline,label='CTD within ADCP coastal current: '+str('{:5.2f}'.format(shipfresh[kk]))+' mSv')
            plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk]]).mean(axis=1),salprs[kk],color='orange',linestyle=kline,label='CTD within mooring range: '+str('{:5.2f}'.format(shipfresh_moorloc[kk]))+' mSv')
            plot(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k',linestyle=kline,
            label='Moorings, daily: '+str('{:5.2f}'.format(cc['fresh'][datechoose[kk]].values))+' mSv ('+str('{:5.2f}'.format(cc['fresh plus'][datechoose[kk]].values))+' mSv)')
            if kk==key14:
                plot(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline,
                label='Moorings, monthly: '+str('{:5.2f}'.format(mean(cc['fresh'][:30].values)))+' mSv ('+str('{:5.2f}'.format(mean(cc['fresh plus'][:30].values)))+' mSv)')
            if kk==key16:
                plot(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline,
                label='Moorings, monthly: '+str('{:5.2f}'.format(mean(cc['fresh'][-30:].values)))+' mSv ('+str('{:5.2f}'.format(mean(cc['fresh plus'][:30].values)))+' mSv)')
            ylabel('depth [m]')
            ylim([190,0])
            xlabel('mean coastal current salinity')
            legend()
            subplot(122)
            plot(mean(salvel[kk][:len(salprs[kk]),:minind_fresh[kk]],axis=1),salprs[kk],'r',linestyle=kline,label='CTD within ADCP coastal current: '+str('{:5.2f}'.format(shipfresh[kk]))+' mSv')
            plot(mean(salvel[kk][:len(salprs[kk]),moordistind[kk]],axis=1),salprs[kk],color='orange',linestyle=kline,label='CTD within mooring range: '+str('{:5.2f}'.format(shipfresh_moorloc[kk]))+' mSv')
            plot(daily['across track velocity'].isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k',linestyle=kline,
            label='Moorings, daily: '+str('{:5.2f}'.format(cc['fresh'][datechoose[kk]].values))+' mSv ('+str('{:5.2f}'.format(cc['fresh plus'][datechoose[kk]].values))+' mSv)')
            if kk==key14:
                plot(daily['across track velocity'].where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline,
                label='Moorings, monthly: '+str('{:5.2f}'.format(mean(cc['fresh'][:30].values)))+' mSv ('+str('{:5.2f}'.format(mean(cc['fresh plus'][:30].values)))+' mSv)')
            if kk==key16:
                plot(daily['across track velocity'].where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline,
                label='Moorings, monthly: '+str('{:5.2f}'.format(mean(cc['fresh'][-30:].values)))+' mSv ('+str('{:5.2f}'.format(mean(cc['fresh plus'][:30].values)))+' mSv)')
            gca().set_yticklabels('')
            ylim([190,0])
            tight_layout()
            xlabel('mean coastal current velocity [m/s]')

salvelcomp_paper()

def salvelcomp_paper():
            figure(figsize=(9,5))
            subplot(121)
            kk=key14
            plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'-',color='grey',label='Deployment, 2014')
            plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'--',color='grey',label='Recovery, 2016')
            plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'r',label='Full Shipboard')
            kk=key16
            plot(squeeze(sal[kk][:len(salprs[kk]),:minind_fresh[kk]]).mean(axis=1),salprs[kk],'r--')
            kk=key14
            plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk]]).mean(axis=1),salprs[kk],color='orange',label='Shipboard within mooring range')
            kk=key16
            plot(squeeze(sal[kk][:len(salprs[kk]),moordistind[kk]]).mean(axis=1),salprs[kk],color='orange',linestyle='--',label='')
            kk=key14
            plot(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k',label='Mooring daily average')
            kk=key16
            plot(daily.salinity.isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k--',label='')
            plot(daily.salinity.where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').mean(dim='distance').T,daily.depth,'b',label='Mooring monthly average')
            plot(daily.salinity.where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').mean(dim='distance').T,daily.depth,'b--',label='')
            ylabel('depth [m]')
            ylim([190,0])
            xlabel('mean coastal current salinity')
            legend()
            subplot(122)
            for kk in [key14,key16]:
                if kk==key14:
                    kline='-'
                elif kk==key16:
                    kline='--'
                plot(mean(salvel[kk][:len(salprs[kk]),:minind_fresh[kk]],axis=1),salprs[kk],'r',linestyle=kline)
                plot(mean(salvel[kk][:len(salprs[kk]),moordistind[kk]],axis=1),salprs[kk],color='orange',linestyle=kline)
                plot(daily['across track velocity'].isel(date=datechoose[kk]).where(daily.distance<distvec[1]).mean(dim='distance').T,daily.depth,'k',linestyle=kline,)
                if kk==key14:
                    plot(daily['across track velocity'].where(daily.distance<distvec[1])[:,:,:30].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline)
                if kk==key16:
                    plot(daily['across track velocity'].where(daily.distance<distvec[1])[:,:,-30:].mean(dim='date').mean(dim='distance').T,daily.depth,'b',linestyle=kline)
                gca().set_yticklabels('')
                ylim([190,0])
                tight_layout()
                xlabel('mean coastal current velocity [m/s]')
            savefig('../figures/paperfigs/salvelcomp_coastal.pdf',bbox_inches='tight')


salvelcomp_paper()

psf(cc['fresh'],ccol,0,100,'Coastal current freshwater flux')
psf(cc['fresh plus'],'grey',0,100,'Coastal current freshwater flux')

psf(daily.salinity.isel(depth=30,distance=1),'b',32,34.5,'Salinity at CF1 (60m recon)')

psf(daily.salinity.isel(depth=75,distance=1),'b',33,34.5,'Salinity at CF1 (150m)')

psf(daily.salinity.isel(depth=30,distance=5),'b',32,34.5,'Salinity at CF2 (60m)')

psf(daily.salinity.isel(depth=75,distance=5),'b',33,34.5,'Salinity at CF2 (150m)')

transtest
shipfresh
shipfresh_moorloc

cc['fresh'][0]
cc['fresh'][-1]
mean(cc['fresh'][:30])
mean(cc['fresh'][-30:])


#################################################################################
#### Plots comparing shipboard and mooring
#################################################################################

# minimal framework for adcp contouring
def contadcp(axspec,distveca,adic):
    axspec.contourf(distveca,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(distveca,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]

#same for mooring snap
def contmoor_pnt(axspec,datechoice):
    im=axspec.contourf(daily.distance,daily.depth,daily['across track velocity'][:,:,datechoice].T,arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(daily.distance,daily.depth,daily['across track velocity'][:,:,datechoice].T,[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return im


#mooring version for longer mean!

def contmoor_ave(axspec,date1,date2):
    im=axspec.contourf(daily.distance,daily.depth,daily['across track velocity'].sel(date=slice(date1,date2)).mean(dim='date').T,arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(daily.distance,daily.depth,daily['across track velocity'].sel(date=slice(date1,date2)).mean(dim='date').T,[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    return im

#adapted for geostrophic vel
def contgeo(axspec,key):
    axspec.contourf(dist[key][:-1],prs_ctd,geovel[key],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r,extend='both');
    axspec.contour(dist[key][:-1],prs_ctd,geovel[key],[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]


def mkbox(axc,x1,x2,y1,y2):
    axc.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],color='purple',linewidth=3)


def plt_mooradcp():
    fs=11
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8,6))
    contadcp(ax1,adcp_dist,adcp_14)
    ax1.set_xlim([-15,110])
    ax1.set_ylim([800,0])
    im=contmoor_pnt(ax2,0)
    contadcp(ax3,adcp_dist[:-3],adcp_16)
    im=contmoor_pnt(ax4,-1)
    ax1.set_title('Vessel mounted ADCP',fontsize=fs+2)
    ax2.set_title('Daily snap from moorings',fontsize=fs+2)
    ax1.set_ylabel('August, 2014 \n \n depth [m]',fontsize=fs)
    ax3.set_ylabel('July, 2016 \n \n depth [m]',fontsize=fs)
    mkbox(ax1,-12,adcp_dist[minind14],150,5)
    mkbox(ax2,-5,daily.distance[minind[0]],150,5)
    mkbox(ax3,-5,adcp_dist[minind16],150,5)
    mkbox(ax4,-5,daily.distance[minind[-1]],150,5)
    ax3.set_xlabel('distance [km]',fontsize=fs)
    ax4.set_xlabel('distance [km]',fontsize=fs)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax,label='across-track velocity [m/s]')

    savefig('../figures/compare_shipmoor/Coastalcurr_adcpmoorcomp.pdf')
    savefig('../figures/compare_shipmoor/Coastalcurr_adcpmoorcomp.png')


plt_mooradcp()

geovel[key14][1,:minind_geo[key14]]

def plt_moorgeo():
    fs=11
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8,6))
    contgeo(ax1,key14)
    ax1.set_xlim([-15,110])
    ax1.set_ylim([800,0])
    im=contmoor_ave(ax2,'2014-08-07','2014-09-07')
    contgeo(ax3,key16)
    im=contmoor_ave(ax4,'2016-06-28','2016-07-28')
    ax1.set_title('Geostrophic velocity',fontsize=fs+2)
    ax2.set_title('Monthly average from moorings',fontsize=fs+2)
    ax1.set_ylabel('August, 2014 \n \n depth [m]',fontsize=fs)
    ax3.set_ylabel('July, 2016 \n \n depth [m]',fontsize=fs)
    mkbox(ax1,-12,middist[key14][minind_geo[key14]-1],150,5)
    mkbox(ax2,-5,daily.distance[int(mean(minind[:30]))],150,5)
    mkbox(ax3,-12,middist[key16][minind_geo[key16]-1],150,5)
    mkbox(ax4,-5,daily.distance[int(mean(minind[-30:]))],150,5)
    ax3.set_xlabel('distance [km]',fontsize=fs)
    ax4.set_xlabel('distance [km]',fontsize=fs)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax,label='across-track velocity [m/s]')

    savefig('../figures/compare_shipmoor/Coastalcurr_geomoorcomp.pdf')
    savefig('../figures/compare_shipmoor/Coastalcurr_geomoorcomp.png')

plt_moorgeo()

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
