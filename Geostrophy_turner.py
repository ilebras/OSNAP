from aux_funcs import *

datin=pickle.load(open('../pickles/xarray/CF_xarray_notid_SAtheta.pickle','rb'))

fcor=gsw.f(60)

geoshear=(-diff(datin['potential density'],axis=0).T*9.8/fcor/1028/diff(datin.distance).T/1e3).T

plot(mean(geoshear,axis=2).T);

geovel=nancumsum(geoshear[:,::-1,:],axis=1)[:,::-1,:]*2

prsvec=range(0,2113,2)

SA=datin['salinity'].copy()
for ii in range(8):
    SA[ii,:,:]=gsw.SA_from_SP(datin['salinity'][ii,:,:].T,prsvec,CFlon[ii],CFlat[ii]).T

CT=gsw.CT_from_pt(SA,datin['temperature'])

SAint=SA[:-1,:,:]+diff(SA,axis=0)
CTint=CT[:-1,:,:]+diff(CT,axis=0)
alphaint=CTint.copy()
betaint=SAint.copy()
for ii in range(7):
    alphaint[ii,:,:]=gsw.alpha(SAint[ii,:,:].T,CTint[ii,:,:].T,prsvec).T
    betaint[ii,:,:]=gsw.beta(SAint[ii,:,:].T,CTint[ii,:,:].T,prsvec).T

datin

avelint=datin['across track velocity'][:-1,:,:]+diff(datin['across track velocity'],axis=0)/2

alphaT=-alphaint*diff(datin['temperature'],axis=0)
betaS=betaint*diff(datin['salinity'],axis=0)

dendiff=diff(datin['potential density'],axis=0)

R=alphaT/betaS

turner=arctan2((alphaT-betaS)*dendiff/abs(dendiff),(alphaT+betaS)*dendiff/abs(dendiff))*180/pi

turner_noden=arctan2((alphaT-betaS),(alphaT+betaS))*180/pi

veldiff=avelint-geovel

# Then, establish min in abs vel to differentiate EGCC
minind=zeros(len(datin.date))
for tt,na in enumerate(datin.date):
    minind[tt]=where(max(datin['across track velocity'][:5,0,tt])==datin['across track velocity'][:5,0,tt])[0][0]+1

minind=[int(mm) for mm in minind]



dat=xr.Dataset({'geostrophic velocity': (['distance', 'depth', 'date'],  geovel),
                'absolute velocity': (['distance', 'depth', 'date'],  avelint),
                'velocity difference': (['distance', 'depth', 'date'],  veldiff),},
                # 'turner angle': (['distance', 'depth', 'date'],  turner),
                # 'turner angle noden': (['distance', 'depth', 'date'],  turner_noden),
                # 'density ratio': (['distance', 'depth', 'date'],  R),
                # 'alpha DeltaT': (['distance', 'depth', 'date'],  alphaT),
                # 'beta': (['distance', 'depth', 'date'],  betaint),
                # 'alpha': (['distance', 'depth', 'date'],  alphaint),
                # 'DeltaS': (['distance', 'depth', 'date'], diff(datin['salinity'],axis=0) ),
                # 'DeltaT': (['distance', 'depth', 'date'], diff(datin['temperature'],axis=0) ),
                # 'beta DeltaS': (['distance', 'depth', 'date'],  betaS),
                # 'density gradient': (['distance', 'depth', 'date'],  dendiff)},
                coords={'distance': datin.distance[:-1]+diff(datin.distance)/2,
                        'depth': datin.depth,
                        'date': datin.date})

cc_absvel=dat['absolute velocity'].copy()
for tt,mm in enumerate(minind):
        cc_absvel[mm:,:,tt]=NaN

cc_geovel=dat['geostrophic velocity'].copy()
for tt,mm in enumerate(minind):
        cc_geovel[mm:,:,tt]=NaN

eg_absvel=dat['absolute velocity'].copy()
for tt,mm in enumerate(minind):
        eg_absvel[:mm,:,tt]=NaN

eg_geovel=dat['geostrophic velocity'].copy()
for tt,mm in enumerate(minind):
        eg_geovel[:mm,:,tt]=NaN

mid_dist=hstack((0,(diff(dat.distance)[:-1]+diff(dat.distance)[1:])/2,0))
middistmat=transpose((tile(mid_dist,[len(dat.depth)-1,len(dat.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(dat.depth),[len(dat.distance),len(dat.date),1])),(0,2,1))


cc={}
cc['abs trans']=(cc_absvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
cc['geo trans']=(cc_geovel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

eg={}
eg['abs trans']=(eg_absvel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')
eg['geo trans']=(eg_geovel[:,:-1,:]*depthdiffmat*middistmat/1e3).sum('depth').sum('distance')

cc['abs trans'].plot()


def velcomplot(dic,colc,tit):
    figure(figsize=(12,3))
    dic['abs trans'].plot(color=colc,linestyle='-',label='absolute transport',)
    dic['geo trans'].plot(color=colc,linestyle='--',label='gestrophic transport',)
    legend()
    title(tit)

velcomplot(eg,egicol,'Slope Current')
ylim([-2.5,0])

velcomplot(cc,ccol,'Coastal Current')
ylim([-2.5,0])


def onecont(field,tit,vrange,coloor,hlevs,nomoorlines=0):
    ax1=contourf(dat.distance,dat.depth,field,vrange,cmap=coloor,extend='both')
    ax2=contour(dat.distance,dat.depth,field,levels=hlevs,colors='k')
    clabel(ax2)
    fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance (km)')
    ylabel('depth (m)')
    xlim([-5,100])
    ylim([2200,0])
    title(tit)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return ax1,ax2

def meanplot(field):
    figure()
    ax1,ax2=onecont(dat[univec[field][0]].mean(dim='date').T,'Mean '+univec[field][0],univec[field][1],univec[field][2],univec[field][3])
    colorbar(ax1,label=univec[field][4])
    savefig('../figures/sections/fullmean_'+field+'_fulldpth.png',bbox_inches='tight')
    savefig('../figures/sections/fullmean_'+field+'_fulldpth.pdf',bbox_inches='tight')


def stdplot(field,limoverride):
    figure()
    ax1,ax2=onecont(dat[univec[field][0]].std(dim='date').T,'Standard deviation of '+univec[field][0],limoverride,cm.Purples,limoverride[::5])
    colorbar(ax1,label=univec[field][4])
    savefig('../figures/sections/fullstd_'+field+'_fulldpth.png',bbox_inches='tight')
    savefig('../figures/sections/fullstd_'+field+'_fulldpth.pdf',bbox_inches='tight')

titvec=['August - October','November - January','February - April','May - July']

def fourplot(field,secondyear=0):
    fig=figure(figsize=(8,6))
    for num in range(4):
        if secondyear==1:
            datenum=num+4
            datename='second'
        elif secondyear==0:
            datenum=num
            datename='first'
        else:
            datename='both'

        subplot(2,2,num+1)
        if datename=='both':
            elfieldo=(dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,num]+dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,num+4]).T/2
        else:
            elfieldo=dat[univec[field][0]].resample('3M',how='mean',dim='date')[:,:,datenum].T
        ax1,ax2=onecont(elfieldo,titvec[num],univec[field][1],univec[field][2],
        univec[field][3],nomoorlines=1)

    fig.subplots_adjust(hspace=0.15,wspace=0.1)
    subplot(221)
    gca().set_xticklabels('')
    xlabel('')
    subplot(222)
    gca().set_xticklabels('')
    gca().set_yticklabels('')
    xlabel('')
    ylabel('')
    subplot(224)
    gca().set_yticklabels('')
    ylabel('')

    fig.subplots_adjust(right=0.8)
    cbaxes = fig.add_axes([0.85, 0.15, 0.025, 0.6])
    cbar=colorbar(ax1, cax = cbaxes,ticks=univec[field][3],label=univec[field][4])
    suptitle(univec[field][0],fontsize=20)
    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year.png')
    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year.pdf')
    # figure(figsize=(12,6))
    for nn in range(1,5):
        subplot(2,2,nn)
        ylim([400,0])
        clabel(ax2)

    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year_zoom.png')
    savefig('../figures/sections/seasonal_'+field+'_'+datename+'year_zoom.pdf')

def plot12all(field):
    fourplot(field)
    fourplot(field,1)
    fourplot(field,4)


def seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='both',tu=0):
    if ychoice=='both':
        sind=0
        fin=8
    elif ychoice=='first':
        sind=0
        fin=4
    elif ychoice=='second':
        sind=4
        fin=8
    figure()
    for ii in range(sind,fin):
        if atype=='1depth':
            dathere=dat[field].resample('3M',dim='date',how='mean')[:,int(dpchoice/2),ii]
            disthere=dat.distance
        elif atype=='diff':
            dathere=diff(dat[field].resample('3M',dim='date',how='mean')[:,int(dpchoice/2),ii],axis=0)[1:-1]
            disthere=dat.distance[1:-2]+diff(dat.distance[1:-1])/2
        elif atype=='mean':
            dathere=mean(dat[field].resample('3M',dim='date',how='mean')[:,:int(dpchoice/2),ii],axis=1)
            disthere=dat.distance
        if ii>=4:
            lsty='--'
            lab=''
        else:
            lsty='-'
            lab=titvec[ii]
        if ychoice=='second':
            lab=titvec[ii-4]
        plot(disthere,dathere,label=lab,linestyle=lsty,marker='o',color=colvec[mod(ii,4)])
    plot(-5,mean(dathere),'k-',label='First year, 2014-2015')
    plot(-5,mean(dathere),'k--',label='Second year, 2015-2016')
    ylabel(ylab)
    xlabel('distance (km)')
    title(field+tit)
    legend(loc=(1.05,0.2))
    axhline(0,color='k',alpha=0.5)
    xlim([-5,100])
    if tu==1:
        [axhline(ii,color='k',alpha=0.5) for ii in range(-90,100,45)]
        yticks(range(-90,100,45))
    savefig('../figures/seasonal/seasonline_'+savename+'_'+ychoice+'.png',bbox_inches='tight')


colvec=['#e41a1c','#ff7f00','#377eb8','#4daf4a']

def seasonline_all(field,dpchoice,ylab,tit,atype,savename,ychoice='both',tu=0):
    seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='both',tu=tu)
    seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='first',tu=tu)
    seasonline(field,dpchoice,ylab,tit,atype,savename,ychoice='second',tu=tu)


def seasonline_comp(field,dpchoice,ylab,tit,savename,add,ychoice='first',tu=0):
    if ychoice=='first':
        sind=0
        fin=4
    elif ychoice=='second':
        sind=4
        fin=8
    figure()
    for ii in range(sind,fin):
        dathere=mean(dat[field].resample('3M',dim='date',how='mean')[:,:int(dpchoice/2),ii],axis=1)
        disthere=dat.distance
        datcomp=mean(dat['absolute velocity'].resample('3M',dim='date',how='mean')[:,:int(dpchoice/2),ii],axis=1)
        distcomp=dat.distance
        lab=titvec[mod(ii,4)]
        lsty='-'
        plot(disthere,dathere-add,label=lab,linestyle=lsty,marker='o',color=colvec[mod(ii,4)])
        lsty='--'
        plot(distcomp,datcomp,label='',linestyle=lsty,marker='d',color=colvec[mod(ii,4)])
    plot(-5,mean(dathere),'k-',label='Geostrophic velocity')
    plot(-5,mean(dathere),'k--',label='Absolute velocity')
    ylabel(ylab)
    xlabel('distance (km)')
    title(tit)
    legend(loc=(1.05,0.2))
    axhline(0,color='k',alpha=0.5)
    xlim([-5,100])
    if tu==1:
        [axhline(ii,color='k',alpha=0.5) for ii in range(-90,100,45)]
        yticks(range(-90,100,45))
    savefig('../figures/seasonal/seasonline_both_'+savename+'_'+ychoice+'.png',bbox_inches='tight')

seasonline_comp('geostrophic velocity',200,'[m/s]','Mean over top 200m, first year','vel200_pos02',0.2,ychoice='first')

seasonline_comp('geostrophic velocity',200,'[m/s]','Mean over top 200m, second year','vel200_pos02',0.2,ychoice='second')


dat

meanplot('geostrophic velocity')
stdplot('geostrophic velocity',arange(0,0.6,0.02))
plot12all('geostrophic velocity')
seasonline_all('geostrophic velocity',200,'[m/s]',' averaged over top 200m','mean','mgeovel200')

univec['velocity difference']=['velocity difference',arange(-0.3,0.305,0.02),cm.RdBu_r,arange(-0.3,0.305,0.1),'[m/s]']


meanplot('velocity difference')
stdplot('velocity difference',arange(0,0.2,0.01))
plot12all('velocity difference')
seasonline_all('velocity difference',200,'[m/s]',' averaged over top 200m','mean','veldiff200')


seasonline_all('absolute velocity',200,'[m/s]',' averaged over top 200m','mean','mabsvel200')


meanplot('turner angle')
stdplot('turner angle',arange(0,190,15))

plot12all('turner angle')



seasonline_all('turner angle',100,'[$^\circ$]',' at 100m','1depth','turner100',tu=1)

seasonline_all('density ratio',200,'',' mean over top 200m','mean','R200mean')

seasonline_all('density gradient',100,'',' at 100m','1depth','1d_dendiff100')

seasonline_all('density gradient',100,'',' mean over top 100m','mean','mdendiff100')

seasonline_all('density gradient',200,'',' mean over top 200m','mean','mdendiff200')

seasonline_all('beta DeltaS',100,'[$m^{-1}$]',' mean over top 100m','mean','mbetaS100')

seasonline_all('beta',200,'',' mean over top 200m','mean','beta200mean')

seasonline_all('alpha',200,'',' mean over top 200m','mean','alpha200mean')

seasonline_all('DeltaS',200,'',' mean over top 200m','mean','mDeltaS200')

seasonline_all('DeltaT',200,'',' mean over top 200m','mean','mDeltaTS200')

seasonline_all('beta DeltaS',200,'[$m^{-1}$]',' mean over top 200m','mean','mbetaS200')

seasonline_all('alpha DeltaT',100,'[$m^{-1}$]',' mean over top 100m','mean','malphaT100')

seasonline_all('alpha DeltaT',200,'[$m^{-1}$]',' mean over top 200m','mean','malphaT200')

seasonline_all('turner angle noden',100,'[$^\circ$]',' mean over top 100m','mean','turner100mean_noden',tu=1)

90+45

seasonline_all('turner angle noden',200,'[$^\circ$]',' mean over top 200m','mean','turner200mean_noden',tu=1)

seasonline_all('turner angle',200,'[$^\circ$]',' mean over top 200m','mean','turner200mean',tu=1)


# scatter(betaS[:,0,:],alphaT[:,0,:],c=turner[:,0,:],cmap=cm.RdBu_r,vmin=-90,vmax=90);
# colorbar(ticks=range(-90,100,45))
# axvline(0,color='k',alpha=0.5)
# axhline(0,color='k',alpha=0.5)
# plot(range(-1,2),range(-1,2),color='k',alpha=0.5)
# plot(range(-1,2),range(2,-1,-1),color='k',alpha=0.5)
# xlim([-0.002,0.002])
# ylim([-0.002,0.002])
