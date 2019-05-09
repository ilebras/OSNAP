from aux_funcs import *

dat=xr.open_dataset('../data/aux_data/ERA_1804/tau_hflux_180411.nc')

dat.history

era=dat.rename({'inss':'tauy','iews':'taux','longitude':'lon','latitude':'lat','time':'date'})
era['lon']=era['lon']-360

from mpl_toolkits.basemap import Basemap

lat_start=55
lat_end  =70
lon_start=-50

lon_end  =-20
lon_0= - (abs(lon_end)+abs(lon_start))/2.0

map = Basemap(llcrnrlat=lat_start,urcrnrlat=lat_end,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='l',projection='cyl',
                lat_1=lat_start,
                suppress_ticks=False)

lonmat,latmat=meshgrid(era.lon,era.lat)
x, y = map(lonmat,latmat)

# map.contourf(x,y,sqrt(era['taux']**2+era['tauy']**2).mean(dim='date'))

skiparr=2
def contmap(date1,date2):
    contit=map.contourf(x,y,sqrt(era['taux']**2+era['tauy']**2).sel(date=slice(date1,date2)).mean(dim='date'),
                        51,cmap=cm.inferno,vmin=0,vmax=0.75,extend='both')
    # map.colorbar(ticks=arange(0,0.8,0.1),label='Wind Stress magnitude [N m$^{-2}$]')
    map.plot(CFlon,CFlat,'w-',latlon=True,linewidth=8)
    map.plot(CFlon,CFlat,'k-',latlon=True,linewidth=2)
    map.quiver(x[::skiparr,::skiparr],y[::skiparr,::skiparr],era['taux'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],
               era['tauy'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],color='white',scale=2)
    map.fillcontinents(color='lightgrey')
    return contit,map
    # savefig('../../confschools/1802_oceansciences/presentation/figures/'+date1+'-'+date2+'_windstress.pdf',bbox_inches='tight')

labvec=[0,0,0,1]


contmap('2014-10-1','2015-1-1')

def plot3stress():
    xlp=-49.5
    tfs=18
    fig, (ax11, ax22, ax33) = plt.subplots(3,1, sharex=True, sharey=True,figsize=(4,9))
    subplot(311)
    c2,map=contmap('2014-10-1','2015-1-1')
    map.drawmeridians(range(-45,-10,10),linewidth=0.01)
    map.drawparallels(range(55,80,10),labels=[1,0,0,0],linewidth=0.01)
    text(xlp,67,'Fall',fontsize=tfs)
    cbaxes = fig.add_axes([1.0, 0.3, 0.06, 0.45])
    colorbar(c2,ticks=arange(0,1,0.1),label='Wind Stress magnitude [N m$^{-2}$]',cax=cbaxes)
    subplot(312)
    c3,map=contmap('2015-1-01','2015-4-1')
    map.drawmeridians(range(-45,-10,10),linewidth=0.01)
    map.drawparallels(range(55,80,10),labels=[1,0,0,0],linewidth=0.01)
    text(xlp,67,'Winter',fontsize=tfs)
    subplot(313)
    c1,map=contmap('2015-04-01','2015-09-01')
    map.drawmeridians(range(-45,-10,10),labels=labvec,linewidth=0.01)
    map.drawparallels(range(55,80,10),labels=[1,0,0,0],linewidth=0.01)
    text(xlp,64.5,'Summer/ \nSpring',fontsize=tfs)
    # plt.tight_layout()
    savefig('../figures/paperfigs/ERAstress_seas_daily.pdf',bbox_inches='tight')

plot3stress()

def wint2():
    title('Not as big a tip jet the second winter!')
    contmap('2015-12-01','2016-3-1')

wint2()

### Get along-flow wind stress near the array and plot time series

era['tau along']=era['taux']*cos(theta)+era['tauy']*sin(theta)

def pltseries():
    figure(figsize=(10,4))
    era['tau along'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(label='At CF moorings',color='grey')
    era['tau along'].sel(lat=slice(65,60)).sel(lon=slice(-40,-35)).mean(dim='lon').mean(dim='lat').plot(label='Upstream',color='k')
    xlim([datetime.datetime(2014,8,15),datetime.datetime(2016,8,15)])
    legend()
    axhline(0,color='k')
    colorstripes()
    title('Tip jet feature dominates, none second winter')

pltseries()


### Get wind stress curl
dlon=vstack((gsw.distance(lonmat,latmat)[:,0],gsw.distance(lonmat,latmat).T)).T
dlat=vstack((gsw.distance(lonmat.T,latmat.T)[:,0],gsw.distance(lonmat.T,latmat.T).T))
curl=gradient(era.tauy)[2]/dlon+gradient(era.taux)[1]/dlat #note sign change because lat is decreasing...
era['curl'] = era.taux.copy()
era['curl'][:]=curl*1e6


def curlmap(date1,date2):
    lon1=-40
    lon2=-30
    lat1=60
    lat2=65
    contit=map.contourf(x,y,era['curl'].sel(date=slice(date1,date2)).mean(dim='date'),51,cmap=cm.RdBu_r,vmin=-3,vmax=3,extend='both')
    # map.colorbar(ticks=arange(0,0.8,0.1),label='Wind Stress magnitude [N m$^{-2}$]')
    map.plot(CFlon,CFlat,'w-',latlon=True,linewidth=8)
    map.plot(CFlon,CFlat,'k-',latlon=True,linewidth=2)
    # map.plot([lon1,lon1,lon2,lon2,lon1],[lat1,lat2,lat2,lat1,lat1],color='k')
    map.quiver(x[::skiparr,::skiparr],y[::skiparr,::skiparr],era['taux'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],
                era['tauy'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],color='k',scale=2)
    map.fillcontinents(color='lightgrey')
    return contit,map

def plot3curl():
    xlp=-49.5
    tfs=18
    fig, (ax11, ax22, ax33) = plt.subplots(3,1, sharex=True, sharey=True,figsize=(5,8))
    subplot(311)
    c2,map=curlmap('2014-10-1','2015-1-1')
    gca().set_xticks([])
    gca().set_yticks(range(55,70,5))
    gca().set_yticklabels(['55$^\circ$N','60$^\circ$N','65$^\circ$N','70$^\circ$N'])
    text(xlp,67,'Fall',fontsize=tfs)
    cbaxes = fig.add_axes([1.01, 0.3, 0.06, 0.45])
    colorbar(c2,label='Wind stress curl [ x $10^{-6}$ N m$^{-3}$]',cax=cbaxes)
    subplot(312)
    c3,map=curlmap('2015-1-01','2015-4-1')
    gca().set_xticks([])
    gca().set_yticks(range(55,70,5))
    gca().set_yticklabels(['55$^\circ$N','60$^\circ$N','65$^\circ$N','70$^\circ$N'])
    text(xlp,67,'Winter',fontsize=tfs)
    subplot(313)
    c1,map=curlmap('2015-04-01','2015-09-01')
    gca().set_xticks(range(-45,-20,10))
    gca().set_xticklabels(['45$^\circ$W','35$^\circ$W','25$^\circ$W'])
    gca().set_yticks(range(55,70,5))
    gca().set_yticklabels(range(55,70,5))
    gca().set_yticklabels(['55$^\circ$N','60$^\circ$N','65$^\circ$N','70$^\circ$N'])
    text(xlp,64.5,'Summer/ \nSpring',fontsize=tfs)
    plt.tight_layout()
    savefig('../figures/paperfigs/ERAcurl_seas_daily.pdf',bbox_inches='tight')


plot3curl()




erasub=era.sel(date=slice('2014-08-17','2016-07-28'))

[cc,egic,eg,ic]=pickle.load(open('../pickles/transdic.pickle','rb'))

def corrcalc(f1,f2str):
    filt=5
    N  = 2    # Filter order
    Wn = 1./filt # Cutoff frequency (5 days)
    BS, AS = sig.butter(N, Wn, output='ba')

    corrmat=zeros((len(era.lon),len(era.lat)))
    for ii,lonn in enumerate(era.lon):
        for jj,latt in enumerate(era.lat):
            corrmat[ii,jj]=corrcoef(sig.filtfilt(BS,AS,f1),sig.filtfilt(BS,AS,erasub[f2str].sel(lon=lonn).sel(lat=latt)))[0,1]
    # corrmat[abs(corrmat)<0.232]=nan
    return corrmat

corrmat={}
corrmat['egictrans-curl']=corrcalc(egic['trans'],'curl')
corrmat['egictrans-stress']=corrcalc(egic['trans'],'tau along')
corrmat['cctrans-stress']=corrcalc(cc['trans'],'tau along')
corrmat['cctrans-curl']=corrcalc(cc['trans'],'curl')

def corrmaplot(fstr,xx,yy,tit):
    colref=map.contourf(x,y,corrmat[fstr].T,31,vmin=-0.8,vmax=0.8,cmap=cm.RdBu_r)
    cc=map.contour(x,y,corrmat[fstr].T,levels=[-0.6,-0.4,0.4,0.6],colors='k')
    map.contourf(x, y, corrmat[fstr].T, [-0.232,0.232], colors='none',hatches=['.'],)
    clabel(cc,fmt='%1.1f')
    map.fillcontinents(color='lightgrey',zorder=2)
    xlp=-49.5
    tfs=18
    text(xlp,67.5,tit,fontsize=tfs,zorder=3)
    if xx==1:
        gca().set_xticks(range(-45,-20,10))
        gca().set_xticklabels(['45$^\circ$W','35$^\circ$W','25$^\circ$W'])
    else:
        gca().set_xticks([])
    if yy==1:
        gca().set_yticks(range(55,70,5))
        gca().set_yticklabels(['55$^\circ$N','60$^\circ$N','65$^\circ$N','70$^\circ$N'])
    else:
        gca().set_yticks([])
    return colref

def corrfig():
    fs=14
    fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(8.5,5))
    subplot(221)
    colref=corrmaplot('cctrans-stress',0,1,'a)')
    title('Coastal current',fontsize=fs+1)
    ylabel('Along-stream\nwind stress',fontsize=fs)
    subplot(222)
    colref=corrmaplot('egictrans-stress',0,0,'b)')
    title('Slope current',fontsize=fs+1)
    cbaxes = fig.add_axes([1.01, 0.3, 0.03, 0.45])
    colorbar(colref,label='Correlation coefficent',cax=cbaxes,ticks=arange(-0.8,0.9,0.2))
    subplot(223)
    colref=corrmaplot('cctrans-curl',1,1,'c)')
    ylabel('Wind stress curl \n',fontsize=fs)
    subplot(224)
    colref=corrmaplot('egictrans-curl',1,0,'d)')
    plt.tight_layout()
    savefig('../figures/paperfigs/ERA_corrmap.pdf',bbox_inches='tight')

corrfig()

CFcurl=era['curl'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon')
UPcurl=era['curl'].sel(lat=slice(65,60)).sel(lon=slice(-40,-30)).mean(dim='lon').mean(dim='lat')

def pltcurlt():
        figure(figsize=(12,4))
        CFcurl.plot(color='k',alpha=0.5,label='')
        plot(era.date,sig.filtfilt(B,A, CFcurl),label='At CF moorings',color='k')
        xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,8,15)])
        ylabel('Wind stress curl [ x $10^{-6}$ N m$^{-3}$]')
        title('')
        colorstripes()
        ylim([-2,8])
        savefig('../figures/paperfigs/ERAcurl_tseries_daily.pdf',bbox_inches='tight')

pltcurlt()

CFrol=pd.rolling_std(CFcurl.values,3)
CFvar=sig.filtfilt(B,A,CFrol[~isnan(CFrol)])

UProl=pd.rolling_std(UPcurl.values,3)
UPvar=sig.filtfilt(B,A,UProl[~isnan(CFrol)])

def pltcurlt():
        figure(figsize=(10,4))
        # CFcurl.plot(color='k',alpha=0.5,label='')
        plot(era.date,sig.filtfilt(B,A, CFcurl),label='Curl magnitude',color='k')
        plot(era.date,sig.filtfilt(B,A, UPcurl),label='Curl magnitude upstream',color='grey')
        plot(era.date[~isnan(CFrol)],CFvar,label='Rolling variance',color='r')
        plot(era.date[~isnan(UProl)],UPvar,label='Rolling variance upstream',color='orange')
        xlim([datetime.datetime(2014,8,15),datetime.datetime(2016,8,15)])
        ylabel('Wind stress curl [ x $10^{-6}$ N m$^{-3}$]')
        title('Curl seasonality')
        colorstripes()
        legend()
        savefig('../figures/paperfigs/ERAcurl_tseries_daily.pdf',bbox_inches='tight')

pltcurlt()


figure(figsize=(12,4))
era['sshf'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(color='k',label='Sensible heat flux')
era['slhf'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(color='grey',label='Latent heat flux')
ylabel('[$W/m^2$]')
colorstripes()


cc['trans'].date[-1]


CFcurlsub=erasub['curl'].sel(lat=62.25).sel(lon=slice(-35,-25)).mean(dim='lon')
CFstressub=erasub['tau along'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon')
shape(egic['trans'])

cohere(CFcurlsub.values,egic['trans'].values)

field=egic['trans']

def cohplot(field1,field2,labit):
    f, Cxy=signal.coherence(field1.values,field2.values,nperseg=len(erasub.date)/4)
    plot(f,Cxy,label=labit)

CFcurlsub

figure(figsize=(12,4))
cohplot(CFcurlsub,egic['trans'],'slope')
cohplot(CFcurlsub,cc['trans'],'coastal')
legend()

1/0.02

figure(figsize=(12,4))
cohplot(CFcurlsub,egic['freshb'],'slope')
cohplot(CFcurlsub,cc['freshb'],'coastal')
legend()


figure(figsize=(12,4))
cohplot(CFstressub,egic['trans'],'slope')
cohplot(CFstressub,cc['trans'],'coastal')
legend()

plot(CFstressub)
