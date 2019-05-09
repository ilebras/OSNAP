from aux_funcs import *
from scipy.stats import t as ttest

# def load_daily_file(fieldfile):
#     files = sort(hstack((glob.glob('/home/isabela/DATA/ERA5/2014/*'+fieldfile+'*'),glob.glob('/home/isabela/DATA/ERA5/2015/*'+fieldfile+'*'),glob.glob('/home/isabela/DATA/ERA5/2016/*'+fieldfile+'*'))))
#     datev=array([])
#     for num,file in enumerate(files):
#         dat=io.loadmat(file)
#         #average each month of data to daily
#         latv=dat['lat'].flatten()
#         lonv=dat['lon'].flatten()
#         dayv=dat['day_ts'].flatten()
#         daysm=unique(dat['day_ts'])
#         mnth_tmp=zeros((len(daysm),len(latv),len(lonv)))
#         date_tmp=array([datetime.datetime(2000,1,1)]*len(daysm))
#         field=fieldfile+'_ts'
#         for ii in daysm:
#             mnth_tmp[ii-1,:,:]=mean(dat[field][dayv==ii,:,:],axis=0)
#             date_tmp[ii-1]=datetime.datetime(unique(dat['year_ts'])[0],unique(dat['month_ts'])[0],ii)
#         if num==0:
#             dailyfield=mnth_tmp
#         else:
#             dailyfield=vstack((dailyfield,mnth_tmp))
#         datev=hstack((datev,date_tmp))
#
#     return dailyfield,lonv,latv,datev
#
# taux,lon,lat,date=load_daily_file('tx')
# tauy,lon,lat,date=load_daily_file('ty')
# icec,lon,lat,date=load_daily_file('icec')
#
# shape(icec)
# shape(taux)
# shape(lon)
# shape(lat)
# shape(date)
#
# era=xr.Dataset({'taux': (['date','lat','lon'],  taux),
#                 'tauy': (['date','lat','lon'],  tauy),
#                 'icec': (['date','lat','lon'],  icec)},
#                 coords={'date':date.astype('M'),
#                         'lat': lat,
#                         'lon': lon})
#
# ### Get along-flow wind stress near the array and plot time series
# era['tau along']=-(era['taux']*cos(theta)+era['tauy']*sin(theta))
#
#
# ### Get wind stress curl
# dlon=vstack((gsw.distance(lonmat,latmat)[:,0],gsw.distance(lonmat,latmat).T)).T
# dlat=vstack((gsw.distance(lonmat.T,latmat.T)[:,0],gsw.distance(lonmat.T,latmat.T).T))
# curl=gradient(era.tauy)[2]/dlon+gradient(era.taux)[1]/dlat #note sign change because lat is decreasing...
# era['curl'] = era.taux.copy()
# era['curl'][:]=curl*1e6
#
# Get correlation with transports
# erasub=era.sel(date=slice('2014-08-17','2016-07-28'))
#
# [cc,egic,eg,ic]=pickle.load(open('../pickles/transdic_1810JHcal.pickle','rb'))
#
# def corrcalc(f1,f2str):
#     filt=5
#     N  = 2    # Filter order
#     Wn = 1./filt # Cutoff frequency (5 days)
#     BS, AS = sig.butter(N, Wn, output='ba')
#
#     corrmat=zeros((len(era.lon),len(era.lat)))
#     for ii,lonn in enumerate(era.lon):
#         for jj,latt in enumerate(era.lat):
#             corrmat[ii,jj]=corrcoef(sig.filtfilt(BS,AS,f1),sig.filtfilt(BS,AS,erasub[f2str].sel(lon=lonn).sel(lat=latt)))[0,1]
#     # corrmat[abs(corrmat)<0.232]=nan
#     return corrmat
#
# corrmat={}
# corrmat['egictrans-curl']=corrcalc(egic['trans'],'curl')
# corrmat['egictrans-stress']=corrcalc(egic['trans'],'tau along')
# # corrmat['egtrans-curl']=corrcalc(eg['trans'],'curl')
# # corrmat['egtrans-stress']=corrcalc(eg['trans'],'tau along')
# corrmat['cctrans-stress']=corrcalc(cc['trans'],'tau along')
# corrmat['cctrans-curl']=corrcalc(cc['trans'],'curl')
# # corrmat['egictrans-south']=corrcalc(egic['trans'],'tauy')
# # corrmat['cctrans-south']=corrcalc(cc['trans'],'tauy')
#
# pickle.dump([era,corrmat],open('../pickles/wind_era5.pickle','wb'))


era,corrmat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/wind_era5.pickle','rb'))

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


lat_start=55
lat_end  =70
lon_start=-50
lon_end  =-20
lon_0= - (abs(lon_end)+abs(lon_start))/2.0
lon_0= - (abs(lon_end)+abs(lon_start))/2.0
lat_0= - (abs(lat_end)+abs(lat_start))/2.0

map = Basemap(lat_0=lat_0,lon_0=lon_0,
                llcrnrlat=lat_start,urcrnrlat=lat_end+1,
                llcrnrlon=lon_start,urcrnrlon=lon_end,
                resolution='l',projection='stere')
                # suppress_ticks=False)

lonmat,latmat=meshgrid(era.lon,era.lat)
x, y = map(lonmat,latmat)

skiparr=4
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

# labvec=[0,0,0,1]

contmap('2014-8-1','2016-8-1')



lonlab=arange(-45,-10,10)
latlab=arange(55,73,5)
[xlab,ylab]=map(lonlab,latlab)
[xlp,ylp]=map(-48,68.5)
[xlq,ylq]=map(-48,66)
ylp2=map(-48,67.5)[1]

def curlmap(date1,date2,ice=1):
    lon1=-40
    lon2=-30
    lat1=60
    lat2=65
    curlplot=era['curl'].sel(date=slice(date1,date2)).mean(dim='date')
    curlim=3
    contit=map.contourf(x,y,curlplot,arange(-3,3,0.25),cmap=cm.RdBu_r,extend='both')
    if ice==1:
        iceplot=era['icec'].sel(date=slice(date1,date2)).mean(dim='date')
        map.contour(x,y,iceplot,[0.1,20],colors='limegreen',linewidths=3)
    # map.colorbar(ticks=arange(0,0.8,0.1),label='Wind Stress magnitude [N m$^{-2}$]')
    map.plot(CFlon,CFlat,'w-',latlon=True,linewidth=10)
    map.plot(CFlon,CFlat,'k-',latlon=True,linewidth=3)
    # map.plot([lon1,lon1,lon2,lon2,lon1],[lat1,lat2,lat2,lat1,lat1],color='k')
    skiparr=5
    qq=map.quiver(x[::skiparr,::skiparr],y[::skiparr,::skiparr],era['taux'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],
                era['tauy'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],color='k',scale=2)
    map.fillcontinents(color='lightgrey')
    ylim(ylab[0],ylab[-1])
    return contit,map,qq

def plot3curl(iceq=1):
    tfs=20
    fig, (ax11, ax22, ax33) = plt.subplots(1,3, sharex=True, sharey=True,figsize=(10,4))
    subplot(131)
    c2,map,qq=curlmap('2014-10-1','2015-1-1',ice=iceq)
    gca().set_yticks(ylab)
    gca().set_yticklabels(['55$^\circ$N','60$^\circ$N','65$^\circ$N','70$^\circ$N'])
    gca().set_xticks(xlab[:-1])
    gca().set_xticklabels(['45$^\circ$W','35$^\circ$W','25$^\circ$W'])
    text(xlp,ylp,'Fall',fontsize=tfs)

    cbaxes = fig.add_axes([1.01, 0.1, 0.02, 0.65])
    colorbar(c2,label='Wind stress curl [ x $10^{-6}$ N m$^{-3}$]',cax=cbaxes,ticks=arange(-3,3,1))
    subplot(132)
    c3,map,qq=curlmap('2015-1-01','2015-4-1',ice=iceq)
    gca().set_xticks(xlab[:-1])
    gca().set_xticklabels(['45$^\circ$W','35$^\circ$W','25$^\circ$W'])
    text(xlp,ylp,'Winter',fontsize=tfs)
    subplot(133)
    c1,map,qq=curlmap('2015-04-01','2015-10-01',ice=iceq)
    plt.quiverkey(qq,1.2,0.85,0.2,'Wind stress \n 0.2 N m$^{-2}$',
                    labelpos='N',coordinates='axes',angle=10)
    text(xlp,ylp2,'Spring/ \nSummer',fontsize=tfs)
    gca().set_xticks(xlab[:-1])
    gca().set_xticklabels(['45$^\circ$W','35$^\circ$W','25$^\circ$W'])
    plt.tight_layout()

plot3curl()
savefig('../figures/paperfigs/ERAcurl_seas_daily.pdf',bbox_inches='tight')
#
# era,corrmat=pickle.load(open('../pickles/wind_era5_dailycorrs.pickle','rb'))

plot3curl(0)
savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/ERAcurl_seas_daily_noice.pdf',bbox_inches='tight')

## Do FDR statistical test
Neff=70
def FDR(fstr,xlimi,ylimi):
    rvals=abs(corrmat[fstr].flatten())
    tvals=rvals/(sqrt((1-rvals**2)/Neff))
    pvals=sort(ttest.sf(abs(tvals),Neff)*2)
    plot(pvals)
    fdrline=0.1*arange(len(pvals))/len(pvals)
    plot(fdrline)
    xlim(0,xlimi)
    ylim(0,ylimi)
    pfr=pvals[(pvals>fdrline)&(pvals>1e-4)][0]#stands for "p for real"
    rfr=sort(rvals)[::-1][(pvals>fdrline)&(pvals>1e-4)][0]
    title('Correct p value is: '+str(pfr)+'\n Minimum r vale is: '+str(rfr))
    return rfr

rfr={}
rfr['egictrans-curl']=FDR('egictrans-curl',2e4,0.03)


rfr['egictrans-stress']=FDR('egictrans-stress',2.5e4,0.05)
rfr['cctrans-curl']=FDR('cctrans-curl',5e4,0.5)
rfr['cctrans-stress']=FDR('cctrans-stress',2e3,0.008)

rfr

def corrmaplot(fstr,xx,yy,tit):
    # corrmat_nonan=corrmat[fstr]
    # corrmat_nonan[abs(corrmat_nonan)<0.232]=NaN
    colref=map.contourf(x,y,corrmat[fstr].T,arange(-0.8,0.81,0.05),cmap=cm.RdBu_r)
    if 'egic' in fstr:
        cc=map.contour(x,y,corrmat[fstr].T,levels=[-0.7,-0.6,-0.5,0.5,0.6,0.7],colors='k')
        map.contour(x,y,corrmat[fstr].T,levels=[-0.4,0.4],colors='k')
    elif 'cc' in fstr:
        cc=map.contour(x,y,corrmat[fstr].T,levels=[-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7],colors='k')
    map.contourf(x, y, corrmat[fstr].T, [-rfr[fstr],rfr[fstr]], colors='none',hatches=['.'],)
    clabel(cc,fmt='%1.1f')
    map.fillcontinents(color='lightgrey',zorder=2)
    map.plot(CFlon,CFlat,'w-',latlon=True,linewidth=10)
    map.plot(CFlon,CFlat,'k-',latlon=True,linewidth=3)

    xlp=-49.5
    tfs=18
    [xtit,ytit]=map(-48,68.8)
    text(xtit,ytit,tit,fontsize=tfs,zorder=3)

    if xx==1:
        gca().set_xticks(xlab[:-1])
        gca().set_xticklabels(['45$^\circ$W','35$^\circ$W','25$^\circ$W'])
    else:
        gca().set_xticks([])
    if yy==1:
        gca().set_yticks(ylab[1:])
        gca().set_yticklabels(['60$^\circ$N','65$^\circ$N','70$^\circ$N'])
    else:
        gca().set_yticks([])
    ylim(ylab[0],ylab[-1])
    return colref

def corrfig():
    fs=14
    fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(7,7))
    subplot(221)
    colref=corrmaplot('cctrans-stress',0,1,'a)')
    title('Coastal current',fontsize=fs+1)
    ylabel('Along-stream\nwind stress',fontsize=fs)
    subplot(222)
    colref=corrmaplot('egictrans-stress',0,0,'b)')
    title('Slope current',fontsize=fs+1)
    subplot(223)
    colref=corrmaplot('cctrans-curl',1,1,'c)')
    ylabel('Wind stress curl \n',fontsize=fs)
    subplot(224)
    colref=corrmaplot('egictrans-curl',1,0,'d)')
    cbaxes = fig.add_axes([1.01, 0.3, 0.04, 0.45])
    colorbar(colref,label='Correlation coefficent',cax=cbaxes,ticks=arange(-0.8,0.9,0.2))
    plt.tight_layout()
    savefig('../figures/paperfigs/ERA_corrmap.pdf',bbox_inches='tight')


corrfig()

def singlecorr(var):
    # figure(figsize=(4,3))
    cref=corrmaplot(var,1,1,'')
    colorbar(cref,label='Correlation coefficent')
    savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/ERAcorr_'+var+'.pdf',bbox_inches='tight')

singlecorr('cctrans-stress')

singlecorr('egictrans-curl')



savefig('/home/isabela/Documents/applications/2018_WHOI/interview/seminar/figures/ERAcorr_slopecurl.pdf',bbox_inches='tight')

# def corrfig2():
#     fs=14
#     fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(7,7))
#     subplot(221)
#     colref=corrmaplot('cctrans-stress',0,1,'a)')
#     title('Coastal current',fontsize=fs+1)
#     ylabel('Along-stream\nwind stress',fontsize=fs)
#     subplot(222)
#     colref=corrmaplot('egtrans-stress',0,0,'b)')
#     title('Irminger Current',fontsize=fs+1)
#     subplot(223)
#     colref=corrmaplot('cctrans-curl',1,1,'c)')
#     ylabel('Wind stress curl \n',fontsize=fs)
#     subplot(224)
#     colref=corrmaplot('egtrans-curl',1,0,'d)')
#     cbaxes = fig.add_axes([1.01, 0.3, 0.04, 0.45])
#     colorbar(colref,label='Correlation coefficent',cax=cbaxes,ticks=arange(-0.8,0.9,0.2))
#     plt.tight_layout()
#     savefig('../figures/paperfigs/ERA_corrmap_IC.pdf',bbox_inches='tight')
#
# corrfig2()
#
# CFcurl=era['curl'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon')
# UPcurl=era['curl'].sel(lat=slice(65,60)).sel(lon=slice(-40,-30)).mean(dim='lon').mean(dim='lat')
#
# def pltcurlt():
#         figure(figsize=(12,4))
#         CFcurl.plot(color='k',alpha=0.5,label='')
#         plot(era.date,sig.filtfilt(B,A, CFcurl),label='At CF moorings',color='k')
#         xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,8,15)])
#         ylabel('Wind stress curl [ x $10^{-6}$ N m$^{-3}$]')
#         title('')
#         colorstripes()
#         ylim([-2,8])
#         savefig('../figures/paperfigs/ERAcurl_tseries_daily.pdf',bbox_inches='tight')
#
# pltcurlt()
#
# CFrol=pd.rolling_std(CFcurl.values,3)
# CFvar=sig.filtfilt(B,A,CFrol[~isnan(CFrol)])
#
# UProl=pd.rolling_std(UPcurl.values,3)
# UPvar=sig.filtfilt(B,A,UProl[~isnan(CFrol)])
#
# def pltcurlt():
#         figure(figsize=(10,4))
#         # CFcurl.plot(color='k',alpha=0.5,label='')
#         plot(era.date,sig.filtfilt(B,A, CFcurl),label='Curl magnitude',color='k')
#         plot(era.date,sig.filtfilt(B,A, UPcurl),label='Curl magnitude upstream',color='grey')
#         plot(era.date[~isnan(CFrol)],CFvar,label='Rolling variance',color='r')
#         plot(era.date[~isnan(UProl)],UPvar,label='Rolling variance upstream',color='orange')
#         xlim([datetime.datetime(2014,8,15),datetime.datetime(2016,8,15)])
#         ylabel('Wind stress curl [ x $10^{-6}$ N m$^{-3}$]')
#         title('Curl seasonality')
#         colorstripes()
#         legend()
#         savefig('../figures/paperfigs/ERAcurl_tseries_daily.pdf',bbox_inches='tight')
#
# pltcurlt()
#
#
# figure(figsize=(12,4))
# era['sshf'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(color='k',label='Sensible heat flux')
# era['slhf'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon').plot(color='grey',label='Latent heat flux')
# ylabel('[$W/m^2$]')
# colorstripes()
#
#
# cc['trans'].date[-1]
#
#
# CFcurlsub=erasub['curl'].sel(lat=62.25).sel(lon=slice(-35,-25)).mean(dim='lon')
# CFstressub=erasub['tau along'].sel(lat=60).sel(lon=slice(-43,-40)).mean(dim='lon')
# shape(egic['trans'])
#
# cohere(CFcurlsub.values,egic['trans'].values)
#
# field=egic['trans']
#
# def cohplot(field1,field2,labit):
#     f, Cxy=signal.coherence(field1.values,field2.values,nperseg=len(erasub.date)/4)
#     plot(f,Cxy,label=labit)
#
# CFcurlsub
#
# figure(figsize=(12,4))
# cohplot(CFcurlsub,egic['trans'],'slope')
# cohplot(CFcurlsub,cc['trans'],'coastal')
# legend()
#
# 1/0.02
#
# figure(figsize=(12,4))
# cohplot(CFcurlsub,egic['freshb'],'slope')
# cohplot(CFcurlsub,cc['freshb'],'coastal')
# legend()
#
#
# figure(figsize=(12,4))
# cohplot(CFstressub,egic['trans'],'slope')
# cohplot(CFstressub,cc['trans'],'coastal')
# legend()
#
# plot(CFstressub)
