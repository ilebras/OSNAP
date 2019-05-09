from aux_funcs import *

CF=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_M_2014-2016_hourlyTSD_1903.pickle','rb'))

dat=io.loadmat(datadir+'OSNAP2016recovery/LS/LSgridded_TS.mat')

dat

eddy_stats=io.loadmat(datadir+'OSNAP2016recovery/Eddies/LS_eddy_stats.mat')

eddy_stats

eddy_date = array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in eddy_stats['time'][0]])

loc_cycl=eddy_stats['loc_cycl']

shape(loc_cycl)

dat.keys()

dat['date']=array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in dat['time'][0]])


ii=5

nnind=(~isnan(dat['D'][:,ii,3000]))

datemat=tile(dat['date'],[sum(nnind),1])

dlist=[datetime.datetime(2014,9,1),datetime.datetime(2015,2,1),datetime.datetime(2015,9,1),datetime.datetime(2016,3,1),datetime.datetime(2016,8,1)]

matplotlib.rcParams['ps.fonttype'] = 42

rcParams['mathtext.fontset'] = 'cm'

Dmat=dat['D'][nnind,ii,:]
Dbot=Dmat[-1,:]
Dmat[-1,isnan(Dbot)]=nanmean(Dbot)
sum(isnan(Dmat[-1,:]))


figure(figsize=(16,4))
contourf(datemat,-Dmat,dat['T'][nnind,ii,:],cmap=cm.RdYlBu_r,levels=arange(1.5,6.1,0.1),extend='both')
colorbar(ticks=range(0,7,1),label='Potential temperature [$^\circ$ C]')
contour(datemat,-Dmat,dat['R'][nnind,ii,:],colors='k',levels=arange(27.6,27.9,0.05))
plot(eddy_date[(loc_cycl[-2,:]>0)==True],0*ones(sum(loc_cycl[-2,:]>0)),'ko')
ylabel('depth [m]')
title('LS'+str(ii+1))
for dd in range(len(dlist[:-1])):
    xlim(dlist[dd],dlist[dd+1])
    savefig('../figures/Eddies/hov/Tdenhov_abs_LS'+str(ii+1)+'_d'+str(dd)+'.png',bbox_inches='tight',dpi=300)
    savefig('../figures/Eddies/hov/Tdenhov_abs_LS'+str(ii+1)+'_d'+str(dd)+'.pdf',bbox_inches='tight')
    # savefig('../figures/Eddies/hov/Tdenhov_abs_LS'+str(ii+1)+'_d'+str(dd)+'.ps',bbox_inches='tight')

# def np64ToDatetime(DA):
#     return [datetime.datetime.utcfromtimestamp((dd-np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for dd in DA]
#
# dcheck=array(np64ToDatetime(dat[6].date))
#
# datemat=tile(dcheck,[8,1])

# for ii in range(6,9):
#         #
#         # figure(figsize=(16,4))
#         # for jj in range(len(dat[ii].dpvec)-2):
#         #     scatter(dcheck,dat[ii].depth[jj,:].T,c=(dat[ii].ptmp[jj,:].values-mean(dat[ii].ptmp[jj,:]).values),s=2,vmin=-1,vmax=1,cmap=cm.RdYlBu_r);
#         # plot(eddy_date[(loc_cycl[ii-5,:]>0)==True],-2000*ones(sum(loc_cycl[ii-5,:]>0)),'k.')
#         # colorbar(label='Potential temperature anomaly [$^\circ$ C]')
#         # title('CF'+str(ii))
#         # savefig('../figures/Eddies/hov/Thov_anom_CF'+str(ii)+'.png',bbox_inches='tight',dpi=300)
#         #
#         # figure(figsize=(16,4))
#         # for jj in range(len(dat[ii].dpvec)-2):
#         #     scatter(dcheck,dat[ii].depth[jj,:].T,c=(dat[ii].sigma0[jj,:].values-mean(dat[ii].sigma0[jj,:]).values),s=2,vmin=-0.1,vmax=0.1,cmap=cm.BrBG)#,vmin=-5,vmax=5,cmap=cm.RdBu_r);
#         # plot(eddy_date[(loc_cycl[ii-5,:]>0)==True],-2000*ones(sum(loc_cycl[ii-5,:]>0)),'k.')
#         # colorbar(label='$\sigma_0$ anomaly [kg m$^{-3}$]')
#         # title('CF'+str(ii))
#         # savefig('../figures/Eddies/hov/Dhov_anom_CF'+str(ii)+'.png',bbox_inches='tight',dpi=300)
#         #
#         # figure(figsize=(16,4))
#         # for jj in range(len(dat[ii].dpvec)-2):
#         #     scatter(dcheck,dat[ii].depth[jj,:].T,c=(dat[ii].ptmp[jj,:].values),s=2,vmin=0,vmax=7,cmap=cm.RdYlBu_r);
#         # plot(eddy_date[(loc_cycl[ii-5,:]>0)==True],-2000*ones(sum(loc_cycl[ii-5,:]>0)),'k.')
#         # colorbar(label='Potential temperature [$^\circ$ C]')
#         # title('CF'+str(ii))
#         # for dd in range(len(dlist[:-1])):
#         #     xlim(dlist[dd],dlist[dd+1])
#         #     savefig('../figures/Eddies/hov/Thov_abs_CF'+str(ii)+'_d'+str(dd)+'.png',bbox_inches='tight',dpi=300)
#
#         # figure(figsize=(16,4))
#         # for jj in range(len(dat[ii].dpvec)-2):
#         #     scatter(dcheck,dat[ii].depth[jj,:].T,c=(dat[ii].sigma0[jj,:].values),s=2,cmap=cm.YlGnBu,vmin=27.3,vmax=28)#,vmin=-5,vmax=5,cmap=cm.RdBu_r);
#         # plot(eddy_date[(loc_cycl[ii-5,:]>0)==True],-2000*ones(sum(loc_cycl[ii-5,:]>0)),'k.')
#         # colorbar(label='$\sigma_0$ [kg m$^{-3}$]')
#         # title('CF'+str(ii))
#         # savefig('../figures/Eddies/hov/Dhov_abs_CF'+str(ii)+'.png',bbox_inches='tight',dpi=300)
#
#
#
# eddy_date[0]
#
#
#
#
#
# tvec=arange(0,len(dcheck)/24,1./24)
# ttest=dat[6].ptmp[-1,:]
# [t_fit,t_std,t_period]=fitsin(tvec,ttest,mean(ttest).values,60,std(ttest).values,365.25)
#
# ttest.plot()
# plot(dcheck,t_fit)
#
# plot(dcheck,ttest-t_fit)
# axhline(0)
#
# XXXXXXXXXXXXXXXXXXXXXXXXXX
#
#
# daily=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1809lpfilt_noextrap_wMLPV.pickle','rb'))
#
#
#
# for ii in range(5,8):
#     # figure(figsize=(14,3))
#     # pcolor(daily.date,daily.depth[250:],daily.temperature[ii,250:,:])
#     # plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
#     # ylabel('depth [m]')
#     # ylim(2e3,500)
#     # colorbar(label='potential temperature [$^\circ$C]')
#     # title('CF'+str(ii+1))
#     # savefig('../figures/Eddies/hov/Thov_abs_dailygridded_CF'+str(ii+1)+'.png',bbox_inches='tight',dpi=300)
#     #
#     # figure(figsize=(14,3))
#     # pcolor(daily.date,daily.depth[250:],daily.temperature[ii,250:,:]-daily.temperature[ii,250:,:].mean(dim='date'),vmin=-1,vmax=1,cmap=cm.RdYlBu_r)
#     # plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
#     # ylabel('depth [m]')
#     # ylim(2e3,500)
#     # colorbar(label='potential temperature anomaly [$^\circ$C]')
#     # title('CF'+str(ii+1))
#     # savefig('../figures/Eddies/hov/Thov_anom_dailygridded_CF'+str(ii+1)+'.png',bbox_inches='tight',dpi=300)
#
#     figure(figsize=(14,3))
#     pcolor(daily.date,daily.depth[250:],daily['potential density'][ii,250:,:],cmap=cm.YlGnBu)
#     plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
#     ylabel('depth [m]')
#     ylim(2e3,500)
#     colorbar(label='$\sigma_0$ [kg m$^{-3}$]')
#     title('CF'+str(ii+1))
#     savefig('../figures/Eddies/hov/Dhov_abs_dailygridded_CF'+str(ii+1)+'.png',bbox_inches='tight',dpi=300)
#
#     # figure(figsize=(14,3))
#     # pcolor(daily.date,daily.depth[250:],daily['potential density'][ii,250:,:]-daily['potential density'][ii,250:,:].mean(dim='date'),vmin=-0.1,vmax=0.1,cmap=cm.RdYlBu_r)
#     # plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
#     # ylabel('depth [m]')
#     # ylim(2e3,500)
#     # colorbar(label='$\sigma_0$ anomaly [kg m$^{-3}$]')
#     # title('CF'+str(ii+1))
#     # savefig('../figures/Eddies/hov/Dhov_anom_dailygridded_CF'+str(ii+1)+'.png',bbox_inches='tight',dpi=300)
#
#
# eddy_date
#
#
# XXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ii=5
# dd=1
#
#
# pcolor(daily.date,daily.depth[250:],daily.temperature[ii,250:,:],zorder=1)
#
# help(contour)
# contour(daily['potential density'][ii,250:,:],colors='k')
#
# figure(figsize=(14,3))
# pcolor(daily.date,daily.depth[250:],daily.temperature[ii,250:,:],zorder=1)
# colorbar(label='potential temperature [$^\circ$C]')
# contour(daily.date,daily.depth[250:],daily['potential density'][ii,250:,:],colors='k',zorder=2)
# plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
# ylabel('depth [m]')
# ylim(2e3,500)
# title('CF'+str(ii+1))
# for dd in range(len(dlist[:-1])):
#     xlim(dlist[dd],dlist[dd+1])
#     savefig('../figures/Eddies/hov/Thov_abs_dailygridded_CF'+str(ii)+'_d'+str(dd)+'.png',bbox_inches='tight',dpi=300)
#
# for ii in range(5,8):
#     figure(figsize=(14,3))
#     pcolor(daily.date,daily.depth[250:],daily.temperature[ii,250:,:],zorder=1)
#     colorbar(label='potential temperature [$^\circ$C]')
#     contour(daily.date,daily.depth[250:],daily['potential density'][ii,250:,:],colors='k',zorder=2)
#     plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
#     ylabel('depth [m]')
#     ylim(2e3,500)
#     title('CF'+str(ii+1))
#     for dd in range(len(dlist[:-1])):
#         xlim(dlist[dd],dlist[dd+1])
#         savefig('../figures/Eddies/hov/Thov_abs_dailygridded_CF'+str(ii)+'_d'+str(dd)+'.png',bbox_inches='tight',dpi=300)
#     # savefig('../figures/Eddies/hov/Thov_abs_dailygridded_CF'+str(ii+1)+'.png',bbox_inches='tight',dpi=300)
#     #
#     # figure(figsize=(14,3))
#     # pcolor(daily.date,daily.depth[250:],daily.temperature[ii,250:,:]-daily.temperature[ii,250:,:].mean(dim='date'),vmin=-1,vmax=1,cmap=cm.RdYlBu_r)
#     # plot(eddy_date[(loc_cycl[ii-4,:]>0)==True],600*ones(sum(loc_cycl[ii-4,:]>0)),'k.')
#     # ylabel('depth [m]')
#     # ylim(2e3,500)
#     # colorbar(label='potential temperature anomaly [$^\circ$C]')
#     # title('CF'+str(ii+1))
#     # savefig('../figures/Eddies/hov/Thov_anom_dailygridded_CF'+str(ii+1)+'.png',bbox_inches='tight',dpi=300)
