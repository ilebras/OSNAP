from aux_funcs import *

winddat=pickle.load(open('../pickles/wind/NARR_linet_newtheta.pickle','rb'))

## wind is from NARR: North American Regional Reanalysis Project
## done by US National Weather Service (NCEP)
## wind is at 10m in units of m/s


winddaily=winddat.resample('D',dim='date',how='mean')
windweekly=winddat.resample('W',dim='date',how='mean')
windmonthly=winddat.resample('M',dim='date',how='mean')

colors=pal.cubehelix.perceptual_rainbow_16.get_mpl_colormap()

plot(CFlon,CFlat,'o-',label='CF mooring positions')
plot(winddat.longitude,winddat.latitude,'o-',label='Wind speed positions')

legend()

figure(figsize=(12,3))
for ii in range(5):
    plot(winddaily['date'],winddaily['along track wind speed'][:,ii],alpha=0.2,label='',color=colors(ii*50))
    plot(windweekly['date'],windweekly['along track wind speed'][:,ii],label=str(float(winddat.distance[ii])),color=colors(ii*50))
    # plot(windmonthly['date'],windmonthly['across track wind speed'][:,ii],label=str(float(winddat.distance[ii])),color=colors(ii*50),linewidth=3)
legend(loc=(1.05,0.2))
axhline(-5)
title('along track wind speed')
savefig('../figures/wind/all_alongwind.png',bbox_inches='tight')


figure(figsize=(12,3))
for ii in range(5):
    plot(winddaily['date'],winddaily['across track wind speed'][:,ii],alpha=0.2,label='',color=colors(ii*50))
    plot(windweekly['date'],windweekly['across track wind speed'][:,ii],label=str(float(winddat.distance[ii])),color=colors(ii*50))
    # plot(windmonthly['date'],windmonthly['across track wind speed'][:,ii],label=str(float(winddat.distance[ii])),color=colors(ii*50),linewidth=3)
legend(loc=(1.05,0.2))
axhline(-5)
title('across track wind speed')
savefig('../figures/wind/all_acrosswind.png',bbox_inches='tight')




mydateticks=[datetime.datetime(2014,8,1)+relativedelta(months=3*ii) for ii in range(9)]

figure(figsize=(12,3))
ii=2
plot(winddaily['date'],winddaily['along track wind speed'][:,ii],alpha=0.5,label='',color='red')
ax=plot(windweekly['date'],windweekly['along track wind speed'][:,ii],color='red')
# gca().set_xticks(mydateticks)
grid('on')
axhline(-5,color='k')
title('central along track wind speed')
savefig('../figures/wind/center_alongwind.png',bbox_inches='tight')

figure(figsize=(12,3))
ii=2
plot(winddaily['date'],winddaily['across track wind speed'][:,ii],alpha=0.5,label='',color='purple')
ax=plot(windweekly['date'],windweekly['across track wind speed'][:,ii],color='purple')
# gca().set_xticks(mydateticks)
grid('on')
axhline(-5,color='k')
title('central across track wind speed')
savefig('../figures/wind/center_acrosswind.png',bbox_inches='tight')

figure(figsize=(12,3))
ii=2
# ax=plot(windmonthly['date'],windmonthly['across track wind speed'][:,ii],color='purple')
# ax=plot(windmonthly['date'],-windmonthly['along track wind speed'][:,ii],color='red')
ax=plot(winddat['date'],winddat['across track wind speed'][:,ii],color='purple',alpha=0.5,label='')
ax=plot(winddaily['date'],-winddaily['along track wind speed'][:,ii],color='blue',alpha=0.5,label='')
ax=plot(windweekly['date'],windweekly['across track wind speed'][:,ii],color='purple')
ax=plot(windweekly['date'],-windweekly['along track wind speed'][:,ii],color='darkblue')
# gca().set_xticks(mydateticks)
grid('on')
axhline(-5,color='k')
axhline(5,color='k')
ylim([-20,20])
ylabel('[m/s]')
legend()
savefig('../figures/wind/center_alongandacrosswind.png',bbox_inches='tight')
savefig('../figures/wind/center_alongandacrosswind.pdf',bbox_inches='tight')
