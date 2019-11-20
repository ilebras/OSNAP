from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'

dat={}
for ii in range(1,8):
    dat[ii]=xr.open_dataset(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(ii)+'_2018recovery_hourlymerged.nc')


[date,month,prs,salfinal,tmp]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/TSdailydic/TS_daily_dic_wcorr_15min.pickle','rb'))

for ii in range(1,8):
    figure(figsize=(12,4))
    for kk in date[ii]:
        plot(date[ii][kk],prs[ii][kk])
    gca().invert_yaxis()
    plot(dat[ii].TIME,dat[ii].PRES.T)
    title('CF'+str(ii))
    ylabel('pressure [db]')
    savefig(figdir+'pressure_overview/prs_synth_CF'+str(ii)+'.png')
