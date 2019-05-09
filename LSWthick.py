
from aux_funcs import *
### Take a look at LSW layer thickness: sigma_0= 27.74 - 27.8

dat=pickle.load(open('../pickles/CF_xarray_notid_pdenfix.pickle','rb'))


## Get time series at each mooring
lswmin=27.74
lswmax=27.8

denlim=pd.DataFrame(index=dat.date,columns=dat.distance)


def find_den(denbnd):
    denlim=pd.DataFrame(index=dat.date,columns=dat.distance[3:8])
    for ii in range(5):
        print(ii)
        for tt in enumerate(dat.date):
            wherevec=where(dat['potential density'][ii+3,:,tt[0]]>denbnd)[0]
            if (len(wherevec)>1):
                denlim.iloc[tt[0],ii]=float(dat.depth[wherevec[0]])
            elif (sum(~isnan(dat['potential density'][ii+3,:,tt[0]]))==0):
                denlim.iloc[tt[0],ii]=nan
            else:
                denlim.iloc[tt[0],ii]=float(max(dat.depth[~isnan(dat['potential density'][ii+3,:,tt[0]])]))

    denlim=denlim.apply(pd.to_numeric)

    return denlim


mind_lsw=find_den(lswmin)

mind_lsw.plot()

maxd_lsw=find_den(lswmax)

maxd_lsw.plot()

ax1=subplot(111)
mind_lsw.plot(ax=ax1)
maxd_lsw.plot(ax=ax1)
mind_lsw.resample('M',how='mean').plot(ax=ax1)
maxd_lsw.resample('M',how='mean').plot(ax=ax1)
savefig('../figures/LSWthick/LSWthickness.pdf')

figure(figsize=(10,3.5))
ax1=subplot(111)
(maxd_lsw-mind_lsw).iloc[:,2:].resample('W').mean().plot(ax=ax1)
title('Labrador Sea Water layer thickness')
ylabel('[m]')
grid('on')
savefig('../figures/LSWthick/LSWthickness_weekly.pdf')

pickle.dump([mind_lsw,maxd_lsw],open('../pickles/aux/lswbnds.pickle','wb'))
