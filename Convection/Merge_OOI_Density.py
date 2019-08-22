from aux_funcs import *

ooi_all=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_all_props.pickle','rb'))

dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1808lpfilt.pickle','rb'))


#################################################################################################
######### For now, make a combined OOI density profile
#################################################################################################

# Combine top 90m of surface mooring with fla
# Check how well the merge works
# Compare with the other two moorings systematically

ooi_all['pden'].sel(moor='flb').plot()

plim=100
sfcden=ooi_all['pden'].sel(moor='sfc')[:,ooi_all.prs<=plim]
fladen=ooi_all['pden'].sel(moor='fla')[:,ooi_all.prs>plim]
sfcsal=ooi_all['sal'].sel(moor='sfc')[:,ooi_all.prs<=plim]
flasal=ooi_all['sal'].sel(moor='fla')[:,ooi_all.prs>plim]
sfctmp=ooi_all['ptmp'].sel(moor='sfc')[:,ooi_all.prs<=plim]
flatmp=ooi_all['ptmp'].sel(moor='fla')[:,ooi_all.prs>plim]

denmerge=hstack((sfcden.values,fladen.values))
tmpmerge=hstack((sfctmp.values,flatmp.values))
salmerge=hstack((sfcsal.values,flasal.values))



pcolor(ooi_all['pden'].sel(moor='fla'))
colorbar()



pind=int(plim/5-1)
pind2=int(40/5-1)
pind3=int(550/5-1)
ooiden=xr.Dataset({ 'pden': (['date','prs'],denmerge),'ptmp': (['date','prs'],tmpmerge),'sal': (['date','prs'],salmerge),'depth': (['prs'],-gsw.z_from_p(ooi_all.prs.values,60))},coords={'date': ooi_all.date.values, 'prs': ooi_all.prs.values})
for var in ['pden','ptmp','sal']:
    ooiden[var][:,pind:pind+5]=NaN
    ooiden[var][:,pind2-3:pind2+5]=NaN
    ooiden[var][:,pind3-30:pind3+30]=NaN
ooiden=ooiden.interpolate_na(dim='prs')

for ii in range(14):
    figure()
    plot(ooiden.pden[ii*100:(ii+1)*100,:].T,ooiden.prs);
    title(str(ooiden.date[ii*100].values)[:10]+' -- '+str(ooiden.date[(ii+1)*100].values)[:10])
    axhline(100,color='grey')
    [axvline(dd,color='grey') for dd in arange(d1,d3,0.01)]
    [axvline(dd,color='k') for dd in [d1,d2,d3]]
    ylim(1500,0)
    xlim(27.4,27.85)

for ii in range(14):
    figure()
    plot(ooiden.pden[ii*100:(ii+1)*100,:].T,ooiden.prs);
    title(str(ooiden.date[ii*100].values)[:10]+' -- '+str(ooiden.date[(ii+1)*100].values)[:10])
    # axhline(100,color='grey')
    [axvline(dd,color='grey') for dd in arange(d1,d3,0.01)]
    [axvline(dd,color='k') for dd in [d1,d2,d3]]
    ylim(1500,0)
    xlim(27.4,27.85)

    figure()
    plot(ooiden.ptmp[ii*100:(ii+1)*100,:].T,ooiden.prs);
    title(str(ooiden.date[ii*100].values)[:10]+' -- '+str(ooiden.date[(ii+1)*100].values)[:10])
    ylim(1500,0)

    figure()
    plot(ooiden.sal[ii*100:(ii+1)*100,:].T,ooiden.prs);
    title(str(ooiden.date[ii*100].values)[:10]+' -- '+str(ooiden.date[(ii+1)*100].values)[:10])
    ylim(1500,0)


def compdenstats():
    figure(figsize=(7,8))

    for ii,mm in enumerate(['prf','fla','flb']):
        ii+=1
        plot(ooi_all.pden.sel(moor=mm).mean(dim='date'),ooi_all.prs,color='C'+str(ii),linewidth=3,label=mm)
        plot(ooi_all.pden.sel(moor=mm).mean(dim='date')-ooi_all.pden.sel(moor=mm).std(dim='date'),ooiden.prs,color='C'+str(ii),label='')
        plot(ooi_all.pden.sel(moor=mm).mean(dim='date')+ooi_all.pden.sel(moor=mm).std(dim='date'),ooiden.prs,color='C'+str(ii),label='')

    plot(ooiden.pden.mean(dim='date'),ooiden.prs,'k',linewidth=3,label='merged')
    plot(ooiden.pden.mean(dim='date')-ooiden.pden.std(dim='date'),ooiden.prs,'k',label='')
    plot(ooiden.pden.mean(dim='date')+ooiden.pden.std(dim='date'),ooiden.prs,'k',label='')
    legend()
    # plot(dat['potential density'][-1,:,:].mean(dim='date'),dat.depth,label='M1',linewidth=3,color='red')
    axhline(90,color='grey')
    ylim(2000,50)
    xlim(27.4,27.85)

compdenstats()

oom=ooiden.resample('1D',dim='date')
pickle.dump(ooiden,open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_denmerged_xray.pickle','wb'),protocol=2)
