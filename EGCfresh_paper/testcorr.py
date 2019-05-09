from aux_funcs import *

[cc,egic,eg,ic]=pickle.load(open('../pickles/transdic_1810JHcal.pickle','rb'))

era,corrmat=pickle.load(open('../pickles/wind_era5.pickle','rb'))
erasub=era.sel(date=slice('2014-08-17','2016-07-28'))

def corrcalc(f1,f2str):
    filt=5
    N  = 2    # Filter order
    Wn = 1./filt # Cutoff frequency (5 days)
    BS, AS = sig.butter(N, Wn, output='ba')

    corrmat=zeros((len(era.lon),len(era.lat)))
    for ii,lonn in enumerate(era.lon):
        for jj,latt in enumerate(era.lat):
            # corrmat[ii,jj]=corrcoef(sig.filtfilt(BS,AS,f1),sig.filtfilt(BS,AS,erasub[f2str].sel(lon=lonn).sel(lat=latt)))[0,1]
            corrmat[ii,jj]=corrcoef(f1,erasub[f2str].sel(lon=lonn).sel(lat=latt))[0,1]
    # corrmat[abs(corrmat)<0.232]=nan
    return corrmat

corrmat={}
corrmat['egictrans-curl']=corrcalc(egic['trans'],'curl')
corrmat['egictrans-stress']=corrcalc(egic['trans'],'tau along')
# corrmat['egtrans-curl']=corrcalc(eg['trans'],'curl')
# corrmat['egtrans-stress']=corrcalc(eg['trans'],'tau along')
corrmat['cctrans-stress']=corrcalc(cc['trans'],'tau along')
corrmat['cctrans-curl']=corrcalc(cc['trans'],'curl')
# corrmat['egictrans-south']=corrcalc(egic['trans'],'tauy')
# corrmat['cctrans-south']=corrcalc(cc['trans'],'tauy')

pickle.dump([era,corrmat],open('../pickles/wind_era5_dailycorrs.pickle','wb'))

plot(corrmat['cctrans-stress'])
