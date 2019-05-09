#############################################################################
#############################################################################
############## Identify features of vertical structure ######################
######################### Modal decomp #####################################
#############################################################################

from aux_funcs import *

###########################################
################### Load data ############
##########################################

### Note: Took out flexibility to M1 mooring

moornum=7

moorname='CF'+str(moornum)

aqdlist=sort(glob.glob(datadir+'AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc'))

aqd={}
aqd['u']={}
aqd['v']={}
aqd['p']={}
aqd['date']={}
mindate=[]
maxdate=[]
datelen=[]
for ii,dd in enumerate(aqdlist):
        dat = Dataset(dd, 'r')
        thekey=int(mean(list(dat.variables['PRES'][:].flatten())))
        u=[float(tt) for tt in dat.variables['UCUR'][:].flatten()]
        v=[float(tt) for tt in dat.variables['VCUR'][:].flatten()]
        p=[float(tt) for tt in dat.variables['PRES'][:].flatten()]
        time_tmp=list(dat.variables['TIME'][:])
        #interpolate over nans
        nanind=where(isnan(u)==False)[0]
        print(mean(p),min(nanind),max(nanind),len(time_tmp))

        ufunc=interpolate.interp1d([time_tmp[nn] for nn in nanind],[u[nn] for nn in nanind])
        vfunc=interpolate.interp1d([time_tmp[nn] for nn in nanind],[v[nn] for nn in nanind])
        pfunc=interpolate.interp1d([time_tmp[nn] for nn in nanind],[p[nn] for nn in nanind])
        # mprs_tmp=nanmean(prs_tmp)*ones(len(time_tmp))
        if min(nanind)==1:
            time_tmp=time_tmp[1:]
        aqd['date'][thekey]=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=tt)
                        for tt in time_tmp])
        aqd['u'][thekey]=ufunc(time_tmp)
        aqd['v'][thekey]=vfunc(time_tmp)
        aqd['p'][thekey]=pfunc(time_tmp)
        mindate=hstack((min(aqd['date'][thekey]),mindate))
        maxdate=hstack((max(aqd['date'][thekey]),maxdate))
        datelen=hstack((len(aqd['date'][thekey]),datelen))

varvec=['u','v','p','date']

nomprsvec=[int(mean(aqd['p'][key])) for key in aqd['p']]

for nn in nomprsvec:
    dateind=(aqd['date'][nn]>=max(mindate))&(aqd['date'][nn]<=min(maxdate))
    for vv in varvec:
        aqd[vv][nn]=aqd[vv][nn][dateind]
        # print(vv,nn,len(aqd[vv][nn]))



adat=xr.Dataset({'u': (['nomprs', 'date'],  array([aqd['u'][key] for key in aqd['p']])),
                'v': (['nomprs', 'date'],  array([aqd['v'][key] for key in aqd['p']])),
                'p': (['nomprs', 'date'],  array([aqd['p'][key] for key in aqd['p']])),},
                coords={'nomprs': nomprsvec,
                        'date': aqd['date'][nomprsvec[0]]})

from scipy import signal
N  = 2    # Filter order
Wn = 0.01 # Cutoff frequency
B, A = signal.butter(N, Wn, output='ba')


#############################################################################
############## Plot time series and variance preserving spectra #############
##########################################################################
nomprsvec=sort(nomprsvec)[::-1]

pb=adat['p'].sel(nomprs=max(nomprsvec))

figure(figsize=(15,9))
subplot(311)
plot(adat.date,signal.filtfilt(B,A, pb),'k')
gca().set_xticklabels('')
gca().invert_yaxis()
ylabel('bottom pressure [db]')
title(moorname+' filtered (50 hr cutoff) time series',fontsize='x-large')
subplot(312)
for ii,nn in enumerate(nomprsvec):
    plot(adat.date,signal.filtfilt(B,A, adat['v'].sel(nomprs=nn)),label=str(nn)+' db',color=cm.YlGnBu(ii*30+80))
legend(loc=(1.02,0.2))
ylabel('Meridional velocity [m/s]')
gca().set_xticklabels('')
subplot(313)
for ii,nn in enumerate(nomprsvec):
    plot(adat.date,signal.filtfilt(B,A, adat['u'].sel(nomprs=nn)),label=str(nn)+' db',color=cm.YlGnBu(ii*30+80))
ylabel('Zonal velocity [m/s]')
savefig('../figures/VertModes/from_raw_data/FiltTseries_'+moorname+'.png',bbox_inches='tight')




from Spectrum_funcs import *

adat

freqs={}
ps={}
psd={}
for vv in varvec:
    freqs[vv]={}
    ps[vv]={}
    psd[vv]={}

def calcspec(var,nn):
    freqs[var][nn], ps[var][nn], psd[var][nn] = spectrum4(adat[var].sel(nomprs=nn),nsmooth=5)
    return freqs,ps,psd

varvec=['u','v','p']
for var in varvec:
    for nn in nomprsvec:
            freqs,ps,psd=calcspec(var,nn)

figure(figsize=(10,10))
frqall=freqs['p'][max(nomprsvec)]*48
frqind=[frqall<0.8]
subplot(311)
semilogx(frqall[frqind],frqall[frqind]*psd['p'][max(nomprsvec)][frqind],'k')
gca().set_xticklabels('')
title(moorname+' variance preserving spectra',fontsize='x-large')
subplot(312)
for ii,nn in enumerate(nomprsvec):
    semilogx(frqall[frqind],frqall[frqind]*psd['v'][nn][frqind],label=str(nn)+' db',color=cm.YlGnBu(ii*30+80))
legend(loc=(1.02,0.2))
ylabel('Variance (prop to)')
gca().set_xticklabels('')
subplot(313)
for ii,nn in enumerate(nomprsvec):
    semilogx(frqall[frqind],frqall[frqind]*psd['u'][nn][frqind],label=str(nn)+' db',color=cm.YlGnBu(ii*30+80))
xlabel('CPD')
savefig('../figures/VertModes/from_raw_data/VarSpec_'+moorname+'.png',bbox_inches='tight')


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True,figsize=(10,7))
for ii,nn in enumerate(nomprsvec):
        f, Cxy = signal.coherence(pb, adat['v'].sel(nomprs=nn),nperseg=4096)
        fday=f*48
        fint=fday<0.8
        ax1.semilogx(fday[fint], Cxy[fint],label=str(nn)+' db',color=cm.YlGnBu(ii*30+80))
ax1.legend(loc=(1.02,0.2))
ax1.set_ylabel('with meridional velocity')
for ii,nn in enumerate(nomprsvec):
        f, Cxy = signal.coherence(pb, adat['v'].sel(nomprs=nn),nperseg=4096)
        fday=f*48
        fint=fday<0.8
        ax2.plot(fday[fint], Cxy[fint],label=str(nn)+' db',color=cm.YlGnBu(ii*30+80))
ax2.set_ylabel('with zonal velocity')
ax1.set_title(moorname+' coherence between bottom pressure and velocity',fontsize='x-large')
ax2.set_xlabel('CPD')
plt.tight_layout()
savefig('../figures/VertModes/from_raw_data/Coherence_'+moorname+'.png',bbox_inches='tight')
