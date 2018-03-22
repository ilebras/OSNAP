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

moornum=5

moorname='CF'+str(moornum)


aqdlist=sort(glob.glob(datadir+'AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc'))

aqdlist
aqd={}
aqd['u']={}
aqd['v']={}
aqd['p']={}
aqd['date']={}
for ii,dd in enumerate(aqdlist):
        dat = Dataset(dd, 'r')
        thekey=int(mean(list(dat.variables['PRES'][:].flatten())))
        u=[float(tt) for tt in dat.variables['UCUR'][:].flatten()]
        v=[float(tt) for tt in dat.variables['VCUR'][:].flatten()]
        p=[float(tt) for tt in dat.variables['PRES'][:].flatten()]
        time_tmp=list(dat.variables['TIME'][:])
        #interpolate over nans
        nanind=where(isnan(u)==False)[0]
        ufunc=interpolate.interp1d([time_tmp[nn] for nn in nanind],[u[nn] for nn in nanind],bounds_error='False')
        vfunc=interpolate.interp1d([time_tmp[nn] for nn in nanind],[v[nn] for nn in nanind])
        pfunc=interpolate.interp1d([time_tmp[nn] for nn in nanind],[p[nn] for nn in nanind])
        aqd['u'][thekey]=ufunc(time_tmp)
        aqd['v'][thekey]=vfunc(time_tmp)
        aqd['p'][thekey]=pfunc(time_tmp)
        # mprs_tmp=nanmean(prs_tmp)*ones(len(time_tmp))
        aqd['date'][thekey]=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=tt)
                        for tt in time_tmp])

nomprsvec=[int(mean(aqd['p'][key])) for key in aqd['p']]

adat=xr.Dataset({'u': (['nomprs', 'date'],  array([aqd['u'][key] for key in aqd['p']])),
                'v': (['nomprs', 'date'],  array([aqd['v'][key] for key in aqd['p']])),
                'prs': (['nomprs', 'date'],  array([aqd['p'][key] for key in aqd['p']])),},
                coords={'nomprs': nomprsvec,
                        'date': aqd['date'][nomprsvec[0]]})

N  = 2    # Filter order
Wn = 0.01 # Cutoff frequency
B, A = signal.butter(N, Wn, output='ba')

#Find correlation between bottom pressure and all velocities
figure(figsize=(12,3))
adat['prs'].sel(nomprs=max(nomprsvec)).plot()
plot(adat.date,signal.filtfilt(B,A, adat['prs'].sel(nomprs=max(nomprsvec))))

for nn in nomprsvec:
    figure(figsize=(12,3))
    adat['v'].sel(nomprs=nn).plot()
    plot(adat.date,signal.filtfilt(B,A, adat['v'].sel(nomprs=nn)))


for nn in nomprsvec:
    print('The following are for: '+str(nn))
    # print('u to u: ')
    # print(corrcoef(adat['u'].sel(nomprs=max(nomprsvec)),adat['u'].sel(nomprs=nn))[0,1])
    # print('v to v: ')
    # print(corrcoef(adat['v'].sel(nomprs=max(nomprsvec)),adat['v'].sel(nomprs=nn))[0,1])
    # print('prs to u: ')
    # print(corrcoef(adat['prs'].sel(nomprs=max(nomprsvec)),adat['u'].sel(nomprs=nn))[0,1])
    print('prs to v: ')
    print(corrcoef(adat['prs'].sel(nomprs=max(nomprsvec)),adat['v'].sel(nomprs=nn))[0,1])
    # print(int(sum(isnan(adat['prs'].sel(nomprs=nn)))))
    # figure(figsize=(12,6))
    # subplot(211)
    # plot(adat['u'].sel(nomprs=nn))
    # subplot(212)
    # plot(adat['v'].sel(nomprs=nn))

from scipy import signal

adat

x=adat['prs'].sel(nomprs=max(nomprsvec))
figure(figsize=(8,6))
for nn in nomprsvec:
    f, Cxy = signal.coherence(x, adat['v'].sel(nomprs=nn))
    semilogy(f, Cxy,label=str(nn)+' db')
xlabel('frequency [1/30min]')
ylabel('Coherence')
title('Coherence between v and bottom pressure')
legend(loc=(1.05,0.4))


x=adat['prs'].sel(nomprs=max(nomprsvec))
figure(figsize=(8,6))
for nn in nomprsvec:
    f, Cxy = signal.coherence(x, adat['v'].sel(nomprs=nn))
    semilogy(f, Cxy,label=str(nn)+' db')
xlabel('frequency [1/30min]')
ylabel('Coherence')
title('coherence between v and bottom pressure')
xlim([0,0.03])
ylim([0.02,1])
legend(loc=(1.05,0.4))

figure(figsize=(8,6))
for nn in nomprsvec:
    f, Cxy = signal.coherence(x, adat['u'].sel(nomprs=nn))
    semilogy(f, Cxy,label=str(nn)+' db')
xlabel('frequency [1/30min]')
ylabel('Coherence')
title('Coherence between u and bottom pressure')
legend(loc=(1.05,0.4))

x=adat['v'].sel(nomprs=max(nomprsvec))
figure(figsize=(8,6))
for nn in nomprsvec:
    f, Cxy = signal.coherence(x, adat['v'].sel(nomprs=nn))
    plot(f, Cxy,label=str(nn)+' db')
xlabel('frequency [1/30min]')
ylabel('Coherence')
title('Coherence between v and bottom v')
legend(loc=(1.05,0.4))
xlim([0.03,0.05])
ylim([0.5,1.01])
#
# ~12 hour peak.. tidal

# 1. Use dispersion relation to get horizontal and vertical scales out of the frequency information.
# 
# 2. Is the vertical scale of fall-off consistent with the dropoff of correlation of v and bottom pressure?
#
# 3. Do this pointwise at multiple moorings - perhaps mean flow effects are apparent?

1/0.0438
##### Load ADCP data

## Leave out for now
# dat=io.loadmat(datadir+'ADCP_Data_CF/OSNAP_cf'+str(moornum)+'_Final1_ilebras.mat')
# prsa=dat['z']
# ua=dat['u']
# va=dat['v']
