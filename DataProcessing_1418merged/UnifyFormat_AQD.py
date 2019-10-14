from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'

# In this script I want to assess the current status and lay out the necessary steps to acheive uniformity.

# What does M1 2016-2018 data look like? Might be what I'm trying to match.
m1dir18=glob.glob(datadir+'OSNAP2018recovery/M1_netcdf/*Nortek*')
m1_18=xr.open_dataset(m1dir18[0])

m1dir16=glob.glob(datadir+'OSNAP2016recovery/M1_netcdf/*Nortek*')
m1_16=xr.open_dataset(m1dir16[0])


def quickpresplot(xray):
    figure(figsize=(14,6))
    plot(xray.TIME,xray.PRES);

#data return is every 0.5 hour
#here is a (2nd order butterworth) 40hour lowpass filter:
def lowpassfilt(xray,tstep):
    xray_sm=NaN*xray
    Z,X = sig.butter(2,tstep/40, output='ba')
    for var in ['PRES','UCUR','VCUR']:
        for ii,dd in enumerate(xray.DEPTH):
            nanind=~isnan(xray[var][:,ii])
            xray_sm[var][nanind,ii]=sig.filtfilt(Z,X,xray[var][nanind,ii])
    return xray_sm

m1_18_sm=lowpassfilt(m1_18,0.5)
m1_16_sm=lowpassfilt(m1_16,1)

m1_18.VCUR[:,-1].plot(figsize=(14,3))
m1_18_sm.VCUR[:,-1].plot()

m1_18_daily=m1_18_sm.resample(TIME='1D').mean(dim='TIME')
m1_16_daily=m1_16_sm.resample(TIME='1D').mean(dim='TIME')
