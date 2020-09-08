# In this script, going to look at hourly instrument data, but lowpass filtered

from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/Slantwise_Convection/'

# Start from file prepared for Astrid
mcat=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/hourly/CF_M_hourly_2014-2018_MCATonly.nc')
vel=xr.open_dataset(datadir+'OSNAP_CFgridded_2014-2018/hourly/CF_M_hourly_2014-2018_VELonly.nc')


# filter all data with 25 hour cutoff lowpass filter
def lowpassfilt(xray):
    tstep=1
    xray_sm=NaN*xray
    Z,X = sig.butter(2,tstep/25, output='ba')
    for var in xray:
        for ii,zz in enumerate(xray.DEPTH):
            for jj,dd in enumerate(xray.distance):
                nanind=~isnan(xray[var][ii,jj,:])
                if sum(nanind)>10:
                    xray_sm[var][ii,jj,nanind]=sig.filtfilt(Z,X,xray[var][ii,jj,nanind])
                else:
                    xray_sm[var][ii,jj,:]=NaN
    return xray_sm


vel_sm=lowpassfilt(vel)
mcat_sm=lowpassfilt(mcat)

vel_sm.DEPTH.values

vel_sm.depth_var[32,0,:].plot()
