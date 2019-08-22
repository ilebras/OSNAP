from __future__ import division
import numpy
from matplotlib import pyplot

import pycwt as wavelet
from pycwt.helpers import find

import glob as glob

import xarray as xr

import pandas as pd

figdir='/home/isabela/Documents/projects/OSNAP/figures/wavelet/'

datadir='/home/isabela/Documents/projects/OSNAP/data/OSNAP2016recovery/'


from pylab import *

## Tidal filtering was BOLLOCKS
# import ttide as tt
# date[0].values
#
#
# def ttide_filt(utest,vtest):
#     #Subtract TIDES
#     if sum(~isnan(utest))<3:
#         ufilt=utest
#         vfilt=vtest
#     else:
#         tfit_u = tt.t_tide(utest+1j*vtest,dt=0.5,lat=60);
#         tidepart=tfit_u(arange(0,len(utest)/2,0.5));
#         ufilt=utest-tidepart.real;
#         vfilt=vtest-tidepart.imag;
#     return ufilt,vfilt
#
# ufilt,vfilt=ttide_filt(dataset['UCUR'].values.flatten(),dataset['VCUR'].values.flatten())


title = 'Current meter kinetic energy'
label = 'Kinetic energy'
units = '$m^2 /s^2$'


# don't think I need to de-trend...

## Define wavelet parameters
mother = wavelet.Morlet(6)
s0 = 1  # Starting scale, in hours
dj = 1 / 48  # Number of sub-octaves per octaves
J = 12 / dj  # Number powers of two with dj sub-octaves


def AQDlist(moornum):
    moorname='CF'+str(moornum)
    if moorname=='CF8':
        aqdlist=hstack((glob.glob(datadir+'NOC_M1/nocm1_01_2014/nor/*.edt'),
                        glob.glob(datadir+'NOC_M1/nocm1_02_2015/nor/*.edt')))
    else:
        aqdlist=glob.glob(datadir+'AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc')

    return aqdlist

aqdlist=AQDlist(7)

def takewav_makefig(dd,moornum):

        if moornum==8:
            dt=1
            dat=pd.read_csv(dd,header=12,sep='\s*')
            date=unique([datetime.datetime(int(dat.ix[ii,0]),
                                               int(dat.ix[ii,1]),
                                               int(dat.ix[ii,2]),
                                               int(dat.ix[ii,3])) for ii in range(len(dat))])
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            nomd=int(nanmean(array(dat.ix[:,5])))
            dat=utest**2+vtest**2
            savetit='M1-'+str(nomd)+'m'

        else:
            dataset=xr.open_dataset(dd)
            date=dataset['TIME']
            ke=dataset['UCUR']**2+dataset['VCUR']**2
            dat=ke.values.flatten()
            dt=0.5
            nomd=int(dataset.geospatial_vertical_min)
            savetit=dataset.platform_code[-3:]+'-'+str(nomd)+'m'

        dat[isnan(dat)]=nanmean(dat)
        alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise

        N=len(dat)
        #in hours

        t = numpy.arange(0, N) * dt

        std=dat.std()
        var=std**2
        dat_norm=dat/std


        # The following routines perform the wavelet transform and inverse wavelet transform using the parameters defined above. Since we have normalized our input time-series, we multiply the inverse transform by the standard deviation.
        wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm,dt, dj, s0, J, mother)

        iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

        # We calculate the normalized wavelet and Fourier power spectra, as well as the Fourier equivalent periods for each wavelet scale.

        power = (numpy.abs(wave)) ** 2
        fft_power = numpy.abs(fft) ** 2
        period = 1 / freqs

        # Optionally, we could also rectify the power spectrum according to the suggestions proposed by Liu et al. (2007)[2]

        power /= scales[:, None]

        # We could stop at this point and plot our results. However we are also interested in the power spectra significance test. The power is significant where the ratio power / sig95 > 1.


        signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                                 significance_level=0.95,
                                                 wavelet=mother)


        sig95 = numpy.ones([1, N]) * signif[:, None]
        sig95 = power / sig95

        # Then, we calculate the global wavelet spectrum and determine its significance level.

        glbl_power = power.mean(axis=1)
        dof = N - scales  # Correction for padding at edges
        glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha,
                                                significance_level=0.95, dof=dof,
                                                wavelet=mother)

        # We also calculate the scale average between pmin and pmax, and its significance level.
        f,dx = pyplot.subplots(6,1,figsize=(12,12),sharex=True)
        bands=[1,2,8,16,48,128,512]
        for ii in range(len(bands)-1):
            pmin=bands[ii]
            pmax=bands[ii+1]
            sel = find((period >= pmin) & (period < pmax))
            Cdelta = mother.cdelta
            scale_avg = (scales * numpy.ones((N, 1))).transpose()
            scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
            scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
            scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha,
                                                         significance_level=0.95,
                                                         dof=[scales[sel[0]],
                                                              scales[sel[-1]]],
                                                         wavelet=mother)

            dx[ii].axhline(scale_avg_signif, color='C'+str(ii), linestyle='--', linewidth=1.)
            dx[ii].plot(date, scale_avg, '-', color='C'+str(ii), linewidth=1.5,label='{}--{} hour band'.format(pmin,pmax))
            [dx[ii].axvline(dd,color=clist[jj],linewidth=3) for jj,dd in enumerate(dlist)]
            dx[ii].legend()

        dx[0].set_title('Scale-averaged power: '+savetit)
        dx[3].set_ylabel(r'Average variance [{}]'.format(units))
        if moornum ==8:
            dx[0].set_xlim(date[0],date[-1])
        else:
            dx[0].set_xlim(date[0].values,date[-1].values)
        savefig(figdir+'ScaleSep_'+savetit+'.png',bbox_inches='tight')

        pmin=2
        pmax=24


        sel = find((period >= pmin) & (period < pmax))
        Cdelta = mother.cdelta
        scale_avg = (scales * numpy.ones((N, 1))).transpose()
        scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
        scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
        scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha,
                                                     significance_level=0.95,
                                                     dof=[scales[sel[0]],
                                                          scales[sel[-1]]],
                                                     wavelet=mother)



        figprops = dict(figsize=(11, 8), dpi=72)
        fig = pyplot.figure(**figprops)

        # First sub-plot, the original time series anomaly and inverse wavelet
        # transform.
        ax = pyplot.axes([0.1, 0.75, 0.65, 0.2])
        ax.plot(date, dat, linewidth=1.5, color=[0.5, 0.5, 0.5])
        ax.plot(date, iwave, 'k-', linewidth=1,zorder=100)
        if moornum ==8:
            ax.set_xlim(date[0],date[-1])
        else:
            ax.set_xlim(date[0].values,date[-1].values)
        # ax.set_title('a) {}'.format(title))
        ax.set_ylabel(r'{} [{}]'.format(label, units))
        # Second sub-plot, the normalized wavelet power spectrum and significance
        # level contour lines and cone of influece hatched area. Note that period
        # scale is logarithmic.

        bx = pyplot.axes([0.1, 0.37, 0.65, 0.28])
        levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
        bx.contourf(t, numpy.log2(period), numpy.log2(power), numpy.log2(levels),
                    extend='both', cmap=pyplot.cm.viridis)
        extent = [t.min(), t.max(), 0, max(period)]
        bx.contour(t, numpy.log2(period), sig95, [-99, 1], colors='k', linewidths=2,
                   extent=extent)
        bx.fill(numpy.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                                   t[:1] - dt, t[:1] - dt]),
                numpy.concatenate([numpy.log2(coi), [1e-9], numpy.log2(period[-1:]),
                                   numpy.log2(period[-1:]), [1e-9]]),
                'k', alpha=0.3, hatch='x')
        bx.set_title('{} Wavelet Power Spectrum ({})'.format(label, mother.name))
        bx.set_ylabel('Period (hours)')
        #
        Yticks = 2 ** numpy.arange(numpy.ceil(numpy.log2(period.min())),
                                   numpy.ceil(numpy.log2(period.max())))
        bx.set_yticks(numpy.log2(Yticks))
        bx.set_yticklabels(Yticks)
        bx.set_xticklabels('')
        bx.set_xlim(t.min(),t.max())


        # Third sub-plot, the global wavelet and Fourier power spectra and theoretical
        # noise spectra. Note that period scale is logarithmic.
        cx = pyplot.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
        cx.plot(glbl_signif, numpy.log2(period), 'k--')
        cx.plot(var * fft_theor, numpy.log2(period), '--', color='#cccccc')
        cx.plot(var * fft_power, numpy.log2(1./fftfreqs), '-', color='#cccccc',
                linewidth=1.)
        cx.plot(var * glbl_power, numpy.log2(period), 'k-', linewidth=1.5)
        cx.set_title('Global Wavelet Spectrum')
        cx.set_xlabel(r'Power [({})^2]'.format(units))
        cx.set_xlim([0, glbl_power.max() + var])
        cx.set_ylim(numpy.log2([period.min(), period.max()]))
        cx.set_yticks(numpy.log2(Yticks))
        cx.set_yticklabels(Yticks)
        pyplot.setp(cx.get_yticklabels(), visible=False)

        spowdic={}
        spowdic['sig']=scale_avg_signif
        if moornum==8:
            spowdic['date']=date
        else:
            spowdic['date']=date.values
        spowdic['spow']=scale_avg

        # Fourth sub-plot, the scale averaged wavelet spectrum.
        dx = pyplot.axes([0.1, 0.07, 0.65, 0.2], sharex=ax)
        dx.axhline(scale_avg_signif, color='k', linestyle='--', linewidth=1.)
        dx.plot(date, scale_avg, 'k-', linewidth=1.5)
        dx.set_title('{}--{} hour scale-averaged power'.format(pmin,pmax))
        # [dx.axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]
        # dx.set_xlabel('Time (hours)')
        dx.set_ylabel(r'Average variance [{}]'.format(units))
        if moornum ==8:
            dx.set_xlim(date[0],date[-1])
        else:
            dx.set_xlim(date[0].values,date[-1].values)

        fig.suptitle(savetit)
        savefig(figdir+'Wavelet_'+savetit+'.png',bbox_inches='tight')

        return nomd,spowdic

clist=['#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b']
dlist=['2014-9-20','2015-2-20','2015-4-20','2015-9-20','2016-2-1','2016-4-1']

spowdic={}
for moornum in range(4,9):
    aqdlist=AQDlist(moornum)
    spowdic[moornum]={}
    for dd in aqdlist:
        nomd,spd_tmp=takewav_makefig(dd,moornum)
        spowdic[moornum][nomd]=spd_tmp



import pickle
pickle.dump(spowdic,open(datadir+'pickles/spectral/scale_avg_power_2-24hour.pickle','wb'),protocol=2)


def spowcomp():
    for moor in [5,6,7,8]:
        figure(figsize=(12,2))
        for ii,key in enumerate(spowdic[moor]):
            plot(spowdic[moor][key]['date'],spowdic[moor][key]['spow'],label=key,color='C'+str(ii),alpha=0.5)
            axhline(spowdic[moor][key]['sig'],color='C'+str(ii),label='',linestyle='--',alpha=0.5)
            [axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]
        legend(loc=(1.01,0))
        if (moor==5) | (moor==6):
            ylim(0,0.006)
        else:
            ylim(0,0.0015)
        pyplot.title('CF'+str(moor))
        savefig(figdir+'Spowcomp_CF'+str(moor)+'.png',bbox_inches='tight')


def spowcomp():
    figure(figsize=(18,3))
    for moor in [5,6,7,8]:
        for ii,key in enumerate(spowdic[moor]):
            if (key==500) | (key==516) | (key==549):
                plot(spowdic[moor][key]['date'],spowdic[moor][key]['spow'],label='CF'+str(moor),alpha=0.75)
                # axhline(spowdic[moor][key]['sig'],color='C'+str(ii),label='',linestyle='--',alpha=0.5)
                [axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]
        legend(loc=(1.01,0))
        ylim(0,0.003)
        savefig(figdir+'Spowcomp_500m.png',bbox_inches='tight')



spowcomp()

#### Try just doing the scale-averaged power for a bunch of bands for each instrument to get a better idea for how its distributed
#### Then can compare different instruments
