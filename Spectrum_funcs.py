from pylab import *
import scipy.stats as ss
from numpy import fft


#from pycurrents.num import spectra

def spectrum1(h, dt=1):
    """
    First cut at spectral estimation: very crude.

    Returns frequencies, power spectrum, and
    power spectral density.
    Only positive frequencies between (and not including)
    zero and the Nyquist are output.
    """
    nt = len(h)
    npositive = nt//2
    pslice = slice(1, npositive)
    freqs = fft.fftfreq(nt, d=dt)[pslice]
    ft = fft.fft(h)[pslice]
    psraw = abs(ft) ** 2
    # Double to account for the energy in the negative frequencies.
    psraw *= 2
    # Normalization for Power Spectrum
    psraw /= nt**2
    # Convert PS to Power Spectral Density
    psdraw = psraw * dt * nt  # nt * dt is record length
    return freqs, psraw, psdraw


def spectrum2(h, dt=1, nsmooth=5):
    """
    Add simple boxcar smoothing to the raw periodogram.

    Chop off the ends to avoid end effects.
    """
    freqs, ps, psd = spectrum1(h, dt=dt)
    weights = ones(nsmooth, dtype=float) / nsmooth
    ps_s = convolve(ps, weights, mode='valid')
    psd_s = convolve(psd, weights, mode='valid')
    freqs_s = convolve(freqs, weights, mode='valid')
    return freqs_s, ps_s, psd_s


def detrend(h):
    n = len(h)
    t = arange(n)
    p = polyfit(t, h, 1)
    h_detrended = h - polyval(p, t)
    return h_detrended

def quadwin(n):
    """
    Quadratic (or "Welch") window
    """
    t = arange(n)
    win = 1 - ((t - 0.5 * n) / (0.5 * n)) ** 2
    return win

def spectrum4(h, dt=1, nsmooth=5):
    """
    Detrend and apply a quadratic window.
    """
    n = len(h)

    h_detrended = detrend(h)

    winweights = quadwin(n)
    h_win = h_detrended * winweights

    freqs, ps, psd = spectrum2(h_win, dt=dt, nsmooth=nsmooth)

    # Compensate for the energy suppressed by the window.
    psd *= n / (winweights**2).sum()
    ps *= n**2 / winweights.sum()**2

    return freqs, ps, psd
