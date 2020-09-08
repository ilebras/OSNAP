
# coding: utf-8

# ## Get some useful functions together to use in several scripts

from pylab import *

from scipy import io,interpolate
from scipy.interpolate import griddata,interp1d
from netCDF4 import Dataset
import gsw
import glob
import datetime as dt
import xarray as xr
import seawater as sw
import palettable as pal
import matplotlib.gridspec as gridspec

from scipy import signal
import csv

from seabird import cnv

fcor=gsw.f(60)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

rc('xtick',labelsize='Large')
rc('ytick',labelsize='Large')
rc('axes', labelsize='Large')

datadir='/home/isabela/Documents/projects/OSNAP/data/'
figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'


# theta=arctan(vmean[2]/umean[2]) This is how new theta was determined -- just estimated from depth and time averaged flow at this point.
# Copied value here for now so I don't have to load vel necessarily.
theta=1.164248016423335

##################################################################################
##################################################################################
salvec=linspace(31,36,100)
tmpvec=linspace(-2,10,100)
salmat,tmpmat=meshgrid(salvec,tmpvec)

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),-41,60)
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
sigma1mat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3
        sigma1mat[jj,ii]=gsw.sigma1(SA_vec[ii],CT_vec[jj])


def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap



colors = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]#(237,248,177),,
sal_cmap = make_cmap(colors,position=[0,0.9,0.96,0.99,1],bit=True)#0.9,

colors = [(255,237,160),(127,205,187),(5,112,176),(110,1,107)]
pden_cmap=make_cmap(colors,position=[0,0.666666666666,0.833333333333,1],bit=True)


uni={}
uni['sal']={}
uni['sal']['cmap']=sal_cmap
uni['sal']['vmin']=28
uni['sal']['vmax']=35.2
uni['tmp']={}
uni['tmp']['cmap']=cm.RdYlBu_r
uni['tmp']['vmin']=1
uni['tmp']['vmax']=10
uni['den']={}
uni['den']['cmap']=pden_cmap
uni['den']['vmin']=24.5
uni['den']['vmax']=28


##################################################################################
##################################################################################


years=matplotlib.dates.YearLocator()
months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')


## 2nd order 50 day default filter
import scipy.signal as sig
# Design the Buterworth filter
N  = 2    # Filter order
Wn = 0.02 # Cutoff frequency (50 days)
B, A = sig.butter(N, Wn, output='ba')

#potentially useful ploting scripts

def date_from_matlab(matdate):
    pydate=array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in matdate])
    return pydate


import calendar
def toTime(d):
  return [calendar.timegm(dd.timetuple()) for dd in d]

def np64ToDatetime(DA):
  return [datetime.datetime.utcfromtimestamp((dd-np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for dd in DA]
