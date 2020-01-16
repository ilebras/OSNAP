from pylab import *
get_ipython().magic('matplotlib inline')

from scipy import io,interpolate
from scipy.interpolate import griddata,interp1d
from netCDF4 import Dataset
from dateutil.relativedelta import relativedelta
import pickle
import gsw
import glob
import pandas as pd
import datetime as dt
import xarray as xr
import seawater as sw

from seabird import cnv

datadir='/home/isabela/Documents/projects/OSNAP/data/Shipboard/AR30_2018/'
figdir='/home/isabela/Documents/projects/OSNAP/figures/'

startsec={}
startsec['cal 1']=0
startsec['section 1']=2
startsec['section 1 moorings']=2
startsec['section 2']=35
startsec['section 3']=59
startsec['section 4']=111
startsec['section 5']=113
startsec['section 3b']=148
startsec['section 6']=177
startsec['releases']=[178]
startsec['section 7']=195
startsec['cal 2']=196
startsec['section 8']=198
startsec['cal 3']=219
startsec['section 9']=238
startsec['section 10']=239
startsec['section 11']=273

seclab={}
seclab['cal 1']=[0,1]
seclab['section 1']=hstack((range(2,20),[21,24,25,22,26,55,58]))
seclab['section 1 moorings']=[49,50,55,58]
seclab['section 2']=hstack((range(27,49),range(51,55),range(56,58)))
seclab['section 3']=hstack((range(59,87),range(88,91)))
seclab['section 4']=range(91,112)#remove station 113 (i.e. index 112) strange data issues there
seclab['section 5']=hstack((range(113,132),range(149,154),range(179,182)))
seclab['section 3b']=range(132,149)
seclab['section 6']=range(154,178)
seclab['releases']=[178]
seclab['section 7']=range(182,196)
seclab['cal 2']=[196,197]
seclab['section 8']=range(198,219)
seclab['cal 3']=range(219,222)
seclab['section 9']=range(222,239)
seclab['section 10']=range(239,272)
seclab['section 11']=range(272,296)

def run_ave(vec,rundiv):
    runave=vec.copy()
    for ii in range(int(rundiv/2)):
        runave[ii]=nanmean(vec[:int(ii+rundiv/2)])
    for ii in range(int(rundiv/2),len(vec)-int(rundiv/2)):
        runave[ii]=nanmean(vec[int(ii-rundiv/2):int(ii+rundiv/2)])
    for ii in range(len(vec)-int(rundiv/2),len(vec)):
        runave[ii]=nanmean(vec[int(ii-rundiv/2):])
    # runave=array(runave)
    return runave


def getdist(lon0,lat0,lonex,latex):
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([lat0,latex[ii]],[lon0,lonex[ii]])[0][0]
    return distout


def make_cmap(colors, position=None, bit=False,numbins=256):
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

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,numbins)
    return cmap



colors = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(colors,position=[0,0.9,0.96,0.99,1],bit=True)

colors = [(255,237,160),(127,205,187),(5,112,176),(110,1,107)]
pden_cmap=make_cmap(colors,position=[0,0.666666666666,0.833333333333,1],bit=True)

colors = [(208,28,139),(241,182,218),(184,225,134),(77,172,38)]
turn_map=make_cmap(colors,position=[0,0.33,0.66,1],bit=True,numbins=4)


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
uni['turner']={}
uni['turner']['cmap']=turn_map
uni['turner']['vmin']=-90
uni['turner']['vmax']=90

uni['dTdz']={}
uni['dTdz']['vmin']=-4
uni['dTdz']['vmax']=-1
uni['dTdz']['cmap']=cm.YlGnBu

uni['chi']={}
uni['chi']['vmin']=-12
uni['chi']['vmax']=-5
uni['chi']['cmap']=cm.rainbow

uni['KT']={}
uni['KT']['vmin']=-9
uni['KT']['vmax']=0
uni['KT']['cmap']=cm.rainbow

uni['N2']={}
uni['N2']['vmin']=-8
uni['N2']['vmax']=-2
uni['N2']['cmap']=cm.YlGn

uni['sal_1m']=uni['sal']
uni['tmp_1m']=uni['tmp']
uni['turner_1m']=uni['turner']


uni['o2']={}
uni['o2']['vmin']=4
uni['o2']['vmax']=9
uni['o2']['cmap']=cm.rainbow

salvec=linspace(31,36,103)
tmpvec=linspace(-3,16,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),-44,59.5)
SA_vec_1000=gsw.SA_from_SP(salvec,1e3*ones(len(salvec)),-44,59.5)

CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
sigma1mat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3
        sigma1mat[jj,ii]=gsw.sigma1(SA_vec[ii],CT_vec[jj])
