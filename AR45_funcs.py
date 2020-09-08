from pylab import *
import gsw
from AAoptode_ctdconvert_funcs import *

from seabird import *
import glob
from scipy import io,interpolate
import seawater as sw
import xarray as xr

figdir='/home/isabela/Documents/projects/OSNAP/figures/AR45/'

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

uni['o2']={}
uni['o2']['vmin']=270
uni['o2']['vmax']=325
uni['o2']['cmap']=cm.rainbow

uni['odiff']={}
uni['odiff']['vmin']=-20
uni['odiff']['vmax']=20
uni['odiff']['cmap']=cm.RdBu_r

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
