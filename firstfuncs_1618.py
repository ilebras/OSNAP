## Redoing useful functions for use in the merged 2016/2018 analysis
## Deleting what I deem unimportant/ potentially confusing
## Can refer to aux_funcs.py for the original version used for 2016 recovery

from pylab import *
from scipy import io,interpolate
from scipy.interpolate import griddata,interp1d
from scipy import signal as sig
from netCDF4 import Dataset
import gsw
import glob
import datetime as dt
import pickle
import seawater as sw
import matplotlib.gridspec as gridspec
import xarray as xr
import pandas as pd # there is some panda and pickle dependency below and in a couple scripts, it would be good to phase this out in favor of xarray and netcdf


datadir='/home/isabela/Documents/projects/OSNAP/data/'
figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'

##################################################################################
##### formatting
##################################################################################

years=matplotlib.dates.YearLocator()
months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

rc('xtick',labelsize='Large')
rc('ytick',labelsize='Large')
rc('axes', labelsize='Large')


def date_from_matlab(matdate):
    pydate=array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in matdate])
    return pydate

##################################################################################
##### bathymetry/distance/direction
##################################################################################

bathy=io.loadmat(datadir+'k3_clean.mat')
bathdist=hstack((-5,bathy['dist'][:1000].flatten()))
bathbath=hstack((bathy['bath'][0],bathy['bath'][:1000].flatten()))

fullbathdist=bathy['dist'].flatten()
fullbathbath=bathy['bath'].flatten()

osnap_bathy=io.loadmat(datadir+'OSNAPbathy_Fli.mat')


info=io.loadmat(datadir+'OSNAP2016recovery/allmoorinfo.mat')
CFlon=[info['CF'+str(ii)+'_header'][0][0][3][0][0] for ii in range(1,8)]
CFlon[3]=-42.20529938
CFlat=[info['CF'+str(ii)+'_header'][0][0][1][0][0] for ii in range(1,8)]
#Add on M1 location
CFlon=hstack((CFlon,-41-8.409/60))
CFlat=hstack((CFlat,59+54.244/60))
distvec=cumsum(hstack((0,sw.dist(CFlat,CFlon)[0])))

# lon0=60

# xdist=sw.dist([CFlat[-2],CFlat[-2]],[CFlon[0],CFlon[-2]])[0][0]
# ydist=sw.dist([CFlat[0],CFlat[-2]],[CFlon[-2],CFlon[-2]])[0][0]
# theta=abs(arctan((xdist)/(ydist)))


# theta=arctan(vmean[2]/umean[2]) This is how new theta was determined -- just estimated from depth and time averaged flow at this point.
# Copied value here for now so I don't have to load vel necessarily.
theta=1.164248016423335


##################################################################################
##### SA_CT_PT
##################################################################################

def add_SA_CT_PT(xray):
    if 'PRES' in list(xray.data_vars):
            PRES_out=xray['PRES']
    else:
            PRES_out=-gsw.p_from_z(xray['DEPTH'])
    SA_out=gsw.SA_from_SP(xray['PSAL'],PRES_out,xray.LONGITUDE,xray.LATITUDE)
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['ASAL']=(('TIME','DEPTH'),SA_out)
    xray['PTMP']=(('TIME','DEPTH'),PT_out)
    xray['CTMP']=(('TIME','DEPTH'),CT_out)
    xray['PDEN']=(('TIME','DEPTH'),PD_out)

##################################################################################
##### fit a sin
##################################################################################


from scipy.optimize import leastsq
def fitsin(t,data,guess_mean,guess_phase,guess_std,guess_period):

    first_guess=guess_std*np.sin(2*pi*(t+guess_phase)/guess_period) +guess_mean

    optimize_func = lambda x: x[0]*np.sin(2*pi*(t+x[1])/guess_period) + x[2] - data

    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

    # tgrid=linspace(t[0],t[-1],100)
    est_period=guess_period

    data_fit = est_std*np.sin(2*pi*(t+est_phase)/est_period) + est_mean

    # date=[]
    # for i in range(len(tgrid)):
    #     date.append(dt.date.fromtimestamp(tgrid[i]))

    return data_fit,est_std,est_period


##################################################################################
##### color defintions and a inst pos plotting (to keep the same as 2016 recovery for now)
##################################################################################


egcol='#33a02c'
ccol='#1f78b4'
icol='#e31a1c'
egicol='purple'


def colorstripes():
    axvspan(datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),color=ccol,alpha=0.4)
    # axvspan(datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,9,1),datetime.datetime(2015,12,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),color=egcol,alpha=0.4)
    # axvspan(datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2015,12,1),datetime.datetime(2016,3,1),color=egcol,alpha=0.4)



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

univec={}
univec['pden']=['potential density',array([26.75,27,27.125,27.25,27.375,27.5,27.55,27.6,27.64,27.68,27.695,27.71,27.725,27.74,27.755,27.77,27.785,27.8,27.815,27.83,27.845,27.86,27.875,27.89]),pden_cmap,[27,27.5,27.68,27.74,27.8,27.86],'[kg/m$^3$]']
univec['uacross']=['across track velocity',arange(0,0.6005,0.025),cm.BuPu,arange(0,0.8,0.1),'[m/s]']
univec['sal']=['salinity',array([33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35]),sal_cmap,array([33,34, 34.8,34.92,34.94,34.96,34.98, 35]),'']
univec['tmp']=['temperature',linspace(-1,8,31),cm.RdYlBu_r,[2,4,5,6,7,8],'[$^\circ$C]']
univec['ualong']=['along track velocity',arange(-0.4,0.401,0.02),cm.RdBu_r,arange(-0.4,0.401,0.1),'[m/s]']
univec['geostrophic velocity']=['geostrophic velocity',arange(-0.6,0.605,0.05),cm.RdBu_r,arange(-0.6,0.605,0.2),'[m/s]']
univec['turner angle']=['turner angle',arange(-90,90,1),cm.RdBu_r,arange(-90,95,45),'$^\circ$']


salvec=linspace(31,36,100)
tmpvec=linspace(-2,10,100)
salmat,tmpmat=meshgrid(salvec,tmpvec)

# SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
# CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
# pdenmat=zeros((shape(salmat)))
# for ii in range(len(salvec)):
#     for jj in range(len(tmpvec)):
#         pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])
# pickle.dump(pdenmat,open('../pickles/aux/pdenmat.pickle','wb'))

# pdenmat=pd.read_pickle(open(datadir+'OSNAP2016recovery/pickles/aux/pdenmat.pickle','rb'))

# # Potential density for properties at 750 - more relevant for IIW
# SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
# CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
# pdenmat2=zeros((shape(salmat)))
# for ii in range(len(salvec)):
#     for jj in range(len(tmpvec)):
#         pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3
# pickle.dump(pdenmat2,open(datadir+'OSNAP2016recovery/pickles/aux/pdenmat2.pickle','wb'))

# pdenmat2=pd.read_pickle(open(datadir+'OSNAP2016recovery/pickles/aux/pdenmat2.pickle','rb'))




maxinstdpth=io.loadmat('../data/OSNAP2016recovery/maxinstdpth.mat')


# Keep these around for now, but may want to rewrite - though probably not worth it
# Make sure to reconsider the instrument placement for the 2018 recovery.
#
#
# depths={}
# for mm in range(1,8):
#     depths['CF'+str(mm)]=ones(len(info['CF'+str(mm)+'_sensors'][0]))
#     for ii in range(len(info['CF'+str(mm)+'_sensors'][0])):
#         if info['CF'+str(mm)+'_sensors'][0][ii][1][0][0]=='B':
#             depths['CF'+str(mm)][ii]=info['CF'+str(mm)+'_header']['depth_actual'][0][0][0][0]
#         else:
#             depths['CF'+str(mm)][ii]=array(info['CF'+str(mm)+'_sensors'][0][ii][1][0][0])
#
# inst={}
# for mm in range(1,8):
#     inst['CF'+str(mm)]=[info['CF'+str(mm)+'_sensors'][0][ii][2][0] for ii in range(len(info['CF'+str(mm)+'_sensors'][0]))]
#
# inst['CF5'][0]='lost'
# depths['M1']=array([50,50,290,530,530,765,1000,1000,1250,1500,
#                     1500,1730,1730,1950,1950])
# inst['M1']=array(['AQ','MC','MC','AQ','MC','MC','AQ','MC','MC','AQ',
#                   'MC','AQ','MC','MC','WH'])

# adcpdp=[160,160,160,340,90,90,90]

# def plotinstpos(axchoice):
#         for rr in range(len(adcpdp_new)):
#             if rr==0:
#                 axchoice.plot(distvec_new[rr],adcpdp_new[rr],'^',color='lightgreen',markersize=8,zorder=40,label='ADCP',mec='k')
#             else:
#                 axchoice.plot(distvec_new[rr],adcpdp_new[rr],'^',color='lightgreen',markersize=8,zorder=40,mec='k')
#         mm=0
#         for key in depths:
#             for dd in range(len(depths[key])):
#                     if ('AQ' in inst[key][dd]):
#                         if (mm==0):
#                             axchoice.plot(distvec_new[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='Current meter')
#                         else:
#                             axchoice.plot(distvec_new[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='')
#
#                     if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
#                         if (mm==0) & (dd==0):
#                             axchoice.plot(distvec_new[mm],depths[key][dd],'o',color='#fee090',zorder=38,markersize=4,mec='k',label='T,S recorder')
#                         else:
#                             axchoice.plot(distvec_new[mm],depths[key][dd],'o',color='#fee090',zorder=38,markersize=4,mec='k',label='')
#
#             axchoice.plot([distvec_new[mm],distvec_new[mm]],[depths[key][0],depths[key][-1]],'k')
#
#             mm+=1
#
#
# def plotinstpos(axchoice,savename):
#         if ('uac' in savename) | ('v' in savename):
#             for rr in range(7):
#                 if rr==0:
#                     axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,zorder=40,label='ADCP')
#                 else:
#                     axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,zorder=40)
#         mm=0
#         for key in depths:
#             for dd in range(len(depths[key])):
#                 if ('sal' in savename) | ('tmp' in savename) | ('pden' in savename):
#                     if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
#                             axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=4)
#
#                 elif 'uac' in savename:
#                     if ('AQ' in inst[key][dd]):
#                         if dd==0:
#                             axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='Current meter')
#                         else:
#                             axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6)
#             mm+=1

# simpler, previous version
# def plotinstpos(axchoice,savename):
#         if 'v' in savename:
#             for rr in range(7):
#                 axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=16,zorder=40)
#         mm=0
#         for key in depths:
#             for dd in range(len(depths[key])):
#                 if ('sal' in savename) | ('tmp' in savename) | ('pden' in savename):
#                     if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
#                         axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35)
#                     if 'tmp' in savename:
#                         if ('tidbit' in inst[key][dd]) | ('olo' in inst[key][dd]):
#                              axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35)
#
#                 elif 'v' in savename:
#                     if ('AQ' in inst[key][dd]):
#                         axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35)
#             mm+=1
