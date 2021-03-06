
# coding: utf-8

# ## Get some useful functions together to use in several scripts

from pylab import *
# get_ipython().magic('matplotlib inline')

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
import palettable as pal
import matplotlib.gridspec as gridspec
# import cmocean


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
figdir='/home/isabela/Documents/projects/OSNAP/figures/'


##################################################################################
##################################################################################
##### A couple "convection specific" parameters
##################################################################################
##################################################################################

###densities (note: its critical that this vector be the same as in the gridded density products)
denvec=arange(27.5,27.85,0.005)

d1=27.65
d2=27.73
d3=27.77

p1=argmin(abs(denvec-d1))
p2=argmin(abs(denvec-d2))
p3=argmin(abs(denvec-d3))

dbnds=[d1,d2,d3]


uppercol='#ae017e'
deepcol='#1c9099'

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

egcol='#33a02c'
ccol='#1f78b4'
icol='#e31a1c'
egicol='purple'



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

def plotcontour(field,colormap,vmi,vma,moorname):
    figure(figsize=(20,8))
    contourf(field.columns,field.index,field.values,50,cmap=colormap,
            vmin=vmi,vmax=vma)
    gca().invert_yaxis()
    colorbar()
    xlabel('date')
    ylabel('pressure (db)')
    title(moorname,fontsize=25)



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

def date_from_matlab(matdate):
    pydate=array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in matdate])
    return pydate


# Below is really just a boxcar averager for any time increment
import calendar
def toTime(d):
  return [calendar.timegm(dd.timetuple()) for dd in d]

def np64ToDatetime(DA):
  return [datetime.datetime.utcfromtimestamp((dd-np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for dd in DA]


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


def colorstripes():
    axvspan(datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),color=ccol,alpha=0.4)
    # axvspan(datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,9,1),datetime.datetime(2015,12,1),color=ccol,alpha=0.4)
    axvspan(datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),color=egcol,alpha=0.4)
    # axvspan(datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),color=egcol,alpha=0.4)
    axvspan(datetime.datetime(2015,12,1),datetime.datetime(2016,3,1),color=egcol,alpha=0.4)

#standard averaging constant for turning 15 min microcat into daily.
aveconst=4*24
hrcon=24

def hrly_ave(vec,hourdiv):
    newlen=int(len(vec)/hourdiv+1)
    hrly_vec=zeros(newlen)
    for ii in range(newlen):
        hrly_vec[ii]=nanmean(vec[int(hourdiv*ii):int(hourdiv*(ii+1))])
    hrly_vec=array(hrly_vec)
    return hrly_vec



def interp_overnan(vec):
    fi=interpolate.interp1d(time_hrly[~isnan(vec)],vec[~isnan(vec)],bounds_error=False)
    interpvec=fi(time_hrly)
    return interpvec


def removeline(tvec,svec):
    linefunc=poly1d(polyfit(tvec,svec, 1))
    liney=linefunc(tvec)
    corrected=svec-liney+liney[0]
    return liney,corrected


def pivspline_ML(afield,datevec,prsvec,pdall):
    panpiv=pd.DataFrame(index=datevec,
                        columns=prsvec)
    for ii,dd in enumerate(datevec):
        ind=pdall['date bin']==dd
        x=pdall['pressure'][ind]
        y=pdall[afield][ind]
        Y=array([Y for X, Y in sorted(zip(x, y))])
        if sum(isnan(Y))!=len(Y):
            f1=interpolate.interp1d(hstack((0,sort(x)[~isnan(Y)])),
                                    hstack((Y[~isnan(Y)][0],Y[~isnan(Y)])),
                                    kind='linear',
                                    fill_value="extrapolate")
            pani=f1(prsvec)
            panpiv.iloc[ii]=pd.to_numeric(pani)
    return panpiv.T



maxinstdpth=io.loadmat('../data/OSNAP2016recovery/maxinstdpth.mat')

adcpdp=[160,160,160,340,90,90,90]

def plotinstpos(axchoice,savename):
        if ('uac' in savename) | ('v' in savename):
            for rr in range(7):
                if rr==0:
                    axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,zorder=40,label='ADCP')
                else:
                    axchoice.plot(distvec[rr],adcpdp[rr],'r^',markersize=12,zorder=40)
        mm=0
        for key in depths:
            for dd in range(len(depths[key])):
                if ('sal' in savename) | ('tmp' in savename) | ('pden' in savename):
                    if ('MC' in inst[key][dd]) | ('CTD' in inst[key][dd])  | ('XR-420' in inst[key][dd]):
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=4)

                elif 'uac' in savename:
                    if ('AQ' in inst[key][dd]):
                        if dd==0:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6,label='Current meter')
                        else:
                            axchoice.plot(distvec[mm],depths[key][dd],'ko',zorder=35,markersize=6)
            mm+=1

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





info=io.loadmat(datadir+'OSNAP2016recovery/allmoorinfo.mat')
CFlon=[info['CF'+str(ii)+'_header'][0][0][3][0][0] for ii in range(1,8)]
CFlon[3]=-42.20529938
CFlat=[info['CF'+str(ii)+'_header'][0][0][1][0][0] for ii in range(1,8)]
#Add on M1 location
CFlon=hstack((CFlon,-41-8.409/60))
CFlat=hstack((CFlat,59+54.244/60))
distvec=cumsum(hstack((0,sw.dist(CFlat,CFlon)[0])))

MMlon=[-41.1118,-40.6905,-40.2772,-38.5665,-37.7995]
MMlat=[59.9030,59.8597,59.8145,59.646,59.5788]


###OOI locations
[ooi_lat,ooi_lon]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_locs.pickle','rb'))

ooi_dist=round(sw.dist([CFlat[0],ooi_lat['prf']],[CFlon[0],ooi_lon['prf']])[0][0])
oom_dist=round(sw.dist([CFlat[0],ooi_lat['fla']],[CFlon[0],ooi_lon['fla']])[0][0])


lon0=60

# xdist=sw.dist([CFlat[-2],CFlat[-2]],[CFlon[0],CFlon[-2]])[0][0]
# ydist=sw.dist([CFlat[0],CFlat[-2]],[CFlon[-2],CFlon[-2]])[0][0]
 # theta=abs(arctan((xdist)/(ydist)))


# theta=arctan(vmean[2]/umean[2]) This is how new theta was determined -- just estimated from depth and time averaged flow at this point.
# Copied value here for now so I don't have to load vel necessarily.
theta=1.164248016423335


depths={}
for mm in range(1,8):
    depths['CF'+str(mm)]=ones(len(info['CF'+str(mm)+'_sensors'][0]))
    for ii in range(len(info['CF'+str(mm)+'_sensors'][0])):
        if info['CF'+str(mm)+'_sensors'][0][ii][1][0][0]=='B':
            depths['CF'+str(mm)][ii]=info['CF'+str(mm)+'_header']['depth_actual'][0][0][0][0]
        else:
            depths['CF'+str(mm)][ii]=array(info['CF'+str(mm)+'_sensors'][0][ii][1][0][0])



inst={}
for mm in range(1,8):
    inst['CF'+str(mm)]=[info['CF'+str(mm)+'_sensors'][0][ii][2][0] for ii in range(len(info['CF'+str(mm)+'_sensors'][0]))]

inst['CF5'][0]='lost'
depths['M1']=array([50,50,290,530,530,765,1000,1000,1250,1500,
                    1500,1730,1730,1950,1950])
inst['M1']=array(['AQ','MC','MC','AQ','MC','MC','AQ','MC','MC','AQ',
                  'MC','AQ','MC','MC','WH'])

bathy=io.loadmat(datadir+'k3_clean.mat')
bathdist=hstack((-5,bathy['dist'][:1000].flatten()))
bathbath=hstack((bathy['bath'][0],bathy['bath'][:1000].flatten()))

fullbathdist=bathy['dist'].flatten()
fullbathbath=bathy['bath'].flatten()

osnap_bathy=io.loadmat(datadir+'OSNAPbathy_Fli.mat')

# ## This is a quick rewrite of Dan Torres' matlab script con2sal.m

# % con2sal    Salinity from conductivity, T, P
# %=========================================================================
# % Based on sw_salt by Phil Morgan
# %        $Revision: 1.2 $  $Date: 1994/01/25 05:30:16 $
# %         Copyright (C) CSIRO, Phil Morgan 1993.
# %
# % First convert conductivity to conductivity ratio, then use
# %  the guts of sw_salt to compute salinity
# %
# % CON2SAL.M : Daniel Torres, WHOI, 11/4/96
# %
# % USAGE: S = con2sal(con,T,P) when con = mmho/cm
# % or
# % USAGE: S = con2sal(con.*10,T,P) when con = S/m
# %
# % DESCRIPTION:
# %   Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.
# %
# % INPUT:
# %   con  = conductivity (mmho/cm)
# %   T    = temperature [degree C (IPTS-68)]
# %   P    = pressure    [db]
# %
# % Convert to conductivity ratio by dividing conductivity by 42.9140
# %
# %   cndr = Conductivity ratio     R =  C(S,T,P)/C(35,15,0) [no units]
# %
# % OUTPUT:
# %   S    = salinity    [psu      (PSS-78)]
# %
# % REFERENCES:
# %    Fofonff, P. and Millard, R.C. Jr
# %    Unesco 1983. Algorithms for computation of fundamental properties of
# %    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
# %=========================================================================
#
def con2sal(c,t,p):
    import seawater as sw

    R=c/42.913999999999
    rt=sw.salrt(t)
    rp=sw.salrp(R,t,p)
    ratio=R/(rp*rt)
    S = sw.sals(ratio,t)

    return S

# Parula map

from matplotlib.colors import LinearSegmentedColormap

cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905],
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143],
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952,
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286],
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238,
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571],
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571,
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429],
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667,
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286],
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571,
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429],
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524,
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048,
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667],
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381,
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381],
 [0.0589714286, 0.6837571429, 0.7253857143],
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429],
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429,
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048],
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619,
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667],
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524,
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905],
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476,
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143],
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333],
 [0.7184095238, 0.7411333333, 0.3904761905],
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667,
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762],
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217],
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857,
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619],
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857,
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381],
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857],
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309],
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333,
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333],
 [0.9763, 0.9831, 0.0538]]

parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
