### http://nsidc.org/data/nsidc-0051
###The data can be read with image processing software by specifying a 300-byte header, with an image size of 304 columns x 448 rows for Arctic
###North: 136492 bytes
#### Data are stored as one-byte integers representing sea ice concentration values. The sea ice concentration data are packed into byte format by multiplying the derived fractional sea ice concentration floating-point values (ranging from 0.0 to 1.0) by a scaling factor of 250.

## Scripting help from: https://ocefpaf.github.io/python4oceanographers/blog/2015/04/20/arctic_sea_ice_concentration/

import numpy as np
import numpy.ma as ma

%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
from pyproj import Proj,transform

files=np.sort(glob.glob('../data/seaice/*'))

dx = dy = 25000

x = np.arange(-3850000, +3750000, +dx)
y = np.arange(+5850000, -5350000, -dy)

inProj =Proj(init='epsg:3413')

outProj = Proj(init='epsg:4326')

xmat,ymat=np.meshgrid(x,y)

lonmat,latmat = transform(inProj,outProj,xmat,ymat)

import datetime
datevec=[datetime.datetime(2014,8,15)+datetime.timedelta(days=ii*30) for ii in range(25)]

from aux_funcs import *

minlon=-43
maxlon=-42.5
minlat=60
maxlat=61
xCF,yCF=np.where((lonmat>minlon)&(lonmat<maxlon)&(latmat>minlat)&(latmat<maxlat))
xFS,yFS=np.where((lonmat>-15)&(lonmat<0)&(latmat>78)&(latmat<80))
xDS,yDS=np.where((lonmat>-28)&(lonmat<-24)&(latmat>66)&(latmat<69))

files

files[9]

CFice=np.zeros(len(files))
FSice=np.zeros(len(files))
DSice=np.zeros(len(files))
for ii,infile in enumerate(files):

    with open(infile, 'rb') as fr:
        hdr = fr.read(300)
        ice = np.fromfile(fr, dtype=np.uint8)

    ice = ice.reshape(448, 304)

    ice = ice / 250.

    ice = ma.masked_greater(ice, 1.0)

    CFice[ii]=ma.mean(ice[xCF,yCF])
    DSice[ii]=ma.mean(ice[xDS,yDS])
    FSice[ii]=ma.mean(ice[xFS,yFS])
    # fig = plt.figure(figsize=(9, 9))
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
    #
    # cs = ax.coastlines(resolution='10m', linewidth=0.5)
    #
    #
    # ax.gridlines()
    # ax.set_extent([-90, 10, 60, 90], crs=ccrs.PlateCarree())
    #
    # kw = dict(central_latitude=90, central_longitude=-45, true_scale_latitude=70)
    # cs = ax.pcolormesh(x, y, ice, cmap=plt.cm.Blues,
    #                    transform=ccrs.Stereographic(**kw))
    # yr=infile[18:22]
    # mnth=infile[22:24]
    # ax.set_title(mnth+'/'+yr,fontsize=18)
    # cbaxes = fig.add_axes([0.95, 0.25, 0.04, 0.45])
    # cb=colorbar(cs,cax=cbaxes)
    # cb.set_label('sea ice concentration',fontsize=16)
    # fig.savefig('../../figures/seaice/maps/'+yr+'_'+mnth+'_seaice.png',bbox_inches='tight')
    # fig.savefig('../../figures/seaice/maps/'+yr+'_'+mnth+'_seaice.pdf',bbox_inches='tight')



infile

[cc,egic,eg,ic]=pickle.load(open('../pickles/transdic_1810JHcal.pickle','rb'))
daily=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_1810JHcal.pickle','rb'))

years=matplotlib.dates.YearLocator()
months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('\n %Y')

def compseries(date1,field1,tit1,c1,date2,field2,tit2,c2,ax1=0):
    if ax1==0:
        fig,ax1=subplots(figsize=(10,3))
    ax1.plot(date1,field1,color=c1,linewidth=2)
    ax1.set_ylabel(tit1,color=c1)
    ax2=ax1.twinx()
    ax2.plot(date2,field2,color=c2,linewidth=2)
    ax2.set_ylabel(tit2,color=c2)
    ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_minor_locator(threemonth)
    ax2.xaxis.set_minor_formatter(monthFMT)
    ax2.xaxis.set_major_formatter(yearFMT)


en4=io.loadmat(datadir+'aux_data/EN4_drainage.mat')

def make_paperfig():
    fig=figure(figsize=(8,8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.75, 1])

    ax=subplot(gs[0],projection=ccrs.NorthPolarStereo(central_longitude=0))
    infile=files[9]
    with open(infile, 'rb') as fr:
        hdr = fr.read(300)
        ice = np.fromfile(fr, dtype=np.uint8)
    ice = ice.reshape(448, 304)
    ice = ice / 250.
    ice = ma.masked_greater(ice, 1.0)
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
    cs = ax.coastlines(resolution='10m', linewidth=0.5)
    # ax.gridlines()
    ax.set_extent([-90, 10, 60, 90], crs=ccrs.PlateCarree())
    kw = dict(central_latitude=90, central_longitude=-45, true_scale_latitude=70)
    cs = ax.contourf(x, y, ice,arange(0,1.1,0.1),cmap=plt.cm.Blues,
                       transform=ccrs.Stereographic(**kw))
    ax.contour(x, y, ice,0.2,colors='k',transform=ccrs.Stereographic(**kw))
    lon_start=-43#4
    lon_end=-42#2
    lat_start=60
    lat_end=61
    x_start,y_start=transform(outProj,inProj,lon_start,lat_start)
    x_end,y_end=transform(outProj,inProj,lon_end,lat_end)
    # x_start,y_start=transform(outProj,inProj,lon_start,lat_start)
    # ax.plot([x_start,x_start,x_end,x_end,x_start],[y_start,y_end,y_end,y_start,y_start],linewidth=2,color='r',transform=ccrs.Stereographic(**kw))
    lon_FS1=-18
    lon_FS2=0
    lat_FS=78.5
    x_FS1,y_FS=transform(outProj,inProj,lon_FS1,lat_FS)
    x_FS2,y_FS=transform(outProj,inProj,lon_FS2,lat_FS)
    ax.plot([x_FS1,x_FS2],[y_FS,y_FS],linewidth=4,color='k',transform=ccrs.Stereographic(**kw))
    lon_DS1=-26.5
    lon_DS2=-24
    lat_DS1=68.5
    lat_DS2=67
    x_DS1,y_DS1=transform(outProj,inProj,lon_DS1,lat_DS1)
    x_DS2,y_DS2=transform(outProj,inProj,lon_DS2,lat_DS2)
    ax.plot([x_DS1,x_DS2],[y_DS1,y_DS2],linewidth=4,color='k',transform=ccrs.Stereographic(**kw))
    ax.plot([x_start,x_end],[y_start,y_start],linewidth=4,color='k',transform=ccrs.Stereographic(**kw))
    x_basin,y_basin=transform(outProj,inProj,en4['lon_basin'],en4['lat_basin'])
    ax.plot(x_basin,y_basin,linewidth=2,color='k',transform=ccrs.Stereographic(**kw))
    yr=infile[18:22]
    mnth=infile[22:24]
    # ax.set_title(mnth+'/'+yr,fontsize=14)
    # cbaxes = ax.add_axes([0.95, 0.25, 0.04, 0.45])
    cb=colorbar(cs)#,cax=cbaxes)
    cb.set_label('sea ice concentration\n 05/2015')

    ax1=subplot(gs[1])
    compseries(daily.date,cc['trans filt'],'Coastal current transport [Sv]','k',
                datevec,CFice,'Sea ice concentration','blue',ax1)
    savefig('../figures/paperfigs/seaice_4poster.pdf',bbox_inches='tight')

make_paperfig()


def eachpanel(field,colo,axit,labit='',xr=daily,pnofilt=1,letlab='',ls='-'):
    if pnofilt==1:
        axit.plot(xr.date,field,alpha=0.5,color=colo,label='',linewidth=0.75)
    axit.plot(xr.date,sig.filtfilt(B,A, field),linewidth=2,color=colo,label=labit,linestyle=ls)
    # axit.set_xlabel('')
    # axit.set_ylabel('')
    axit.text(0.01, 0.85,letlab,transform = axit.transAxes,fontsize=15)


def compfresh():
    fig=figure(figsize=(10,3.5))
    ax=subplot(111)
    freshtype='freshb'
    # fill_between(daily.date.values, sig.filtfilt(B,A, cc[freshtype]), sig.filtfilt(B,A,cc[freshtype+' plus']),color=ccol,label='',alpha=0.6)
    eachpanel(cc[freshtype],ccol,ax,pnofilt=0,labit='EGCC (coastal current)')
    eachpanel(eg[freshtype],egcol,ax,pnofilt=0,labit='EGC (slope current, S<34.9)')
    eachpanel(ic[freshtype],icol,ax,pnofilt=0,labit='IC (slope current, S>34.9)')
    axhline(0,color='k')
    # axvline(datetime.datetime(2015,7,1,0),color='grey',linewidth=0.8)
    # pwf(cc[freshtype]+egic[freshtype],'grey',1)
    ylim([-30,120])
    ylabel('Freshwater transport [mSv]')
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,1)])
    colorstripes()
    # fig.autofmt_xdate()
    # dticks=[datetime.datetime(2014,10,1),datetime.datetime(2015,1,1),datetime.datetime(2015,4,1),datetime.datetime(2015,7,1),datetime.datetime(2015,10,1),datetime.datetime(2016,1,1),datetime.datetime(2016,4,1),datetime.datetime(2016,7,1)]
    # gca().set_xticks(dticks)
    xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_minor_locator(threemonth)
    ax.xaxis.set_minor_formatter(monthFMT)
    ax.xaxis.set_major_formatter(yearFMT)
    legend(loc=1)
    labyy=90
    text(datetime.datetime(2014,10,15),labyy,'FALL',color=ccol,fontsize=13.5)
    text(datetime.datetime(2015,1,10),labyy,'WINTER',color=egcol,fontsize=13.5)
    text(datetime.datetime(2015,4,18),labyy,'SPRING/SUMMER',color='grey',fontsize=13.5)
    ax2=ax.twinx()
    ax2.plot(datevec,CFice,'k--',label='sea ice concentration')
    ax2.legend(loc=4)
    ax2.set_ylim(-0.15,0.6)
    ax2.set_xticklabels('')
    ax2.set_ylabel('sea ice concentration')
    ax2.set_yticks([0,0.2,0.4])
    # text(datetime.datetime(2015,7,5),75,'SUMMER',color='grey',fontsize=15)
    savefig('../figures/paperfigs/freshcomp_wseaice.pdf',bbox_inches='tight')

compfresh()

figure(figsize=(10,3))
ax1=subplot(111)
compseries(daily.date,cc['trans filt'],'Coastal current transport [Sv]',ccol,
            datevec,CFice,'Sea ice concentration','k',ax1)
savefig('../figures/paperfigs/seaice_CC.pdf',bbox_inches='tight')

plt.figure(figsize=(12,3))
plt.plot(datevec,FSice/max(FSice),label='Fram Strait')
plt.plot(datevec,DSice/max(DSice),label='Denmark Strait')
plt.plot(datevec,CFice/max(CFice),label='Cape Farewell')
plt.ylabel('normalized sea ice concentration')
plt.legend()
plt.savefig('../../figures/seaice/EG_seaicenorm_nozero.png')

plt.figure(figsize=(12,3))
plt.plot(datevec,FSice,label='Fram Strait')
plt.plot(datevec,DSice,label='Denmark Strait')
plt.plot(datevec,CFice,label='Cape Farewell')
plt.ylabel('sea ice concentration')
plt.legend()
plt.savefig('../../figures/seaice/EG_seaice.png')

plt.figure(figsize=(12,3))
plt.plot(datevec,DSice,label='Denmark Strait')
plt.plot(datevec,CFice,label='Cape Farewell')
plt.ylabel('sea ice concentration')
plt.legend()
plt.savefig('../../figures/seaice/DSCF_seaice.png')

plt.figure(figsize=(12,3))
plt.plot(CFice[:12],label='2014-2015')
plt.plot(CFice[12:-1],label='2015-2016')
plt.xticks(range(0,13,3))
plt.legend()
plt.gca().set_xticklabels(['August','November','February','May'])
plt.title('sea ice concentration at Cape Farewell')
plt.savefig('../../figures/seaice/CF_yearcomp_seaice_nozero.png')

plt.figure(figsize=(12,3))
plt.plot(datevec,CFice)
plt.title('sea ice concentration at Cape Farewell')
plt.savefig('../../figures/seaice/CF_seaice_nozero.png')


#######################################
## Monthly mean surface air temperature Eastern Greenland https://data.giss.nasa.gov/cgi-bin/gistemp/stdata_show.cgi?id=431043600000&dt=1&ds=5
#######################################
#August 2014-August 2016 inclusive:
from pylab import *
Dan=[3.70 ,  -2.80 , NaN,  -17.70,  NaN, -14.60 , -25.30 , NaN,  NaN,  -11.10  ,  2.70 ,   4.70  ,  3.00 ,  -2.90,  -10.30,  -15.00,  -20.90, -18.10,  -21.90  ,-19.90, -17.90 ,  -4.30 ,   2.70 ,   6.00 ,   6.00 ] #(76.77N, 18.67W) Danmarkshavn
Ang=[NaN,NaN, -0.60,   -1.40,   -4.70,-5.30 ,  -8.80  , -3.70 ,  -3.10 ,   0.80 ,   4.70  ,  7.00  ,  6.50 ,   4.80  ,  0.20  , -3.30 ,  -4.90 ,-2.50  , -4.30 ,  -3.00  ,  0.60  ,  1.80  ,  6.80 ,   9.10  , 8.60 ] #(65.60N, 37.63W) Angmassalik
Chr=[8.70  ,  5.70  ,  2.30  , -0.20 ,  -2.80 ,-4.40,   -5.40 ,  -4.80 , -1.30  ,  1.50,   4.30,    5.80 ,  7.70  ,  5.10  ,  1.60,  -1.00 ,  -3.60 ,-3.20,   -2.60,   -2.30 ,   1.50,    2.50,    5.80 ,   7.30  ,  8.00] #(60.05N, 43.17W) Prins Christi
len(Dan)

figure(figsize=(12,3))
plot(datevec,Dan, label='Danmarkshavn, 77$^\circ$N')
plot(datevec,Ang, label='Angmassalik, 66$^\circ$N')
plot(datevec,Chr, label='Prins Christi, 60$^\circ$N')
legend()
axhline(-2,color='k')
ylabel('Surface Air Temperature ($^\circ$ C)')
savefig('../../figures/seaice/surfaceairtmp_EG.png')

figure(figsize=(12,3))
plot(datevec,Chr)
axhline(-2,color='k')
title('Surface Air Temperature at Prins Christi station, 60$^\circ$N [$^\circ$ C]')
savefig('../../figures/seaice/surfaceairtmp_CF.png')

fig, ax1 = plt.subplots()
t = np.arange(0.01, 10.0, 0.01)
s1 = np.exp(t)
ax1.plot(t, s1, 'b-')
ax1.set_xlabel('time (s)')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('exp', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
s2 = np.sin(2 * np.pi * t)
ax2.plot(t, s2, 'r.')
ax2.set_ylabel('sin', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()




compseries(daily.date,cc['trans filt'],'Coastal current transport [Sv]','b',
            datevec,Chr,'Surface Air Temperature [$^\circ$ C]','r')
compseries(datevec,CFice,'Sea ice concentration','r',
            daily.date,egic['trans filt'],'Coastal current transport [Sv]',egicol)

compseries(datevec,CFice,'Sea ice concentration','grey',
            datevec,Chr,'Surface Air Temperature [$^\circ$ C]','r')
