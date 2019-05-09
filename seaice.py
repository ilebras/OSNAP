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



inProj = Proj(init='epsg:3413')
outProj = Proj(init='epsg:4326')
np.shape(x)
xmat,ymat=np.meshgrid(x,y)

lonmat,latmat = transform(inProj,outProj,xmat,ymat)

import datetime
datevec=[datetime.datetime(2014,8,15)+datetime.timedelta(days=ii*30) for ii in range(25)]

minlon=-44
maxlon=-41
minlat=59
maxlat=61
xCF,yCF=np.where((lonmat>minlon)&(lonmat<maxlon)&(latmat>minlat)&(latmat<maxlat))
xFS,yFS=np.where((lonmat>-15)&(lonmat<0)&(latmat>78)&(latmat<80))
xDS,yDS=np.where((lonmat>-28)&(lonmat<-24)&(latmat>66)&(latmat<69))

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

    CFice[ii]=ma.mean(ice[xCF,yCF][ice[xCF,yCF]>0])
    DSice[ii]=ma.mean(ice[xDS,yDS][ice[xDS,yDS]>0])
    FSice[ii]=ma.mean(ice[xFS,yFS][ice[xFS,yFS]>0])
    # fig = plt.figure(figsize=(9, 9))
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
    #
    # cs = ax.coastlines(resolution='10m', linewidth=0.5)
    #
    #
    # ax.gridlines()
    # ax.set_extent([-90, 10, 55, 90], crs=ccrs.PlateCarree())
    #
    # kw = dict(central_latitude=90, central_longitude=-45, true_scale_latitude=70)
    # cs = ax.pcolormesh(x, y, ice, cmap=plt.cm.Blues,
    #                    transform=ccrs.Stereographic(**kw))
    # yr=infile[18:22]
    # mnth=infile[22:24]
    # ax.set_title(mnth+'/'+yr)
    # fig.savefig('../figures/seaice/maps/'+yr+'_'+mnth+'_seaice.png')


plt.figure(figsize=(12,3))
plt.plot(datevec,FSice/max(FSice),label='Fram Strait')
plt.plot(datevec,DSice/max(DSice),label='Denmark Strait')
plt.plot(datevec,CFice/max(CFice),label='Cape Farewell')
plt.ylabel('normalized sea ice concentration')
plt.legend()
plt.savefig('../figures/seaice/EG_seaicenorm_nozero.png')

plt.figure(figsize=(12,3))
plt.plot(datevec,FSice,label='Fram Strait')
plt.plot(datevec,DSice,label='Denmark Strait')
plt.plot(datevec,CFice,label='Cape Farewell')
plt.ylabel('sea ice concentration')
plt.legend()
plt.savefig('../figures/seaice/EG_seaice_nozero.png')

plt.figure(figsize=(12,3))
plt.plot(datevec,DSice,label='Denmark Strait')
plt.plot(datevec,CFice,label='Cape Farewell')
plt.ylabel('sea ice concentration')
plt.legend()
plt.savefig('../figures/seaice/DSCF_seaice_nozero.png')

plt.figure(figsize=(12,3))
plt.plot(CFice[:12],label='2014-2015')
plt.plot(CFice[12:-1],label='2015-2016')
plt.xticks(range(0,13,3))
plt.legend()
plt.gca().set_xticklabels(['August','November','February','May'])
plt.title('sea ice concentration at Cape Farewell')
plt.savefig('../figures/seaice/CF_yearcomp_seaice_nozero.png')

plt.figure(figsize=(12,3))
plt.plot(datevec,CFice)
plt.title('sea ice concentration at Cape Farewell')
plt.savefig('../figures/seaice/CF_seaice_nozero.png')


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
axhline(0,color='k')
ylabel('Surface Air Temperature ($^\circ$ C)')
savefig('../figures/seaice/surfaceairtmp_EG.png')
