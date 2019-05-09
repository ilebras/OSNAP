from aux_funcs  import *

predir='/home/isabela/Documents/projects/bathymetry/etopo1/'



plot(bathdist,-bathbath)
"""Get the etopo1 data"""
etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
etopo1 = Dataset(etopo1name,'r')

lons = etopo1.variables["x"][:]
lats = etopo1.variables["y"][:]

etopo_bath = etopo1.variables["z"]

shape(bathy)
shape(lats)
shape(lons)

## extract bathy from etopo1 along CF mooring line... check against knudsen, then extend

def getdepth(lonpt,latpt):
    lonind=argmin(abs(lons-lonpt))
    latind=argmin(abs(lats-latpt))
    bathout=etopo_bath[latind,lonind]

    return bathout

CFdpth=[getdepth(CFlon[cc],CFlat[cc]) for cc in range(8)]

figure()
plot(bathdist,-bathbath)
plot(distvec,CFdpth,'o')
title('Etopo and Knudsen bathymetry check out')


mean(diff(CFlon))
mean(diff(CFlat))


newlats=[CFlat[0]-mean(diff(CFlat))/4*ii for ii in range(1,10)]
newlons=[CFlon[0]-mean(diff(CFlon))/4*ii for ii in range(1,10)]

diff(lats)

figure()
plot(bathy['knu_proj_lon'],bathy['knu_proj_lat'],'o-')
plot(CFlon,CFlat,'o')
plot(newlons,newlats,'o')

newdpth=[getdepth(newlons[cc],newlats[cc]) for cc in range(len(newlats))]


def getdist(lonex,latex):
    # get distance from CF1
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([CFlat[0],latex[ii]],[CFlon[0],lonex[ii]])[0][0]

    return distout

newdist=-1*getdist(newlons,newlats)

figure()
plot(bathdist,-bathbath,label='Knudsen depths')
plot(distvec,CFdpth,'o',label='CF depths from etopo')
plot(newdist,newdpth,'o-',label='Etopo extension')
axhline(0,color='k')
axvline(-15,color='k')
axvline(-10,color='k')
xlim([-30,40])
ylim([-250,20])
legend()
savefig('../figures/compare_shipmoor/Coastal_shelfextend.pdf',bbox_inches='tight')
