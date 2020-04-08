from map_funcs import *

########################################################################################
################################## 2014 KN221 section ########################################
########################################################################################

prsvec=arange(2,3180,2) #determined in retrospect, found out what maximum pressure was

def loadCTD_cnv(datalist,key):
    lat[key]=zeros(len(datalist))
    lon[key]=zeros(len(datalist))
    date[key]=[1]*len(datalist)
    sal[key]=zeros((len(prsvec),len(datalist)))
    tmp[key]=zeros((len(prsvec),len(datalist)))
    for ii,dd in enumerate(datalist):
        profile = cnv.fCNV(dd)
        lat[key][ii]=profile.attributes['LATITUDE']
        lon[key][ii]=profile.attributes['LONGITUDE']
        date[key][ii]=profile.attributes['datetime']
        salfnc=interpolate.interp1d(profile['PRES'],profile['PSAL'],kind='linear',fill_value='NaN',bounds_error=False)
        tmpfnc=interpolate.interp1d(profile['PRES'],profile['TEMP'],kind='linear',fill_value='NaN',bounds_error=False)
        sal[key][:,ii]=salfnc(prsvec)
        tmp[key][:,ii]=tmpfnc(prsvec)

    return lat,lon,sal,tmp,date


lon={}
lat={}
sal={}
tmp={}
date={}

# ctdlist_2014=sort(glob.glob(datadir+'Shipboard/kn221_2014/cnv_files/*.cnv'))[1:]
# lat,lon,sal,tmp,date=loadCTD_cnv(ctdlist_2014,'2014_full')

folders_2014=sort(glob.glob(datadir+'Shipboard/kn221_2014/cnv_files/k*/'))
for ff,folder in enumerate(folders_2014):
    ctdlist=sort(glob.glob(folder+'*.cnv'))
    lat,lon,sal,tmp,date=loadCTD_cnv(ctdlist,'2014_'+str(ff+1))

for key in lat:
    if '2014_' in key:
        plot(lon[key],lat[key],'o',label=key)
legend(loc=(1.05,0.1))



ctdlist_2016=sort(glob.glob(datadir+'Shipboard/ar07_2016/cnv_files/*.cnv'))[1:]
lat,lon,sal,tmp,date=loadCTD_cnv(ctdlist_2016,'2016')

secs16={}
secs16['1']=where(lat['2016']>61.0)[0]
secs16['2']=where((lat['2016']>60.75) & (lat['2016']<61.0))[0]
secs16['3']=where((lon['2016']>-44) & (lat['2016']>60.4) & (lat['2016']<60.6))[0]
secs16['4']=where((lon['2016']>-44) & (lat['2016']>60.2) & (lat['2016']<60.4))[0]
secs16['5']=range(36,47) #determined by inspection
##note: skipping mooring casts between these two
secs16['6']=range(33) #determined by inspection
secs16['7']=range(92,115) #determined by inspection
secs16['8']=range(47,68) #determined by inspection
secs16['9']=range(222,240)
secs16['10']=range(295,307)
secs16['11']=range(307,315)
secs16['12']=range(282,295)
secs16['13']=range(141,164)
secs16['14']=range(68,92) #determined by inspection
secs16['14b']=range(240,268)
secs16['15']=where(lon['2016']<-48)[0]

secs16

###########################################################################
## SOME MAPS
###########################################################################
# map=makeMap('fullsouth')
# map.plot(lon['2014_full'],lat['2014_full'],'ko',markersize=12,latlon=True)
# title('CTD stations on 2014 OSNAP cruise')
# for key in lat:
#     if '2014' in key:
#         if 'full' not in key:
#             map.plot(lon[key],lat[key],'o',label=key[-1],latlon=True)
# map.plot(CFlon,CFlat,'w*',label='CF moorings',latlon=True,markersize=11)
# legg=legend(loc=(1.02,0.2))
# legg.get_frame().set_facecolor('lightgrey')
# savefig('../figures/Map/CTDmap_2014.pdf',bbox_inches='tight')
#
#
# map=makeMap('fullsouth')
# map.plot(lon['2016'],lat['2016'],'ko',markersize=12,latlon=True)
# title('CTD stations on 2016 OSNAP cruise')
# for key in secs16:
#     if 'b' in key:
#         map.plot(lon['2016'][secs16[key]],lat['2016'][secs16[key]],'wx',label=key,latlon=True)
#     else:
#         map.plot(lon['2016'][secs16[key]],lat['2016'][secs16[key]],'o',label=key,latlon=True)
# map.plot(CFlon,CFlat,'w*',label='CF moorings',latlon=True,markersize=11)
# legg=legend(loc=(1.02,0.2))
# legg.get_frame().set_facecolor('lightgrey')
# savefig('../figures/Map/CTDmap_2016.pdf',bbox_inches='tight')
#
# # I had been loading all points here, wanted to see them all together
# map=makeMap('fullsouth')
# map.plot(lon['2014_full'],lat['2014_full'],'ko',markersize=12,latlon=True,label='2014')
# map.plot(lon['2016'],lat['2016'],'ro',latlon=True,label='2016')
# map.plot(CFlon,CFlat,'w*',latlon=True,label='CF moorings',markersize=10);
# legend()
# savefig('../figures/Map/CTDmap_1416.pdf',bbox_inches='tight')
# savefig('../figures/Map/CTDmap_1416.png',bbox_inches='tight',dpi=500)



def subit(field,key):
    if shape(shape(field['2016']))[0]>1:
        field['L2016_'+key]=field['2016'][:,secs16[key]]
    else:
        field['L2016_'+key]=array(field['2016'])[secs16[key]]
    return field


for key in secs16:
    lat=subit(lat,key)
    lon=subit(lon,key)
    sal=subit(sal,key)
    tmp=subit(tmp,key)
    date=subit(date,key)

io.savemat('../data/Shipboard/2016LON_forLeah.mat',lon)
io.savemat('../data/Shipboard/2016LAT_forLeah.mat',lat)

lon

plot(lon['2014_3'],lat['2014_3'],'o')
plot(lon['2016_6'],lat['2016_6'],'o')

def delit(myDict):
    if '2016' in myDict: del myDict['2016']

delit(lat)
delit(lon)
delit(sal)
delit(tmp)
delit(date)

# Get potential density
pden={}
for key in sal:
    print(key)
    pden[key]=ma.zeros(shape(sal[key]))
    for ll in range(len(lat[key])):
        salnan=~isnan(sal[key][:,ll])
        if sum(salnan)>2:
            SA=gsw.SA_from_SP(sal[key][salnan,ll],prsvec[salnan],[CFlon[3]]*sum(salnan),[CFlat[3]]*sum(salnan))
            pden[key][salnan,ll]=gsw.sigma0(SA,gsw.CT_from_t(SA,tmp[key][salnan,ll],prsvec[salnan]))
    pden[key][pden[key]==0]=ma.masked

def plot_trackdir(yrstr):
    figure(figsize=(10,4))
    subplot(121)
    for key in lat:
        if yrstr in key:
            scatter(lon[key],lat[key],label=key,c=range(len(lat[key])),cmap=cm.hot_r)
    colorbar()
    subplot(122)
    for key in lat:
        if yrstr in key:
            plot(lon[key],lat[key],'o',label=key)
    gca().set_yticklabels('')
    suptitle(yrstr,fontsize=20)
    legend(loc=(1.05,0))



plot_trackdir('2014')
plot_trackdir('2016')

### to order these in a sensible manner, going to define a point in Greenland interior, get distances from first point on each track.
### then make distances negative if they are closer to the interior of Greenland than reference point. Finally, can just add the biggest negative value and distance will be referenced to innermost point.
### but first, going to focus on the tracks close to mine and just get diff from point that makes sense for comparison with moorings.



def getdist(lonex,latex):
    # get distance from longitude of CF1, but projected onto the nearest point in latitude for line
    # this dist is only really useful for the lines close to CF moorings at this point
    distout=zeros(len(lonex))
    lat0=latex[argmin(lonex-CFlon[0])]
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([lat0,latex[ii]],[CFlon[0],lonex[ii]])[0][0]
        if lonex[ii]<CFlon[0]:
            distout[ii]=-1*distout[ii]

    return distout

dist={}
for key in sal:
    dist[key]=getdist(lon[key],lat[key])


## need to order all the dictionaries by distance so that the horizontal gradients make some sort of sense!
for key in dist:
        dsortind=argsort(dist[key])
        lon[key]=lon[key][dsortind]
        lat[key]=lat[key][dsortind]
        date[key]=array(date[key])[dsortind]
        sal[key]=sal[key][:,dsortind]
        tmp[key]=tmp[key][:,dsortind]
        pden[key]=pden[key][:,dsortind]
        dist[key]=sort(dist[key])

dist.keys()



fcor=gsw.f(60)

geovel={}
geoshear={}
for key in sal:
    geovel[key]=ma.zeros(shape(sal[key]))
    geoshear[key]=(-diff(pden[key],axis=1)*9.8/fcor/1028/diff(dist[key])/1e3)
    geovel[key]=nancumsum(geoshear[key][::-1,:],axis=0)[::-1,:]*2


for key in geovel:
    figure()
    plot(pden[key],prsvec)
    gca().invert_yaxis()
    title(key)

for key in geovel:
    figure()
    plot(geovel[key],prsvec)
    gca().invert_yaxis()
    title(key)

#get in between distances, where the geostrophic velocity actually is:
middist={}
for key in dist:
    middist[key]=dist[key][:-1]+diff(dist[key])/2
    figure()
    plot(middist[key])

pickle.dump([prsvec,dist,middist,lon,lat,sal,tmp,pden,geoshear,geovel],open('../pickles/Shipboard/CTD_1416geo.pickle','wb'))
