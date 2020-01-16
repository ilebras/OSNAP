from AR30_funcs import *

datadir

bin_list=sort(glob.glob(datadir+'CTD/binned/*'))

prslen=3164
prsvec=range(prslen)

stanum=len(bin_list)

dat=io.loadmat(bin_list[0])

def loadCTD_bin(datalist):
    # pmax=zeros(len(datalist))
    lat=zeros(len(datalist))
    lon=zeros(len(datalist))
    date=[1]*len(datalist)
    sal=zeros((prslen,len(datalist)))
    smsal=zeros((prslen,len(datalist)))
    SA=zeros((prslen,len(datalist)))
    tmp=zeros((prslen,len(datalist)))
    smtmp=zeros((prslen,len(datalist)))
    CT=zeros((prslen,len(datalist)))
    den=zeros((prslen,len(datalist)))
    o2=zeros((prslen,len(datalist)))
    smden=zeros((prslen,len(datalist)))
    for ii,dd in enumerate(bin_list):
        # print(ii)
        dat=io.loadmat(dd)
        # pmax[ii]=max(dat['datad_1m'][0][0][18]) #to determine max pressure...
        lat[ii]=nanmean(dat['datad_1m'][0][0][14])
        lon[ii]=nanmean(dat['datad_1m'][0][0][13])
        date[ii]=nanmean(dat['datad_1m'][0][0][15])
        prsprof=dat['datad_1m'][0][0][18].flatten()
        salprof=dat['datad_1m'][0][0][4].flatten()
        nanind=~isnan(salprof)
        if len(salprof[nanind])>1:
            salfnc=interpolate.interp1d(prsprof[nanind],salprof[nanind],kind='linear',fill_value='NaN',bounds_error=False)
            tmpfnc=interpolate.interp1d(prsprof[nanind],dat['datad_1m'][0][0][6][nanind].flatten(),kind='linear',fill_value='NaN',bounds_error=False)
            denfnc=interpolate.interp1d(prsprof[nanind],dat['datad_1m'][0][0][8][nanind].flatten(),kind='linear',fill_value='NaN',bounds_error=False)
            o2fnc=interpolate.interp1d(prsprof[nanind],dat['datad_1m'][0][0][10][nanind].flatten(),kind='linear',fill_value='NaN',bounds_error=False)
            sal[:,ii]=salfnc(prsvec)
            smit=5
            smsal[:,ii]=run_ave(sal[:,ii],smit)
            SA[:,ii]=gsw.SA_from_SP(sal[:,ii],prsvec,lon[ii],lat[ii])
            tmp[:,ii]=tmpfnc(prsvec)
            smtmp[:,ii]=run_ave(tmp[:,ii],smit)
            CT[:,ii]=gsw.CT_from_pt(SA[:,ii],tmp[:,ii])
            den[:,ii]=denfnc(prsvec)
            o2[:,ii]=o2fnc(prsvec)
            smden[:,ii]=run_ave(den[:,ii],smit)

    distvec=array([NaN]*stanum)
    for ss in seclab:
        if len(seclab[ss])>1:
            distvec[seclab[ss]]=getdist(lon[startsec[ss]],lat[startsec[ss]],lon[seclab[ss]],lat[seclab[ss]])
        else:
            distvec[seclab[ss]]=zeros(1)



    SAint=SA[:-1,:]+diff(SA,axis=0)
    CTint=CT[:-1,:]+diff(CT,axis=0)

    prsi=prsvec[:-1]+diff(prsvec)/2

    alphaint=gsw.alpha(SAint.T,CTint.T,prsi).T
    betaint=gsw.beta(SAint.T,CTint.T,prsi).T

    alphaT=-alphaint*diff(smtmp,axis=0)
    betaS=betaint*diff(smsal,axis=0)

    R=alphaT/betaS

    dendiff=diff(smden,axis=0)

    # turner_wden=arctan2((alphaT-betaS)*dendiff/abs(dendiff),(alphaT+betaS)*dendiff/abs(dendiff))*180/pi
    turner=arctan2((alphaT-betaS),(alphaT+betaS))*180/pi

    grdat=xr.Dataset({'sal': (['prs','sta'],  sal),
                      'tmp': (['prs','sta'],  tmp),
                      'den': (['prs','sta'],  den-1e3),
                      'o2': (['prs','sta'],  o2),
                      'turner': (['prsi','sta'],  turner),
                      # 'turner_wden': (['prsi','sta'],  turner_wden),
                      'dendiff': (['prsi','sta'], dendiff),
                      'lat': (['sta'],  lat),
                      'lon': (['sta'],  lon),
                      'dist': (['sta'], distvec)},
                    coords={'sta': range(stanum),
                            'prs': prsvec,
                            'prsi': prsi,})

    return grdat

grdat=loadCTD_bin(bin_list)

grdat

grdat.o2.max()

for ss in seclab:
    plot(grdat.lon[seclab[ss]],grdat.lat[seclab[ss]],'o',label=ss)
    plot(grdat.lon[startsec[ss]],grdat.lat[startsec[ss]],'kx',label='')
    legend(loc=(1.05,0))

grdat

# uni['dendiff']={}
# uni['dendiff']['cmap']=turn_map
# uni['dendiff']['vmin']=-0.0001
# uni['dendiff']['vmax']=0.0001

def pcolorsec(var,axi,ss):
    hh=axi.pcolor(sort(grdat.dist[seclab[ss]]),grdat.prs,grdat[var][:,seclab[ss]].sortby(grdat.dist),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
    if var=='turner':
        colorbar(hh,ax=axi,ticks=[-90,-45,0,45,90])
    else:
        colorbar(hh,ax=axi)
    axi.set_title(var)

for ss in seclab:
    f,axx=subplots(1,3,sharex=True,sharey=True,figsize=(20,4))
    pcolorsec('sal',axx[0],ss)
    pcolorsec('tmp',axx[1],ss)
    # f,axx=subplots(1,1)
    pcolorsec('turner',axx[2],ss)
    axx[0].set_ylabel('pressure [db]')
    axx[1].set_xlabel('distance [km]')
    suptitle('Section '+ss)
    gca().invert_yaxis()

#save 1m binned CTD into a netcdf to load elsewhere

grdat.to_netcdf('OSNAP2018cruise/data/CTD_1mbin.nc','w')

XXXXXXXXXXXX
# Keep this around for now
# ctd_list=sort(glob.glob(datadir+'proc/CTD/cruiseproc_Leah/*'))
#
# prsvec=range(2,3228,2) #determined from looping through below and finding max
# def loadCTD_cnv(datalist):
#     lat=zeros(len(datalist))
#     lon=zeros(len(datalist))
#     date=[1]*len(datalist)
#     sal=zeros((len(prsvec),len(datalist)))
#     tmp=zeros((len(prsvec),len(datalist)))
#     for ii,dd in enumerate(datalist):
#         profile = cnv.fCNV(dd);
#         lat[ii]=profile.attributes['LATITUDE']
#         lon[ii]=profile.attributes['LONGITUDE']
#         date[ii]=profile.attributes['datetime']
#         salfnc=interpolate.interp1d(profile['PRES'],profile['PSAL'],kind='linear',fill_value='NaN',bounds_error=False)
#         tmpfnc=interpolate.interp1d(profile['PRES'],profile['TEMP'],kind='linear',fill_value='NaN',bounds_error=False)
#         sal[:,ii]=salfnc(prsvec)
#         tmp[:,ii]=tmpfnc(prsvec)
#
#     grdat=xr.Dataset({'sal': (['prs','sta'],  sal),
#                       'tmp': (['prs','sta'],  tmp),
#                       'section': (['sta'],  seclab),
#                       'lat': (['sta'],  lat),
#                       'lon': (['sta'],  lon),},
#                     coords={'sta': range(len(datalist)),
#                             'prs': prsvec,})
#
#
#     return grdat
#
# dat=loadCTD_cnv(ctd_list)
#
# dat
