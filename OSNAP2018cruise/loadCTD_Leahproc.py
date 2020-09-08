from AR30_funcs import *

from seabird import *

bin_list=sort(glob.glob(datadir+'CTD/cruiseproc_Leah/*'))

profile = fCNV(bin_list[0])
profile.attributes
profile.attributes['datetime']
profile.keys()

# obase=profile['oxygen_ml_L'].data
#
# nanind=obase>1
#
# prsvec=range(2,2*prslen+1,2)
# array(prsvec)[nanind]
# obase[nanind]

prslen=1614
stanum=len(bin_list)
def loadCTD_bin(datalist):
    lat=zeros(len(datalist))
    lon=zeros(len(datalist))
    date=[1]*len(datalist)
    prs=NaN*ones((prslen,len(datalist)))
    sal=NaN*ones((prslen,len(datalist)))
    SA=NaN*ones((prslen,len(datalist)))
    tmp=NaN*ones((prslen,len(datalist)))
    CT=NaN*ones((prslen,len(datalist)))
    den=NaN*ones((prslen,len(datalist)))
    o2=NaN*ones((prslen,len(datalist)))
    prsvec=array(range(2,2*prslen+1,2))
    for ii,dd in enumerate(datalist):
        dat=fCNV(dd)
        lat[ii]=dat.attributes['LATITUDE']
        lon[ii]=dat.attributes['LONGITUDE']
        date[ii]=dat.attributes['datetime']
        pl=len(dat['PRES'].data)
        pst=int((dat['PRES'].data[0]-2)/2)
        prs[pst:pl+pst,ii]=dat['PRES'].data
        sal[pst:pl+pst,ii]=dat['PSAL'].data
        tmp[pst:pl+pst,ii]=dat['TEMP'].data
        SA[pst:pl+pst,ii]=gsw.SA_from_SP(sal[:pl,ii],prs[:pl,ii],lon[ii],lat[ii])
        CT[pst:pl+pst,ii]=gsw.CT_from_pt(SA[:pl,ii],tmp[:pl,ii])
        den[pst:pl+pst,ii]=gsw.sigma0(SA[:pl,ii],CT[:pl,ii])
        obase=dat['oxygen_ml_L'].data
        nanind=obase>1
        prel=prsvec[pst:pl+pst]
        ofunc=interpolate.interp1d(prel[nanind],obase[nanind],kind='linear',fill_value='NaN',bounds_error=False)
        o2[:,ii]=ofunc(prsvec)

    distvec=array([NaN]*stanum)
    for ss in seclab:
        if len(seclab[ss])>1:
            distvec[seclab[ss]]=getdist(lon[startsec[ss]],lat[startsec[ss]],lon[seclab[ss]],lat[seclab[ss]])
        else:
            distvec[seclab[ss]]=zeros(1)


    grdat=xr.Dataset({'sal': (['prs','sta'],  sal),
                      'tmp': (['prs','sta'],  tmp),
                      'SA': (['prs','sta'],  SA),
                      'CT': (['prs','sta'],  CT),
                      'den': (['prs','sta'],  den),
                      'o2': (['prs','sta'],  o2),
                      'lat': (['sta'],  lat),
                      'lon': (['sta'],  lon),
                      'dist': (['sta'], distvec)},
                    coords={'sta': range(stanum),
                            'prs': prsvec})

    return grdat

grdat=loadCTD_bin(bin_list)

for ss in ['section 2']:
    plot(grdat.lon[seclab[ss]],grdat.lat[seclab[ss]],'o',label=ss)
    plot(grdat.lon[startsec[ss]],grdat.lat[startsec[ss]],'kx',label='')
    legend(loc=(1.05,0))
    figure()
    plot(grdat.sta[seclab[ss]],grdat.dist[seclab[ss]],'.')

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



uni['o2']={}
uni['o2']['vmin']=5
uni['o2']['vmax']=9
uni['o2']['cmap']=cm.rainbow

for ss in seclab:
    f,axx=subplots(1,3,sharex=True,sharey=True,figsize=(20,4))
    pcolorsec('sal',axx[0],ss)
    pcolorsec('tmp',axx[1],ss)
    # f,axx=subplots(1,1)
    pcolorsec('o2',axx[2],ss)
    axx[0].set_ylabel('pressure [db]')
    axx[1].set_xlabel('distance [km]')
    suptitle('Section '+ss)
    gca().invert_yaxis()

#save 2m binned CTD into a netcdf to load elsewhere

grdat.to_netcdf('OSNAP2018cruise/data/CTD_2m_Leahproc.nc','w')

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
