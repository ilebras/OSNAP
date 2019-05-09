from aux_funcs import *

nclist=glob.glob('/home/isabela/DATA/NARR/*.nc')


#Try extracting

#Use initial_time0_hours to get time from each, seems easiest to use. Build array that is lon X lat X time.
minlat=55
maxlat=70
minlon=-50
maxlon=-20

dat=Dataset(nclist[0])

dat

for na,dd in enumerate(sort(nclist)):
    dat=Dataset(dd)
    hourtmp=dat['initial_time0_hours'][:]
    date=[datetime.datetime(1800,1,1)+datetime.timedelta(hours=hh) for hh in hourtmp]


    tmp=xr.Dataset({'taux': (['date', 'y','x'],  dat['U_GRD_221_HTGL'][:]),
                    'tauy': (['date', 'y','x'],  dat['V_GRD_221_HTGL'][:]),},
                    coords={'lon':(['y','x'],dat['gridlon_221'][:]),
                            'lat':(['y','x'],dat['gridlat_221'][:]),
                            'date':date})

    tmp_locind=tmp.where(tmp.lon>=minlon).where(tmp.lon<=maxlon).where(tmp.lat>=minlat).where(tmp.lat<=maxlat).dropna(dim='x',how='all').dropna(dim='y',how='all')

    tmpsel=tmp_locind.resample('D',dim='date')

    if na==0:
        windxr=tmpsel
    else:
        windxr=xr.concat((windxr,tmpsel),dim='date')

dat
##Make map series of stress

import cartopy
import cartopy.crs as ccrs

skiparr=2

help(plt.axes)

def contmap(date1,date2,axi):
    cc=pcolormesh(windxr.lon,windxr.lat, sqrt(windxr['taux']**2+windxr['tauy']**2).sel(date=slice(date1,date2)).mean(dim='date'),transform=ccrs.PlateCarree(),cmap=cm.inferno,)
    cc=pcolormesh(windxr.lon,windxr.lat, (windxr['taux']).sel(date=slice(date1,date2)).mean(dim='date'),transform=ccrs.PlateCarree())
    # quiver(windxr.lon[::skiparr,::skiparr],windxr.lat[::skiparr,::skiparr],
            # windxr['taux'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],
            # windxr['tauy'].sel(date=slice(date1,date2)).mean(dim='date')[::skiparr,::skiparr],
            # color='white',transform=ccrs.PlateCarree(),scale=1.5)
    xlim([minlon+1,maxlon-1])
    ylim([minlat+1,maxlat-1])
    axi.coastlines(resolution='50m',zorder=10,color='w')
    return cc

# date1='2015-1-01'
# date2='2015-4-1'
#
# sqrt(windxr['taux']**2+windxr['tauy']**2).sel(date=slice(date1,date2)).mean(dim='date')

def plot3stress():
    xlp=-49.5
    tfs=18
    fig, (ax11, ax22, ax33) = plt.subplots(3,1, sharex=True, sharey=True,figsize=(6,12))
    ax11=subplot(311,projection=ccrs.PlateCarree())
    c2=contmap('2014-10-1','2015-1-1',ax11)
    text(xlp,67,'Fall',fontsize=tfs)
    cbaxes = fig.add_axes([1.0, 0.3, 0.06, 0.45])
    colorbar(c2,cax=cbaxes,label='Wind Stress magnitude [N m$^{-2}$]')#,ticks=arange(0,0.8,0.1))
    ax22=subplot(312,projection=ccrs.PlateCarree())
    c3=contmap('2015-1-01','2015-4-1',ax22)
    text(xlp,67,'Winter',fontsize=tfs)
    ax33=subplot(313,projection=ccrs.PlateCarree())
    c1=contmap('2015-04-01','2015-09-01',ax33)
    text(xlp,64.5,'Summer/ \nSpring',fontsize=tfs)
    suptitle('Its speed, not stress, and it looks weird!')
    # plt.tight_layout()
    # savefig('../figures/paperfigs/NARRstress_seas.pdf',bbox_inches='tight')

plot3stress()


def wint2():
    ax=subplot(111,projection=ccrs.PlateCarree())
    contmap('2015-12-01','2016-3-1',ax)

wint2()


def selloc(tmp,minlon,maxlon,minlat,maxlat):
    tmploc=tmp.where(tmp.lon>=minlon).where(tmp.lon<=maxlon).where(tmp.lat>=minlat).where(tmp.lat<=maxlat)
    return tmploc



CFwind=windxr
