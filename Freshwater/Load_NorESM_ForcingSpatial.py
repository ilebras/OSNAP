from firstfuncs_1618 import *
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'
dat1=xr.open_dataset(glob.glob(datadir+'NorESM/*sources_spatial*0912.nc')[0])
dat2=xr.open_dataset(glob.glob(datadir+'NorESM/*sources_spatial*1812.nc')[0])

dat=xr.concat([dat1,dat2],dim='time',data_vars='minimal').sel(x=slice(50,150)).sel(y=slice(300,370))
for kk in dat:
    if 'NORDIC_' in kk:
        dat=dat.rename({kk:kk[7:]})
startyear=2000
startmonth=1
endyear=2018
endmonth=12
dat['time']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])

dat['pe']=(-dat['evap']+dat['solprec']+dat['liqprec'])*60**2*24
dat['si']=(dat['mltfrz'])*60**2*24
dat['ro']=dat['runoff']*60**2*24

import cartopy.crs as ccrs
import cartopy

central_lon, central_lat = -20, 72
extent = [-45, 5, 55, 84]


############# Map P-E and sea ice meltfreeze next to each other and count up amount in polar and deep limbs
osnap=xr.open_dataset(datadir+'NorESM/NorESM_osnap_xray_1912.nc')
ns=xr.open_dataset(datadir+'NorESM/NorESM_ns_xray_1912.nc')
fs=xr.open_dataset(datadir+'NorESM/NorESM_fs_xray_1912.nc')
bso=xr.open_dataset(datadir+'NorESM/NorESM_bso_xray_1912.nc')

osnap.LONGITUDE.values[:5][::-1]

wbnd_lon=hstack((fs.LONGITUDE.values[:10],-20,-30,osnap.LONGITUDE.values[:5][::-1],-45,fs.LONGITUDE.values[0]))
wbnd_lat=hstack((fs.LATITUDE.values[:10],70,65,osnap.LATITUDE.values[:5][::-1],76,fs.LATITUDE.values[0]))

ebnd_lon=hstack((fs.LONGITUDE.values[9:],bso.LONGITUDE.values,20,ns.LONGITUDE.values[::-1],osnap.LONGITUDE.values[4:][::-1],-30,-20,fs.LONGITUDE.values[9]))
ebnd_lat=hstack((fs.LATITUDE.values[9:],bso.LATITUDE.values,62,ns.LATITUDE.values[::-1],osnap.LATITUDE.values[4:][::-1],65,70,fs.LATITUDE.values[9]))
plot(ebnd_lon,ebnd_lat)

import regionmask

wbnd_tpl=array(tuple((zip(wbnd_lon,wbnd_lat))))
ebnd_tpl=array(tuple((zip(ebnd_lon,ebnd_lat))))

wbox = regionmask.Regions([wbnd_tpl], name='west box')
ebox = regionmask.Regions([ebnd_tpl], name='east box')

#create masks
emask = ebox.mask(dat.plon.values,dat.plat.values)
wmask = wbox.mask(dat.plon.values,dat.plat.values)

dat=dat.rename({'y':'lon_idx','x':'lat_idx'})

dat['lon_idx']=range(len(dat['lon_idx']))
dat['lat_idx']=range(len(dat['lat_idx']))

wmask.plot()

emask.plot()

east_dat=dat.where(emask==0)
west_dat=dat.where(wmask==0)

noresm=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_18yrs_2004.nc')
noresm['pe']=(noresm['solprec']+noresm['liqprec']+noresm['evap'])


dat['pe_tot']=(dat['pe'].where(dat.OSNAP_mask==1)*dat['parea']).sum(dim='lon_idx').sum(dim='lat_idx')/1e3/1e6/(60**2*24)

for xray in [east_dat,west_dat,dat]:
    for var in ['pe','si','ro']:
        xray[var+'_tot']=(xray[var]*xray['parea']).sum(dim='lon_idx').sum(dim='lat_idx')/1e3/1e6/(60**2*24)

"""Get the etopo1 data"""
from map_funcs import *
etopo1name=predir+'ETOPO1_Ice_g_gmt4.grd'
etopo1 = Dataset(etopo1name,'r')

lons = etopo1.variables["x"][:]
lats = etopo1.variables["y"][:]
res = findSubsetIndices(extent[2]-20,extent[3]+20,extent[0]-80,extent[1]+80,lats,lons)
lon,lat=np.meshgrid(lons[int(res[0]):int(res[1])],lats[int(res[2]):int(res[3])])
bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
bathySmoothed = laplace_filter(bathy,M=None)

def map_FW():
    fig, axx = plt.subplots(1,3, figsize=(9,5*9/13), subplot_kw={'projection': ccrs.Orthographic(central_lon, central_lat)},constrained_layout=True)
    data_crs = ccrs.PlateCarree()
    clvec=[3,3,20]
    labs=['Precipitation - Evaporation','Runoff','Sea ice melt and freeze',]
    for ax in axx:
            ax.set_extent(extent)
            ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor=None,color='k')
            bb=ax.contour(lon,lat,bathySmoothed,levels=[-400],colors='grey',transform=data_crs,zorder=2)
            gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1, color='gray', alpha=0.5,linestyle='--')
            for sec in [osnap,ns,fs,bso]:
                ax.plot(wbnd_lon,wbnd_lat,'k-',zorder=100,transform=data_crs,linewidth=3)
                ax.plot(ebnd_lon,ebnd_lat,'k-',zorder=100,transform=data_crs,linewidth=3)

    for ii,var in enumerate(['pe','ro','si']):
        hh=axx[ii].pcolor(dat.plon,dat.plat,dat[var].mean(dim='time'),transform=data_crs,cmap=cm.BrBG,vmin=-clvec[ii],vmax=clvec[ii])#levels=41,,extend='both'
        ypos=0
        ywi=0.04
        if ii==0:
            colorbar(hh,ax=axx[ii],cax=fig.add_axes([0.075,ypos,0.5,ywi]),orientation='horizontal',
            ticks=arange(-clvec[ii],clvec[ii]+0.2,1),
            label='[mm day$^{-1}$ m$^{-2}$]',ticklocation='bottom')
        elif ii==2:
            colorbar(hh,ax=axx[ii],cax=fig.add_axes([0.72,ypos,0.25,ywi]),orientation='horizontal',
            ticks=arange(-clvec[ii],clvec[ii]+1,10),
            label='[mm day$^{-1}$ m$^{-2}$]',ticklocation='bottom')
        axx[ii].set_title(labs[ii],fontsize=14)
        axx[ii].text(-46,72,str(int(mean(west_dat[var+'_tot'])*1e3+0.5))+' mSv',transform=data_crs,fontsize=14,color='white',zorder=300)
        axx[ii].text(-7,72,str(int(mean(east_dat[var+'_tot'])*1e3+0.5))+' mSv',transform=data_crs,fontsize=14,color='black',zorder=300)
    savefig(figdir_paper+'EastWest_Fresh_Map.png',bbox_inches='tight')
    savefig(figdir_paper+'EastWest_Fresh_Map.pdf',bbox_inches='tight')

map_FW()

########################### save


# def map_FW_2x2test():
#     fig, axx = plt.subplots(2,2, figsize=(11,12), subplot_kw={'projection': ccrs.Orthographic(central_lon, central_lat)},constrained_layout=True)
#     data_crs = ccrs.PlateCarree()
#     clvec=[3,20,3]
#     cdiv=[3,4,3]
#     ypos=0.9
#     ywi=0.02
#     clabs=['Precipitation - Evaporation\n[mm day$^{-1}$ m$^{-2}$]','Sea ice melt and freeze\n[mm day$^{-1}$ m$^{-2}$]','Runoff\n[mm day$^{-1}$ m$^{-2}$]']
#     caxvec=[fig.add_axes([0.15,ypos,0.3,ywi]),fig.add_axes([0.575,ypos,0.3,ywi]),fig.add_axes([0.15,ypos-0.4,0.3,ywi])]
#     for axi in axx:
#         for ax in axi:
#             ax.set_extent(extent)
#             ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor=None,color='k')
#             bb=ax.contour(lon,lat,bathySmoothed,levels=[-400],colors='grey',transform=data_crs,zorder=2)
#             for sec in [osnap,ns,fs,bso]:
#                 ax.plot(wbnd_lon,wbnd_lat,'k-',zorder=100,transform=data_crs,linewidth=3)
#                 ax.plot(ebnd_lon,ebnd_lat,'k-',zorder=100,transform=data_crs,linewidth=3)
#     for ii,var in enumerate(['pe','si','ro']):
#         if ii<2:
#             jj=0
#         else:
#             jj=1
#         hh=axx[jj,mod(ii,2)].pcolor(dat.plon,dat.plat,dat[var].mean(dim='time'),transform=data_crs,cmap=cm.BrBG,vmin=-clvec[ii],vmax=clvec[ii])#levels=41,,extend='both'
#         colorbar(hh,ax=axx[jj,mod(ii,2)],cax=caxvec[ii],orientation='horizontal',ticks=arange(-clvec[ii],clvec[ii]+clvec[ii]/6,clvec[ii]/cdiv[ii]),label=clabs[ii],ticklocation='top')
#     axx[1,1].cla()
#     savefig(figdir_paper+'EastWest_Fresh_Map.png',bbox_inches='tight')
#     savefig(figdir_paper+'EastWest_Fresh_Map.pdf',bbox_inches='tight')

# fig, axx = plt.subplots(1,2, figsize=(11,3),sharey=True, sharex=True)
# lili=2
# wmean=west_dat.groupby('time.month').mean(dim='time')
# wstd=west_dat.groupby('time.month').std(dim='time')
# emean=east_dat.groupby('time.month').mean(dim='time')
# estd=east_dat.groupby('time.month').std(dim='time')
# alf=1
# axx[0].fill_between((wmean.pe_tot-wstd.pe_tot).values,(wmean.pe_tot+wstd.pe_tot).values,color='grey',alpha=alf)
# axx[0].fill_between(emean.pe_tot-estd.pe_tot,emean.pe_tot+estd.pe_tot,color='grey',alpha=alf)
# wmean.pe_tot.plot(label='',ax=axx[0],color='k',linewidth=lili,linestyle='--')
# emean.pe_tot.plot(label='',ax=axx[0],color='k',linewidth=lili)
# axx[1].fill_between(wmean.si_tot-wstd.si_tot,wmean.si_tot+wstd.si_tot,color='grey',alpha=alf)
# axx[1].fill_between(emean.si_tot-estd.si_tot,emean.si_tot+estd.si_tot,color='grey',alpha=alf)
# wmean.si_tot.plot(label='Western region',ax=axx[1],color='k',linewidth=lili,linestyle='--')
# emean.si_tot.plot(label='Eastern region',ax=axx[1],color='k',linewidth=lili)
# for axi in axx:
#     axi.set_xlabel('')
#     axi.set_ylabel('')
# axx[0].set_ylabel('Fresh water transport [Sv]')
# axx[0].set_xlim(1,12)
# axx[1].legend(loc=(-1.1,0.7))
# savefig(figdir_paper+'EastWest_Fresh_2.png',bbox_inches='tight')
# savefig(figdir_paper+'EastWest_Fresh_2.pdf',bbox_inches='tight')
