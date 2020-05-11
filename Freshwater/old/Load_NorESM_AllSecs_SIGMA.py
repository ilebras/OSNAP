from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'

################################################################################################################################
#sigma2 levels copied over from rho.txt:
sig2lev=array([27.2200000000000, 27.7200000000000,28.2020000000000,28.6810000000000,29.1579999999999,29.6320000000001,
30.1020000000001,30.5670000000000,31.0260000000001,31.4770000000001,31.9200000000001,
32.3520000000001,32.7719999999999,33.1759999999999,33.5640000000001,33.9320000000000,34.2790000000000,
34.6020000000001,34.9000000000001,35.1720000000000,35.4169999999999,35.6369999999999,35.8320000000001,
36.0029999999999,36.1530000000000,36.2840000000001,36.3979999999999,36.4970000000001,
36.5840000000001,36.6600000000001,36.7280000000001,36.7890000000000,36.8430000000001,
36.8930000000000,36.9390000000001,36.9820000000000,37.0219999999999,37.0599999999999,
37.0960000000000,37.1310000000001,37.1659999999999,37.1990000000001,37.2310000000000,
37.2639999999999,37.2950000000001,37.3270000000000,37.3580000000000,37.3879999999999,
37.4190000000001,37.4500000000000,37.4800000000000,37.5799999999999,37.8000000000000])

plot(sig2lev,'.')

################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################
# figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM_testing/'
#
# for whichone in ['FS','NS','OSNAP','BSO']:
#     noresm1=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*_sigma_2000*.nc')[0])
#     noresm2=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*_sigma_2010*.nc')[0])
#     for var in ['salt','temp','uvel','vvel']:
#         figure()
#         noresm1[var].isel(time=10).plot()
#         title(whichone+', 2000-')
#         savefig(figdir+whichone+'_'+var+'_2000.png',bbox_inches='tight')
#         figure()
#         noresm2[var].isel(time=10).plot()
#         title(whichone+', 2010-')
#         savefig(figdir+whichone+'_'+var+'_2010.png',bbox_inches='tight')

#### load noresm and force formatting to be the same
def load_noresm(whichone):
    noresm2=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*_sigma_2010*.nc')[0])
    noresm1=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*_sigma_2000*.nc')[0])
    noresm=xr.concat([noresm1,noresm2],dim='time',data_vars='minimal')
    noresm=noresm.rename({'time': 'TIME','dz':'DZ','cell':'LONGITUDE'})
    noresm=noresm.rename({'salt':'PSAL','temp':'PTMP'})
    noresm=noresm.assign_coords(LONGITUDE=(noresm1.lon.values))
    noresm=noresm.assign_coords(LATITUDE=(noresm1.lat.values))
    noresm=noresm.drop('lat').drop('lon')
    noresm['DEPTH']=noresm['deltaz'].cumsum(dim='DZ')
    startyear=2000
    startmonth=1
    endyear=2018
    endmonth=12
    noresm['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return noresm

osnap=load_noresm('OSNAP')

fs=load_noresm('FS')


bso_redundant=load_noresm('BSO')

ns=load_noresm('NS')

def fix_bso(bso1): #
    lonind=[0,2,4,5,7,9,10,12,14,16,17,19,21,22,23]
    bso2=xr.Dataset({'PTMP':(['TIME','DZ','LONGITUDE'],bso1.PTMP.values[:,:,lonind]),
                     'PSAL':(['TIME','DZ','LONGITUDE'],bso1.PSAL.values[:,:,lonind]),
                     'uvel':(['TIME','DZ','LONGITUDE'],bso1.uvel.values[:,:,lonind]),
                     'vvel':(['TIME','DZ','LONGITUDE'],bso1.vvel.values[:,:,lonind]),
                     'vertices_lat':(['corner','LONGITUDE'],bso1.vertices_lat[:,lonind]),
                     'vertices_lon':(['corner','LONGITUDE'],bso1.vertices_lon[:,lonind]),
                     'deltaz':(['TIME','DZ','LONGITUDE'],bso1.deltaz.values[:,:,lonind]),
                     'DEPTH':(['TIME','DZ','LONGITUDE'],bso1.DEPTH.values[:,:,lonind]),
                     'phi':(['LONGITUDE'],bso1.phi[lonind]),
                     'cell_depth':(['LONGITUDE'],bso1.cell_depth[lonind]),
                     'dx':(['LONGITUDE'],bso1.dx[lonind]),
                     'dy':(['LONGITUDE'],bso1.dy[lonind]),},
                    coords={'LONGITUDE':bso1.LONGITUDE.values[lonind],'LATITUDE':bso1.LATITUDE.values[lonind],'TIME':bso1.TIME.values,'DZ':bso1.DZ.values})
    return bso2

bso=fix_bso(bso_redundant)


def add_PDEN(xray):
    PRES_out=gsw.p_from_z(-xray['DEPTH'],60)
    SA_out=NaN*xray['PSAL'].copy()
    for ii in range(shape(PRES_out)[1]):
        for jj,ll in enumerate(xray.LONGITUDE):
            SA_out[:,ii,jj]=gsw.SA_from_SP(xray['PSAL'][:,ii,jj],PRES_out[:,ii,jj],ll,xray.LATITUDE[jj])
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['PDEN']=(('TIME','DZ','LONGITUDE'),PD_out)

add_PDEN(osnap)
add_PDEN(fs)
add_PDEN(bso)
add_PDEN(ns)
################################################################################################################################
################################################################################################################################
###################################     GET A HANDLE ON VEL DIRECTIONS   #######################################################
################################################################################################################################
################################################################################################################################
# import cartopy.crs as ccrs
#
# def main():
#     fig = plt.figure(figsize=(8, 10))
#
#     x, y, u, v, vector_crs = sample_data(shape=(50, 50))
#     ax1 = fig.add_subplot(2, 1, 1, projection=ccrs.PlateCarree())
#     ax1.coastlines('50m')
#     ax1.set_extent([-45, 55, 20, 80], ccrs.PlateCarree())
#     ax1.quiver(x, y, u, v, transform=vector_crs)
#
#     ccrs.NorthPolarStereo

bso['VFILTER']=(['LONGITUDE'],array([1,1,1,0,1,1,0,1,1,0,1,1,0,1,0]))
osnap['VFILTER']=(['LONGITUDE'],hstack((array([0]*7),array([1]*10),array([0]*19))))
osnap['VFILTER_PLUS']=(['LONGITUDE'],hstack((array([0]*7),array([1,0,0,1,0,0,1,0,0,1]),array([0]*19))))

def plot_map_veldir(xrw,xx='na'):
    if xx=='bso':
        f,[ax1,ax2]=subplots(1,2,figsize=(12,8))
        ax1.plot(xrw.vertices_lon,xrw.vertices_lat,'ko-')
        ax1.plot(xrw.LONGITUDE[xrw.VFILTER.values==1],xrw.LATITUDE[xrw.VFILTER.values==1],'o',markersize=30,zorder=2)
        phi=ax1.scatter(xrw.LONGITUDE,xrw.LATITUDE,c=xrw.phi,zorder=3)
        colorbar(phi,ax=ax1)
        # ax1.quiver(bso.LONGITUDE.values,bso.LATITUDE.values,bso_ureal.sum(dim='DEPTH').mean(dim='TIME').values,bso_vreal.sum(dim='DEPTH').mean(dim='TIME').values,zorder=4,color='r')
    else:
        f,[ax1,ax2]=subplots(2,1,figsize=(12,8))
        ax1.plot(xrw.LONGITUDE,xrw.LATITUDE,'o-')
    # if xx=='osnap':
    #     ax1.plot(xrw.LONGITUDE[xrw.VFILTER.values==1],xrw.LATITUDE[xrw.VFILTER.values==1],'o',markersize=30,zorder=2)
    #     ax1.plot(xrw.LONGITUDE[xrw.VFILTER_PLUS.values==1],xrw.LATITUDE[xrw.VFILTER_PLUS.values==1],'o',markersize=20,zorder=3)
    ax1.plot(xrw.vertices_lon.values[xrw.vertices_lon.values<180],xrw.vertices_lat.values[xrw.vertices_lon.values<180],'ko')
    ax1.plot(xrw.vertices_lon.values[xrw.vertices_lon.values>180]-360,xrw.vertices_lat.values[xrw.vertices_lon.values>180],'ko')
    ax1.quiver(xrw.LONGITUDE.values,xrw.LATITUDE.values,xrw.uvel.sum(dim='DZ').mean(dim='TIME').values,xrw.vvel.sum(dim='DZ').mean(dim='TIME').values,zorder=4)
    ax2.plot(xrw.LONGITUDE,xrw.phi,'o-')
    ax2.axhline(0,color='k',alpha=0.2)
    ax2.set_ylim(-pi/2,pi/2,)


plot_map_veldir(osnap,'osnap')


plot_map_veldir(bso,'bso')



plot_map_veldir(fs)
plot_map_veldir(ns)


################################################################################################################################
################################################################################################################################
###########################################      PLOT TS SPACE    #######################################################
################################################################################################################################
################################################################################################################################

def TS_hex_plot(dat,titi):
    figure()
    hexbin(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),cmap='hot_r',mincnt=1,bins='log')
    colorbar(label='# of Observations')
    xlabel('salinity')
    ylabel('pot.temperature')
    title(titi)
    xlim(32,36)
    ylim(-2.5,14)
    savefig(figdir+'TShex_'+titi+'.png',bbox_inches='tight')


TS_hex_plot(osnap,'NorESM_OSNAP')
TS_hex_plot(fs,'NorESM_FS')
TS_hex_plot(bso,'NorESM_BSO')
TS_hex_plot(ns,'NorESM_NS')

################################################################################################################################
################################################################################################################################
#####################     Calc perpendicular velocity and Transport   #######################################################
################################################################################################################################
################################################################################################################################

fs['VELO']=fs['vvel'].copy()

ns['VELO']=ns['vvel'].copy()

################################################################################################################################
###########################################      Calculate BSO velocity carefully!    #######################################################
################################################################################################################################
# first, use phi to project all velocities to the true north-south velocity:
# BUT, only use "v" velocities that are not adjacent to other gridcells
# from morven:
# ur=us.*cos(phi)-vs.*sin(phi);
# vr=us.*sin(phi)+vs.*cos(phi);
bso_vvel_corr=bso.vvel.copy()
bso_vvel_corr[:,:,bso.VFILTER==0]=0
bso_ureal=bso.uvel*cos(bso.phi)-bso_vvel_corr*sin(bso.phi)
bso_vreal=bso.uvel*sin(bso.phi)+bso_vvel_corr*cos(bso.phi)
#now project onto direction perpendicular to track
lon1=bso.LONGITUDE[0]
lat1=bso.LATITUDE[0]
lon2=(bso.LONGITUDE[-2]+bso.LONGITUDE[-1])/2
lat2=(bso.LATITUDE[-2]+bso.LATITUDE[-1])/2
bso_angle=float(arctan((lat2-lat1)/(lon2-lon1)))

sin(bso_angle)
bso['VELO']=bso_vreal*abs(cos(bso_angle))+bso_ureal*abs(sin(bso_angle))


osnap['VELO']=osnap['vvel'].copy()
################################################################################################################################
#############################   Correct OSNAP velocities that are at an angle    #############################################
################################################################################################################################
olon1=osnap.LONGITUDE[where(osnap.VFILTER==1)[0][0]]
olon2=osnap.LONGITUDE[where(osnap.VFILTER==1)[0][-1]+1]
olat1=osnap.LATITUDE[where(osnap.VFILTER==1)[0][0]]
olat2=osnap.LATITUDE[where(osnap.VFILTER==1)[0][-1]+1]
osnap_angle=float(arctan((olat2-olat1)/(olon2-olon1)))
osnap_angle

#only use the u velocity of the corner gridcells
osnap_uvel_corr=osnap.uvel.copy()
osnap_uvel_corr[:,:,osnap.VFILTER_PLUS==0]=0
osnap_ureal=osnap_uvel_corr*cos(osnap.phi)-osnap.vvel*sin(osnap.phi)
osnap_vreal=osnap_uvel_corr*sin(osnap.phi)+osnap.vvel*cos(osnap.phi)

osnap['VELO'][:,:,osnap.VFILTER==1]=osnap_vreal[:,:,osnap.VFILTER==1]*abs(cos(osnap_angle))+osnap_ureal[:,:,osnap.VFILTER==1]*abs(sin(osnap_angle))

osnap_dist=array([57382.5094665826, 57367.6635038391, 57353.5409692236, 57340.1361679362, 57327.4434986111, 57315.4573352359, 57304.1721310823, 57293.5824492824,
57957.9063952388, 57948.4302559112, 57939.6484924917, 58619.6339361550, 58612.0555633080, 58605.1744353929, 59301.2369172294, 59295.6417676100, 59290.7490131978,
60003.1944856841, 59999.6566352716, 59996.8270689983, 59994.7052065002, 59993.2908518308, 59992.5837196858, 59992.5837196858, 59993.2908518308,59994.7052065002,
59996.8270689983, 59999.6566352716, 60003.1944856841, 60007.4409760742, 60012.3967978497, 60018.0626167572, 60024.4390999926, 60031.5270237768, 60039.3271940887, 60047.8404139799])
osnap_dist[osnap.VFILTER==1]=osnap_dist[osnap.VFILTER==1]*abs(cos(osnap_angle))

nor_depthgrid=array([0, 2.50000000000000, 7.50000000000000, 12.5000000000000, 17.5000000000000, 22.5000000000000, 27.5000000000000, 35, 45, 56.2000000000000, 68.8000000000000, 81.2000000000000, 93.8000000000000, 106.200000000000, 118.800000000000, 131.200000000000, 143.800000000000, 162.500000000000, 187.500000000000, 212.500000000000, 237.500000000000, 262.500000000000, 287.500000000000, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1562.50000000000, 1687.50000000000, 1812.50000000000, 1937.50000000000, 2125, 2375, 2625, 2875, 3125, 3375, 3625, 3875, 4125, 4375, 4625, 4875, 5125, 5375, 5625, 5875, 6125, 6375, 6625, 8000])
nor_depthdiff=diff(nor_depthgrid)

bso_dx_corr=bso.dx.copy()
bso_dx_corr[bso.VFILTER==0]=0
bso_dist=bso.dy*abs(cos(bso_angle))+bso_dx_corr*abs(sin(bso_angle))

fs_dist=fs.dx.values

ns_dist=ns.dx.values

osnap

# osnap['DZ']=(['DEPTH','LONGITUDE'],tile(osnap_dist/osnap_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(osnap_dist),1]).T)
# fs['DZ']=(['DEPTH','LONGITUDE'],tile(fs_dist/fs_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(fs_dist),1]).T)
# ns['DZ']=(['DEPTH','LONGITUDE'],tile(ns_dist/ns_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(ns_dist),1]).T)
# bso['DZ']=(['DEPTH','LONGITUDE'],tile(bso_dist/bso_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(bso_dist),1]).T)

# #correct DZ based on cell_depth
# def corr_DZ(xrw):
#     for ii,ff in enumerate(xrw.cell_depth.values):
#         dgind=where(nor_depthgrid>ff)[0][0]
#         xrw['DZ'][dgind-1,ii]=ff-nor_depthgrid[dgind-1]
#         xrw['DZ'][dgind:,ii]=NaN
#     return xrw
#
# for xrw in [ns,fs,bso,osnap]:
#     xrw=corr_DZ(xrw)

osnap['DIFFDIST']=(['LONGITUDE'],osnap_dist)
fs['DIFFDIST']=(['LONGITUDE'],fs_dist)
ns['DIFFDIST']=(['LONGITUDE'],ns_dist)
bso['DIFFDIST']=(['LONGITUDE'],bso_dist)


def get_TRANS(xrw):
    xrw['AREA']=xrw['deltaz']*xrw['DIFFDIST']
    xrw['TRANS']=xrw['AREA']*xrw['VELO']/1e6
    return xrw

for ff in [osnap,fs,bso,ns]:
    ff=get_TRANS(ff)

################################################################################################################################
###########################   correct transport based on model diagnostics   #######################################################
################################################################################################################################
startyear=2000
startmonth=1
endyear=2018
endmonth=12

def load_allyears(varlab):
    tmp1=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_200001-200912.nc')[0])
    tmp2=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_201001-201812.nc')[0])
    varout=xr.concat([tmp1,tmp2],dim='time')
    varout=varout.rename({'time':'TIME'})
    varout['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return varout

vts=load_allyears('volumetransports')
hts=load_allyears('heattransports')


def corr_trans(var,netvar):
    var.TRANS.sum(dim='DZ').sum(dim='LONGITUDE').plot(label='before')
    transdiff=(var.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')-vts[netvar])
    velcorr=transdiff/var.AREA.sum(dim='DZ').sum(dim='LONGITUDE')
    var['TRANS']=var['TRANS']-velcorr*var.AREA
    var.TRANS.sum(dim='DZ').sum(dim='LONGITUDE').plot(label='after')
    vts[netvar].plot(label='daignostic')
    legend(loc=(1.01,0.2))
    return var

osnap=corr_trans(osnap,'net_vt_OSNAP')
fs=corr_trans(fs,'net_vt_FS')
bso=corr_trans(bso,'net_vt_BSO')
ns=corr_trans(ns,'net_vt_NS')

cp=3850
rhow=1000

fs.PTMP.mean(dim='TIME').plot()

fs.TRANS.mean(dim='TIME').plot()


### maybe this could work if I used area instead of transport inthe correction...??
# def corr_temp(var,netvar):
#     (var.TRANS*var.PTMP).sum(dim='DZ').sum(dim='LONGITUDE').plot(label='before')
#     (hts[netvar]/cp/rhow/1e6).plot(label='diagnostic')
#     heatdiff=(var.TRANS*var.PTMP).sum(dim='DZ').sum(dim='LONGITUDE')-hts[netvar]/cp/rhow/1e6
#     tempcorr=heatdiff/var.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')
#     var['PTMP']=var['PTMP']-tempcorr
#     (var.TRANS*var.PTMP).sum(dim='DZ').sum(dim='LONGITUDE').plot(label='after')
#     legend(loc=(1.01,0.2))
#     return var
#
# osnap=corr_temp(osnap,'net_ht_OSNAP')


osnap.to_netcdf(datadir+'NorESM/NorESM_osnap_xray_18yrs_2004_sigma_transcorr.nc','w')
# fs.to_netcdf(datadir+'NorESM/NorESM_fs_xray_18yrs_2004_sigma_transcorr.nc','w')
# ns.to_netcdf(datadir+'NorESM/NorESM_ns_xray_18yrs_2004_sigma_transcorr.nc','w')
# bso.to_netcdf(datadir+'NorESM/NorESM_bso_xray_18yrs_2004_sigma_transcorr.nc','w')
