from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM/'

################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################

#### load noresm and force formatting to be the same
def load_noresm(whichone):
    noresm=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*new.nc')[0])
    print(noresm.time)
    noresm=noresm.rename({'time': 'TIME','depth':'DEPTH','cell':'LONGITUDE'})
    noresm=noresm.rename({'salt':'PSAL','temp':'PTMP'})
    noresm=noresm.assign_coords(LONGITUDE=(noresm.lon.values))
    noresm=noresm.assign_coords(LATITUDE=(noresm.lat.values))
    noresm=noresm.drop('lat').drop('lon')
    startyear=2010
    startmonth=1
    endyear=2018
    endmonth=12
    noresm['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    # noresm=noresm.sel(TIME=slice(osnap.TIME[0],osnap.TIME[-1]))
    return noresm

osnap=load_noresm('OSNAP')

osnap
fs=load_noresm('FS')
bso_redundant=load_noresm('BSO')
ns=load_noresm('NS')

ns.LATITUDE.values

ns.LONGITUDE.values

def fix_bso(bso1): #
    lonind=[0,2,4,5,7,9,10,12,14,16,17,19,21,22,23]
    bso2=xr.Dataset({'PTMP':(['TIME','DEPTH','LONGITUDE'],bso1.PTMP.values[:,:,lonind]),
                     'PSAL':(['TIME','DEPTH','LONGITUDE'],bso1.PSAL.values[:,:,lonind]),
                     'uvel':(['TIME','DEPTH','LONGITUDE'],bso1.uvel.values[:,:,lonind]),
                     'vvel':(['TIME','DEPTH','LONGITUDE'],bso1.vvel.values[:,:,lonind]),
                     'vertices_lat':(['corner','LONGITUDE'],bso1.vertices_lat[:,lonind]),
                     'vertices_lon':(['corner','LONGITUDE'],bso1.vertices_lon[:,lonind]),
                     'phi':(['LONGITUDE'],bso1.phi[lonind]),
                     'cell_depth':(['LONGITUDE'],bso1.cell_depth[lonind]),
                     'dx':(['LONGITUDE'],bso1.dx[lonind]),
                     'dy':(['LONGITUDE'],bso1.dy[lonind]),},
                    coords={'LONGITUDE':bso1.LONGITUDE.values[lonind],'LATITUDE':bso1.LATITUDE.values[lonind],'TIME':bso1.TIME.values,'DEPTH':bso1.DEPTH.values})
    return bso2


bso=fix_bso(bso_redundant)

def add_PDEN(xray):
    if 'PRES' in list(xray.data_vars):
            PRES_out=xray['PRES']
    else:
            PRES_out=gsw.p_from_z(-xray['DEPTH'],60)
    SA_out=NaN*xray['PSAL'].copy()
    for ii,pp in enumerate(PRES_out):
        for jj,ll in enumerate(xray.LONGITUDE):
            SA_out[:,ii,jj]=gsw.SA_from_SP(xray['PSAL'][:,ii,jj],pp,ll,xray.LATITUDE[jj])
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['PDEN']=(('TIME','DEPTH','LONGITUDE'),PD_out)

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
    if xx=='osnap':
        ax1.plot(xrw.LONGITUDE[xrw.VFILTER.values==1],xrw.LATITUDE[xrw.VFILTER.values==1],'o',markersize=30,zorder=2)
        ax1.plot(xrw.LONGITUDE[xrw.VFILTER_PLUS.values==1],xrw.LATITUDE[xrw.VFILTER_PLUS.values==1],'o',markersize=20,zorder=3)
    ax1.plot(xrw.vertices_lon.values[xrw.vertices_lon.values<180],xrw.vertices_lat.values[xrw.vertices_lon.values<180],'ko')
    ax1.plot(xrw.vertices_lon.values[xrw.vertices_lon.values>180]-360,xrw.vertices_lat.values[xrw.vertices_lon.values>180],'ko')
    ax1.quiver(xrw.LONGITUDE.values,xrw.LATITUDE.values,xrw.uvel.sum(dim='DEPTH').mean(dim='TIME').values,xrw.vvel.sum(dim='DEPTH').mean(dim='TIME').values,zorder=4)
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


salvec=linspace(31,36,103)
tmpvec=linspace(-3,16,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])

SA_vec_1000=gsw.SA_from_SP(salvec,1e3*ones(len(salvec)),CFlon[3],CFlat[4])

CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
sigma1mat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3
        sigma1mat[jj,ii]=gsw.sigma1(SA_vec[ii],CT_vec[jj])


# def TS_comp_plot():
#     figure()
#     dat=fs
#     plot(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),'o',label='FS',alpha=0.2)
#     dat=bso
#     plot(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),'o',label='BSO',alpha=0.2)
#     dat=osnap
#     plot(dat.PSAL.values.flatten(),dat.PTMP.values.flatten(),'o',label='OSNAP',alpha=0.2)
#     xlabel('salinity')
#     ylabel('pot.temperature')
#     xlim(31,36)
#     ylim(-3,14)
#     contour(salvec,tmpvec,pdenmat,colors='k',levels=arange(25,29,0.2),zorder=5)
#     legend(fontsize=16)
#     savefig(figdir+'TS_NorESM_comp_all.png',bbox_inches='tight')
#
#
# TS_comp_plot()

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


osnap['DZ']=(['DEPTH','LONGITUDE'],tile(osnap_dist/osnap_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(osnap_dist),1]).T)
fs['DZ']=(['DEPTH','LONGITUDE'],tile(fs_dist/fs_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(fs_dist),1]).T)
ns['DZ']=(['DEPTH','LONGITUDE'],tile(ns_dist/ns_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(ns_dist),1]).T)
bso['DZ']=(['DEPTH','LONGITUDE'],tile(bso_dist/bso_dist,[len(nor_depthdiff),1])*tile(nor_depthdiff,[len(bso_dist),1]).T)

#correct DZ based on cell_depth
def corr_DZ(xrw):
    for ii,ff in enumerate(xrw.cell_depth.values):
        dgind=where(nor_depthgrid>ff)[0][0]
        xrw['DZ'][dgind-1,ii]=ff-nor_depthgrid[dgind-1]
        xrw['DZ'][dgind:,ii]=NaN
    return xrw

for xrw in [ns,fs,bso,osnap]:
    xrw=corr_DZ(xrw)

osnap['DIFFDIST']=(['LONGITUDE'],osnap_dist)
fs['DIFFDIST']=(['LONGITUDE'],fs_dist)
ns['DIFFDIST']=(['LONGITUDE'],ns_dist)
bso['DIFFDIST']=(['LONGITUDE'],bso_dist)


def get_TRANS(xrw):
    xrw['AREA']=xrw['DZ']*xrw['DIFFDIST']
    xrw['TRANS']=xrw['AREA']*xrw['VELO']/1e6
    return xrw

for ff in [osnap,fs,bso,ns]:
    ff=get_TRANS(ff)

osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
(osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot()

bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE').plot()
(bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot()


(osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot()
(bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot()


(osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).plot()
(osnap['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')+ns['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-bso['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')-fs['TRANS'].sum(dim='DEPTH').sum(dim='LONGITUDE')).mean()


osnap.to_netcdf(datadir+'NorESM/NorESM_osnap_xray_1912.nc','w')
fs.to_netcdf(datadir+'NorESM/NorESM_fs_xray_1912.nc','w')
ns.to_netcdf(datadir+'NorESM/NorESM_ns_xray_1912.nc','w')
bso.to_netcdf(datadir+'NorESM/NorESM_bso_xray_1912.nc','w')
