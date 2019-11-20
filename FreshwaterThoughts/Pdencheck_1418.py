from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/'

#############################################
### Load data, establish some parameters
#############################################

dat_mat=io.loadmat(datadir+'OAproduct/0716_OSNAP_All_OA_201407-201807_ArgoMooringGlider.mat')
###[('time', 'O'), ('lon', 'O'), ('lat', 'O'), ('dist', 'O'), ('depth', 'O'), ('t', 'O'), ('s', 'O'), ('d', 'O'), ('pd', 'O'), ('pt', 'O'), ('depth_d', 'O'), ('t_err', 'O'), ('s_err', 'O'), ('t_clim', 'O'), ('s_clim', 'O'), ('nprof', 'O'), ('CREATION_DATE', 'O'), ('MATLAB_VERSION', 'O'), ('DATA_FILE', 'O')])}

mat_time=dat_mat['ts_analysis_mm'][0][0][0].flatten()
mat_lon=dat_mat['ts_analysis_mm'][0][0][1].flatten()
mat_lat=dat_mat['ts_analysis_mm'][0][0][2].flatten()
mat_dpth=dat_mat['ts_analysis_mm'][0][0][4].flatten()
mat_sal=dat_mat['ts_analysis_mm'][0][0][6]
mat_pden=dat_mat['ts_analysis_mm'][0][0][8]
mat_ptmp=dat_mat['ts_analysis_mm'][0][0][9]
mat_nprof=dat_mat['ts_analysis_mm'][0][0][-4]

python_time = array([datetime.datetime.fromordinal(int(mm)) + datetime.timedelta(days=mm%1) - datetime.timedelta(days = 366) for mm in mat_time])


dat=xr.Dataset({'PSAL':(['LONGITUDE','DEPTH','TIME'],mat_sal),
                'PTMP':(['LONGITUDE','DEPTH','TIME'],mat_ptmp),
                'PDEN':(['LONGITUDE','DEPTH','TIME'],mat_pden),
                'NPROF':(['LONGITUDE','TIME'],mat_nprof),
                'LATITUDE':(['LONGITUDE'],mat_lat),},
                coords={'TIME': python_time,'LONGITUDE':mat_lon,'DEPTH':mat_dpth})

dat_mat2=io.loadmat(datadir+'OAproduct/0715_OSNAP_All_OA_201407-201807_ArgoOnly.mat')
dat_mat2.keys()
## [('time', 'O'), ('lon', 'O'), ('lat', 'O'), ('dist', 'O'), ('depth', 'O'), ('t', 'O'), ('s', 'O'), ('d', 'O'), ('pd', 'O'), ('pt', 'O'), ('t_err', 'O'), ('s_err', 'O'), ('t_clim', 'O'), ('s_clim', 'O'), ('nprof', 'O'), ('DATA_FILE', 'O'), ('CREATION_DATE', 'O')])}

mat2_time=dat_mat2['ts_analysis_mm'][0][0][0].flatten()
mat2_lon=dat_mat2['ts_analysis_mm'][0][0][1].flatten()
mat2_lat=dat_mat2['ts_analysis_mm'][0][0][2].flatten()
mat2_dpth=dat_mat2['ts_analysis_mm'][0][0][4].flatten()
mat2_sal=dat_mat2['ts_analysis_mm'][0][0][6]
mat2_pden=dat_mat2['ts_analysis_mm'][0][0][8]
mat2_ptmp=dat_mat2['ts_analysis_mm'][0][0][9]
mat2_nprof=dat_mat2['ts_analysis_mm'][0][0][-3]

python2_time = array([datetime.datetime.fromordinal(int(mm)) + datetime.timedelta(days=mm%1) - datetime.timedelta(days = 366) for mm in mat2_time])


dat2=xr.Dataset({'PSAL':(['LONGITUDE','DEPTH','TIME'],mat2_sal),
                'PTMP':(['LONGITUDE','DEPTH','TIME'],mat2_ptmp),
                'PDEN':(['LONGITUDE','DEPTH','TIME'],mat2_pden),
                'NPROF':(['LONGITUDE','TIME'],mat2_nprof),
                'LATITUDE':(['LONGITUDE'],mat2_lat),},
                coords={'TIME': python2_time,'LONGITUDE':mat2_lon,'DEPTH':mat2_dpth})


def plotsec_all(dat,var,var2,vartit,elfolder):
        f,axx=subplots(1,1,figsize=(25,5))
        den=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat[var].mean(dim='TIME').T,univec[var2][1],cmap=univec[var2][2],extend='both')
        axx.contour(dat.LONGITUDE,dat.DEPTH,dat[var].mean(dim='TIME').T,levels=univec[var2][1][::2],colors='k')
        axx.set_ylabel('depth [m]')
        axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
        axx.set_ylim(3750,0)
        axx.set_xlabel('Longitude')
        colorbar(den)
        title('48 month mean '+vartit+' across OSNAP: '+elfolder)
        savefig(figdir+'1418_pden/'+elfolder+'/'+var2+'_mean.png',bbox_inches='tight')

        for ii in range(48):

                f,axx=subplots(1,1,figsize=(25,15))
                den=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat[var][:,:,ii].T,univec[var2][1],cmap=univec[var2][2],extend='both')
                axx.contour(dat.LONGITUDE,dat.DEPTH,dat[var][:,:,ii].T,levels=univec[var2][1][::2],colors='k')
                axx.set_ylabel('depth [m]')
                axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
                axx.set_ylim(3750,0)
                axx.set_xlabel('Longitude')
                colorbar(den)
                title(vartit+': '+str(dat.TIME[ii].values)[:10]+' ('+elfolder+')')
                savefig(figdir+'1418_pden/'+elfolder+'/'+var2+'sec_'+str(ii).zfill(2)+'.png',bbox_inches='tight',dpi=300)

plotsec_all(dat,'PDEN','pden','potential density','ArgoMooringGlider')
plotsec_all(dat2,'PDEN','pden','potential density','ArgoOnly')

58/1e3
