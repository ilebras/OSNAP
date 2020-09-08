from aux_funcs_2020 import *

################################################################################################################################
#################
# # # #### load osnap data

figdir

dat=xr.open_dataset(datadir+'OSNAP_Gridded_Fields_201408_201805_2020.nc')
dat['PDEN']=dat['SIGMA_THETA']-1e3

colors = [(255,237,160),(127,205,187),(5,112,176),(110,1,107)]
pden_cmap=make_cmap(colors,position=[0,0.666666666666,0.833333333333,1],bit=True)

univec={}
univec['pden']=['potential density',array([26.75,27,27.125,27.25,27.375,27.5,27.55,27.6,27.64,27.68,27.695,27.71,27.725,27.74,27.755,27.77,27.785,27.8,27.815,27.83,27.845,27.86,27.875,27.89]),pden_cmap,[27,27.5,27.68,27.74,27.8,27.86],'[kg/m$^3$]']


def plotsec_all(dat,var,var2,vartit):
        f,axx=subplots(1,1,figsize=(25,5))
        den=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat[var].mean(dim='TIME'),univec[var2][1],cmap=univec[var2][2],extend='both')
        axx.contour(dat.LONGITUDE,dat.DEPTH,dat[var].mean(dim='TIME'),levels=univec[var2][1][::2],colors='k')
        axx.set_ylabel('depth [m]')
        # axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k')
        axx.set_ylim(3750,0)
        axx.set_xlabel('Longitude')
        colorbar(den)
        title('48 month mean '+vartit+' across OSNAP')
        savefig(figdir+'SI/osnap_wide/'+var2+'_mean.png',bbox_inches='tight')

        for ii in range(48):
                f,axx=subplots(1,1,figsize=(25,15))
                den=axx.contourf(dat.LONGITUDE,dat.DEPTH,dat[var][ii,:,:],univec[var2][1],cmap=univec[var2][2],extend='both')
                axx.contour(dat.LONGITUDE,dat.DEPTH,dat[var][ii,:,:],levels=univec[var2][1][::2],colors='k')
                axx.contour(dat.LONGITUDE,dat.DEPTH,dat[var][ii,:,:],levels=[27.56],colors='r',linewidth=3)
                axx.set_ylabel('depth [m]')
                axx.set_ylim(3750,0)
                axx.set_xlabel('Longitude')
                colorbar(den)
                title(vartit+': '+str(dat.TIME[ii].values)[:10])
                savefig(figdir+'SI/osnap_wide/'+var2+'_'+str(ii)+'.png',bbox_inches='tight')

plotsec_all(dat,'PDEN','pden','potential density')
