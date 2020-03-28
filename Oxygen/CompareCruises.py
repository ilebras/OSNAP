from AR30_funcs import *
from firstfuncs_1618 import *
##################################################################################################
####################################### LOAD 2018 #################################################
##################################################################################################
ctd_18_all=xr.open_dataset('OSNAP2018cruise/data/CTD_2m_Leahproc.nc')

ctd_18=ctd_18_all.sel(sta=seclab['section 2'])
ctd_18['o2']=ctd_18['o2']*44.661

bathy={}
bathy['section 2']=io.loadmat(datadir+'Shipboard/AR30_2018/ar30_bathy/sect2_long_clean.mat')
#Leave velocity out of it for now, but I DO have LADCP for all past ones, I believe... will have to ask Penny if I use all those.
# vma=xr.open_dataset('OSNAP2018cruise/data/VMADCP_os.nc')
figdirplus='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Oxygen/DO_sections/'
##################################################################################################
####################################### LOAD PAST YEARS ###########################################
##################################################################################################
predat=xr.open_dataset(datadir+'Shipboard/gridded/Shipboard_2014x2_2016_wO2.nc').sel(prs=slice(10,3120))
##################################################################################################
####################### Reformat 2018 data and merge with past data ##############################
##################################################################################################
dat_18=ctd_18.copy()
dat_18

# dat_18=dat_18.drop('turner').drop('dendiff').drop('prsi')
dat_18=dat_18.drop('CT').drop('SA')

dat_18=dat_18.groupby_bins('prs',range(5,3130,10)).mean()
dat_18=dat_18.rename({'prs_bins':'prs'})
dat_18=dat_18.rename({'den':'pden'})

dat_18['prs']=predat['prs']

distcorr=sw.dist([dat_18.lat[0,0],CFlat[0]],[dat_18.lon[0,0],CFlon[0]])[0]

dist18=dat_18.dist[0,:].values-distcorr


dat_18['sta']=dist18
dat_18=dat_18.drop('dist')
dat_18=dat_18.rename({'sta':'dist'})

dat_18=dat_18.sortby('dist')

predat

dat_18=dat_18.groupby_bins('dist',arange(-2.5,725,5)).mean()
dat_18=dat_18.rename({'dist_bins':'dist'})

dat_18['dist']=arange(0,721,5)

dat_18=dat_18.assign_coords({'occ':'2018'})

grdat=xr.concat([predat,dat_18],dim='occ')

##################################################################################################
#################################### PLOT ########################################################
##################################################################################################
d1=1
d2=40
for occii in grdat.occ.values:
    figure()
    plot(grdat.o2.sel(occ=occii).isel(dist=slice(d1,d2)),grdat.prs)
    title(occii)
    ylim(3000,0)
    xlim(250,400)
    ylabel('pressure [db]')
    xlabel('uncalibrated DO concentration [$\mu$M/L]')
    savefig(figdirplus+'DO_conc_profs_'+occii+'.png',bbox_inches='tight')


figdirplus
colors = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(colors,position=[0,0.75,0.9,0.95,1],bit=True)


def plot_var(var,varname,var_cmap,vmini,vmaxi,vlab):
    g=grdat[var][:,d1:d2,:].plot(x='dist', y='prs', col='occ',vmin=vmini,vmax=vmaxi, col_wrap=2,cmap=var_cmap,figsize=(9,6),cbar_kwargs={'label': varname})
    ylim([3e3,0]);
    for ii, ax in enumerate(g.axes.flat):
        ax.set_title(grdat.occ.values[ii],fontsize=14)
        ax.fill_between(fullbathdist,fullbathbath,3000*ones(len(fullbathbath)),color='k',zorder=22)
        kcont=ax.contour(grdat.dist[d1:d2],grdat.prs,grdat[var][:,d1:d2,ii],vlab,colors='k')
        if 'Sal' in varname:
            ax.clabel(kcont,fmt='%1.1f')
        ax.axvline(175,color='k')
    axx=g.axes
    axx[1,0].set_xlabel('distance [km]')
    axx[1,1].set_xlabel('distance [km]')
    axx[0,0].set_ylabel('pressure [db]')
    axx[1,0].set_ylabel('pressure [db]')
    savefig(figdirplus+'Sections_'+var+'.png',bbox_inches='tight')



plot_var('sal','Salinity',sal_cmap,34,35.1,[34.9,35])

plot_var('o2','Oxygen [$\mu$M/L]',cm.rainbow,250,330,[])

plot_var('tmp','Temperature [$^\circ$C]',cm.RdYlBu_r,0,9,[])

plot_var('pden','Density, $\sigma_0$ [kg/m$^3$]',cm.YlGnBu,27.5,28,[])

# Plot the long sections only: 1,2,3,5:
seclist=['section 1','section 2','section 3','section 5']

plot(seclab['section 3'],'o')

plot(seclab['section 5'],'o')

def quickmap():
    for ss in seclist:
        plot(ctd.lon[seclab[ss]],ctd.lat[seclab[ss]],'o',label=ss)
        plot(ctd.lon[startsec[ss]],ctd.lat[startsec[ss]],'kx',label='')
        legend(loc=(1.05,0))
    plot(vma.lon,vma.lat,'k.')
    xlabel('lon')
    ylabel('lat')
    savefig(figdirplus+'Map_LongSecs_Ar30.png',bbox_inches='tight')

quickmap()

ctd[]
for ss in seclist:
    plot(ctd.lon[seclab[ss]],ctd.lat[seclab[ss]],'x',label=ss)


def plot_TS_core():
    f,axx=subplots(1,1,figsize=(5,4))
    distch=80
    for ii,ss in enumerate(seclist[::-1]):
        plot(ctd.sal[:,seclab[ss]].where(ctd.dist[seclab[ss]]<distch).values.flatten(),ctd.tmp[:,seclab[ss]].where(ctd.dist[seclab[ss]]<distch).values.flatten(),'o',label=ss,alpha=0.2,color='C'+str(3-ii))
        legend(loc=(1.05,0))
        den=axx.contour(salvec,tmpvec,pdenmat,colors='k',linewidths=2,levels=arange(27.5,28,0.1))
        axx.set_ylim(1.5,6.5)
        axx.set_xlim(34.825,34.975)
        xlabel('salinity')
        ylabel('temperature')
        title('Properties within '+str(distch)+'km of shore')
        savefig(figdirplus+'TS_intwater_core_LongSecs_Ar30.png',bbox_inches='tight')

plot_TS_core()


def plot_TS():
    f,axx=subplots(1,1,figsize=(5,4))
    for ii,ss in enumerate(seclist[::-1]):
        plot(ctd.sal[:,seclab[ss]].values.flatten(),ctd.tmp[:,seclab[ss]].values.flatten(),'o',label=ss,alpha=0.2,color='C'+str(3-ii))
        legend(loc=(1.05,0))
        den=axx.contour(salvec,tmpvec,pdenmat,colors='k',linewidths=2,levels=arange(27.5,28,0.1))
        axx.set_ylim(1.5,6.5)
        axx.set_xlim(34.825,34.975)
        xlabel('salinity')
        ylabel('temperature')
        savefig(figdirplus+'TS_intwater_LongSecs_Ar30.png',bbox_inches='tight')


plot_TS()


def contourit(axx,var,sect):
    conty=axx.contourf(sort(ctd.dist[seclab[sect]]),ctd.prs,ctd[var][:,seclab[sect]].sortby(ctd.dist),51,vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
    denden=axx.contour(sort(ctd.dist[seclab[sect]]),ctd.prs,ctd['den'][:,seclab[sect]].sortby(ctd.dist),levels=arange(27.5,28,0.1),colors='k',linewidths=2)
    axx.fill_between(bathy[sect]['dist'].flatten(),bathy[sect]['bath'].flatten(),3e3,color='grey')
    axx.text(40,2500,sect,fontsize=16,color='white')
    return conty

def plotvar(var,tit):
    f,axx=subplots(2,2,sharex=True,sharey=True,figsize=(9,6))
    contourit(axx[0,0],var,'section 3')
    contourit(axx[1,0],var,'section 5')
    contourit(axx[0,1],var,'section 2')
    conty=contourit(axx[1,1],var,'section 1')
    caxy=f.add_axes([0.93,0.2,0.02,0.6])
    colorbar(conty,cax=caxy,label=tit)
    axx[1,1].set_ylim(2750,0)
    axx[1,1].set_xlim(35,120)
    f.text(0.5, -0.01, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.01, 0.5, 'pressure [db]', va='center', rotation='vertical',fontsize=16)
    axx[0,0].set_title('West of Greenland',fontsize=16)
    axx[0,1].set_title('East of Greenland',fontsize=16)
    savefig(figdirplus+'Section_'+var+'_LongSecs_Ar30.png',bbox_inches='tight')


uni['o2']={'cmap':cm.rainbow,'vmin':1024,'vmax':1030}

plotvar('o2','Oxygen')

plotvar('sal','Salinity')


plotvar('tmp','Temperature')

def contour_u(axx,var,sect):
    conty=axx.pcolor(sort(dat.dist[seclab[sect]]),dat.prs_u,dat[var][:,seclab[sect]].sortby(dat.dist),vmin=0,vmax=0.5)
    # denden=axx.contour(sort(dat.dist[seclab[sect]]),ctd.prs,ctd['den'][:,seclab[sect]].sortby(ctd.dist),levels=arange(27.5,28,0.1),colors='k',linewidths=2)
    axx.fill_between(bathy[sect]['dist'].flatten(),bathy[sect]['bath'].flatten(),3e3,color='grey')
    axx.text(2,475,sect,fontsize=16,color='white')
    return conty


def plotvar_u(var,tit):
    f,axx=subplots(2,2,sharex=True,sharey=True,figsize=(9,6))
    contour_u(axx[0,0],var,'section 3')
    contour_u(axx[1,0],var,'section 5')
    contour_u(axx[0,1],var,'section 2')
    conty=contour_u(axx[1,1],var,'section 1')
    caxy=f.add_axes([0.93,0.2,0.02,0.6])
    colorbar(conty,cax=caxy,label=tit)
    axx[1,1].set_ylim(500,0)
    axx[1,1].set_xlim(0,120)
    f.text(0.5, -0.01, 'distance [km]', ha='center',fontsize=16)
    f.text(-0.01, 0.5, 'pressure [db]', va='center', rotation='vertical',fontsize=16)
    axx[0,0].set_title('West of Greenland',fontsize=16)
    axx[0,1].set_title('East of Greenland',fontsize=16)
    savefig(figdirplus+'Section_SpeedTest_LongSecs_Ar30.png',bbox_inches='tight')


dat['speed']=sqrt(dat['u']**2+dat['v']**2)

plotvar_u('speed','speed')
