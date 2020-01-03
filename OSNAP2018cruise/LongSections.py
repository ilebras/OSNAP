from AR30_funcs import *
##################################################################################################
####################################### LOAD ######################################################
##################################################################################################

ctd=xr.open_dataset('OSNAP2018cruise/data/CTD_1mbin.nc')

bathy={}
bathy['section 5']=io.loadmat(datadir+'ar30_bathy/sect5_long_clean.mat')
bathy['section 5']['dist']=bathy['section 5']['dist'].max()-bathy['section 5']['dist']
bathy['section 2']=io.loadmat(datadir+'ar30_bathy/sect2_long_clean.mat')
bathy['section 1']=io.loadmat(datadir+'ar30_bathy/sect1_clean.mat')
bathy['section 3']=io.loadmat(datadir+'ar30_bathy/sect3_clean.mat')

ctd

vma=xr.open_dataset('OSNAP2018cruise/data/VMADCP_os.nc')

############ This merging process is not quite working -- will have to fix!!
dat=xr.merge([vma,ctd],compat='override')

figdirplus=figdir+'AR30/1912_LongSecFocus/'
####################################### PLOT ######################################################
##################################################################################################


# Plot the long sections only: 1,2,3,5:
seclist=['section 1','section 2','section 3','section 5']

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


for ss in seclist:
    plot(dat.lon[seclab[ss]],dat.lat[seclab[ss]],'x',label=ss)

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
