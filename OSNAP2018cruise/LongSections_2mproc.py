from AR30_funcs import *
##################################################################################################
####################################### LOAD ######################################################
##################################################################################################
ctd=xr.open_dataset('OSNAP2018cruise/data/CTD_2m_Leahproc.nc')

bathy={}
bathy['section 5']=io.loadmat(datadir+'ar30_bathy/sect5_long_clean.mat')
bathy['section 5']['dist']=bathy['section 5']['dist'].max()-bathy['section 5']['dist']
bathy['section 2']=io.loadmat(datadir+'ar30_bathy/sect2_long_clean.mat')
bathy['section 1']=io.loadmat(datadir+'ar30_bathy/sect1_clean.mat')
bathy['section 3']=io.loadmat(datadir+'ar30_bathy/sect3_clean.mat')

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
    # plot(vma.lon,vma.lat,'k.')
    xlabel('lon')
    ylabel('lat')
    savefig(figdirplus+'Map_LongSecs_Ar30.png',bbox_inches='tight')

quickmap()


def plot_TS_core():
    f,axx=subplots(1,1,figsize=(5,4))
    distch=70
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

uni['o2']['vmin']=6
uni['o2']['vmax']=8

for ss in seclist:
    print(ss,seclab[ss])
distvec

plotvar('o2','Oxygen')


plotvar('sal','Salinity')

plotvar('tmp','Temperature')

##################################################################################################
##################################################################################################
####################################### WANT TO HAVE A CLOSER LOOK #################################
####################################### AT THE LS AREA T-S-O SPACE #################################
######################################## IN EACH SECTION/SPATIAL  #################################
##################################################################################################
##################################################################################################

ctd

distmat,prsmat=meshgrid(ctd.dist.values,ctd.prs.values)

ctd['distmat']=(('prs','sta'),distmat)

help(hexbin)

def hex_TS(sec):
    f,[ax1,ax2]=subplots(1,2,sharex=True,sharey=True,figsize=(14,5))
    p1=ax1.hexbin(ctd.sal.where(ctd.sal>34.8).where(ctd.tmp<6).where(ctd.dist<120)[:,seclab[sec]].values.flatten(),
           ctd.tmp.where(ctd.sal>34.8).where(ctd.tmp<6).where(ctd.dist<120)[:,seclab[sec]].values.flatten(),
           C=ctd.distmat.where(ctd.sal>34.8).where(ctd.tmp<6).where(ctd.dist<120)[:,seclab[sec]].values.flatten(),
           cmap=cm.rainbow,mincnt=1)
    colorbar(p1,ax=ax1,label='distance [km]')
    p2=ax2.hexbin(ctd.sal.where(ctd.sal>34.8).where(ctd.tmp<6)[:,seclab[sec]].values.flatten(),
          ctd.tmp.where(ctd.sal>34.8).where(ctd.tmp<6)[:,seclab[sec]].values.flatten(),
          C=ctd.o2.where(ctd.sal>34.8).where(ctd.tmp<6)[:,seclab[sec]].values.flatten(),
          cmap=cm.RdPu,mincnt=1)
    colorbar(p2,ax=ax2,label='Oxygen [mL/L]')
    suptitle(sec)
    xlabel('salinity')
    ylabel('temperature [$^\circ$C]')
    savefig(figdirplus+'TSO_dist_sec'+sec[-1]+'.png',bbox_inches='tight')


seclist
for sec in seclist:
    hex_TS(sec)
