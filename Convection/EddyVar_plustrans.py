from aux_funcs import *

[uIIW,dIIW,IIW,egic]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_xport/IIWtrans_direct.pickle','rb'))

osnap=pickle.load(open(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_full.pickle','rb'))
uIIW['osnap east trans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=-8).where(osnap.PDEN>=d1).where(osnap.PDEN<d2).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6
dIIW['osnap east trans']=-(osnap['VELO']*osnap['AREA']).where(osnap.LONGITUDE>-43).where(osnap.LONGITUDE<=-8).where(osnap.PDEN>=d2).where(osnap.PDEN<d3).sum(dim='LONGITUDE').sum(dim='DEPTH')/1e6


spowdic=pickle.load(open(datadir+'OSNAP2016recovery/pickles/spectral/scale_avg_power_2-24hour.pickle','rb'))
datadir
psi=xr.open_dataset(datadir+'OSNAP2016recovery/gridded/OSNAP_Streamfunction_201408_201604_2018.nc')

MOC_east=ones(21)
for ii in range(21):
    MOC_east[ii]=psi['T_EAST'][ii,:].max()
cvec=[cf5col,cf6col,m1col,m1col]
# from scipy import signal
#
# N  = 2    # Filter order
# # Wn = 0.03333333/48 # Cutoff frequency (30 days)
# Wn = 0.03333333/20
# B, A = sig.butter(N, Wn, output='ba')

def spowcomp(axx):
    # figure(figsize=(18,3))
    jj=0
    for moor in [5,6,8]:
        for key in spowdic[moor]:
            if (key==500) | (key==516) | (key==549):
                if moor==8:
                    labit='M1'
                else:
                    labit='CF'+str(moor)
                if jj==3:
                    labit=''
                axx.plot(spowdic[moor][key]['date'],spowdic[moor][key]['spow'],alpha=0.8,color=cvec[jj],label=labit)
                # axx.plot(spowdic[moor][key]['date'],sig.filtfilt(B,A,spowdic[moor][key]['spow']),label=labit,linewidth=3,color=cvec[jj])
                # axhline(spowdic[moor][key]['sig'],color=cvec[jj],label='',linestyle='--',alpha=0.5)
                jj+=1
        axx.legend(loc=(1.01,0.35))
        axx.set_ylim(0,0.00325)
        axx.set_yticks(arange(0,0.003,0.001))
        axx.set_yticklabels(['0','1 x 10$^{-3}$','2 x 10$^{-3}$'])
        axx.set_xlim(datetime.datetime(2014,9,15),datetime.datetime(2016,7,15))
        axx.set_ylabel('High frequency\nspeed variance [$m^4/s^4$]')

def plot_trans(axi):
    hfa=0.5
    N  = 2    # Filter order
    Wn = 0.03333333 # Cutoff frequency (30 days)
    B, A = sig.butter(N, Wn, output='ba')
    axi.plot(egic['trans'].date,-egic['trans'],color='k',alpha=hfa)
    axi.plot(egic['trans'].date,-sig.filtfilt(B,A,egic['trans']),color='k',linewidth=3,alpha=0.7)
    axi.set_ylim(8,-37)
    axi.set_ylabel('Boundary current\ntransport [Sv]')
    axx=axi.twinx()
    axx.plot(uIIW['trans'].date,-uIIW['trans'],color=uppercol,alpha=hfa)
    # axx.plot(osnap.TIME,uIIW['osnap east trans'],'o-',color=uppercol,)
    axx.plot(uIIW['trans'].date,-sig.filtfilt(B,A,uIIW['trans']),color=uppercol,linewidth=3)
    # axx.plot(osnap.TIME,dIIW['osnap east trans'],'o-',color=deepcol,)
    axx.plot(dIIW['trans'].date,-dIIW['trans'],color=deepcol,alpha=hfa)
    axx.plot(dIIW['trans'].date,-sig.filtfilt(B,A,dIIW['trans']),color=deepcol,linewidth=3)
    axx.set_ylabel('IIW transport [Sv]')
    axx.set_ylim(0,-19)



    axx.set_xticklabels('')


# grden=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_dengrid/grden_cf3-ooi.pickle','rb'),)

def plot_Fig3_OLD():
    f=figure(figsize=(10,8))

    gs0=gridspec.GridSpec(3,1, hspace = 0.1)
    ax0=plt.subplot(gs0[0])
    ax1=plt.subplot(gs0[1])
    ax2=plt.subplot(gs0[2])

    lw=3
    ax0.plot(psi.TIME,MOC_east,color='k',linewidth=lw)
    axii=ax0.twinx()
    axii.plot(osnap.TIME,uIIW['osnap east trans'],color=uppercol,linewidth=lw)
    axii.plot(osnap.TIME,dIIW['osnap east trans'],color=deepcol,linewidth=lw)
    ax0.text(datetime.datetime(2016,4,1),8,'upper IIW',color=uppercol,fontsize=16)
    ax0.text(datetime.datetime(2016,4,1),-5,'deep IIW',color=deepcol,fontsize=16)
    ax0.set_xticklabels('')
    ax0.set_ylim(-10,25)
    ax0.axhline(0,color='k')
    axii.set_ylim(-5,12.5)
    ax0.set_ylabel('OSNAP East\nOverturning transport [Sv]')
    axii.set_ylabel('IIW contribution\nto overturning [Sv]')

    plot_trans(ax1)
    spowcomp(ax2)

    for axx in [ax0,ax1,ax2]:
        axx.set_xlim([datetime.datetime(2014,7,15),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        if (axx==ax2):
            axx.xaxis.set_minor_formatter(monthFMT)
            axx.xaxis.set_major_formatter(yearFMT)

    # for axx in [ax0,ax1,ax2]:
    #     for yy in [2015,2016]:
    #         axx.axvspan(datetime.datetime(yy-1,11,1),datetime.datetime(yy,5,1),color='C0',alpha=0.2)



plot_Fig3_OLD()



def plot_eddyvar_plustrans_nosec():
    f=figure(figsize=(10,5))

    gs0=gridspec.GridSpec(2,1, hspace = 0.1)
    ax1=plt.subplot(gs0[0])
    ax2=plt.subplot(gs0[1])

    plot_trans(ax1)
    spowcomp(ax2)

    for axx in [ax1,ax2]:
        axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        if (axx==ax2):
            axx.xaxis.set_minor_formatter(monthFMT)
            axx.xaxis.set_major_formatter(yearFMT)

    for axx in [ax1,ax2]:
        for yy in [2015,2016]:
            axx.axvspan(datetime.datetime(yy-1,11,1),datetime.datetime(yy,5,1),color='C0',alpha=0.2)

    savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans_nosec.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans_nosec.png',bbox_inches='tight')


plot_eddyvar_plustrans_nosec()

XXXXXXXXXXXXXXXXX

#
# cfgrid=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))
#
#
# def plotdensec(dl,axx):
#     dmin=datetime.datetime.strptime(dl, '%Y-%m-%d')-datetime.timedelta(days=15)
#     dmax=datetime.datetime.strptime(dl, '%Y-%m-%d')+datetime.timedelta(days=15)
#     velslice=cfgrid['across track velocity'].sel(date=slice(dmin,dmax))
#     vel=axx.contourf(cfgrid.distance,cfgrid.depth,velslice.mean(dim='date').T,13,vmin=-0.6,vmax=0,cmap=cm.Blues_r,extend='both')
#     # axx.contour(cfgrid.distance,cfgrid.depth,velslice.mean(dim='date').T,arange(-0.6,0.1,0.1),colors='w')#arange(-0.6,-0.1,0.1),colors='C0',linewidths=3)
#     axx.contour(cfgrid.distance,cfgrid.depth,cfgrid['potential density'].sel(date=slice(dmin,dmax)).mean(dim='date').T,levels=dbnds,colors='k',linewidths=4)
#     axx.contour(cfgrid.distance,cfgrid.depth,cfgrid['potential density'].sel(date=slice(dmin,dmax)).mean(dim='date').T,levels=arange(27.61,d3,0.015),colors='k')
#     axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
#     axx.axvline(distvec[4],color=cf5col,linewidth=4)
#     axx.axvline(distvec[5],color=cf6col,linewidth=4)
#     axx.axvline(distvec[7],color=m1col,linewidth=4)
#     # axx.axvline(oom_dist,color=ooicol,linewidth=4)
#     axx.set_ylim(1500,0)
#     axx.set_xlim(30,100)
#     # axx.set_xticks(arange(-42.5,-41,0.5))
#     axx.set_yticks(arange(0,1600,500))
#     return velslice.date.values,vel
#
# def plot_eddyvar_plustrans():
#     f=figure(figsize=(12,10))
#     outer = gridspec.GridSpec(2, 1, height_ratios = [2,1],hspace=0.35)
#     gs0=gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec = outer[0], hspace = 0.1)
#     ax1=plt.subplot(gs0[0])
#     ax2=plt.subplot(gs0[1])
#     #make nested gridspec at bottom for examples of vel/den sections when eddy activity high
#     gs1 = gridspec.GridSpecFromSubplotSpec(1, 11, subplot_spec = outer[1])
#     ax3=plt.subplot(gs1[1:5])
#     ax4=plt.subplot(gs1[6:10])
#
#     plot_trans(ax1)
#     spowcomp(ax2)
#
#     for axx in [ax1,ax2]:
#         axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
#         axx.xaxis.set_major_locator(years)
#         axx.xaxis.set_minor_locator(threemonth)
#         if (axx==ax2):
#             axx.xaxis.set_minor_formatter(monthFMT)
#             axx.xaxis.set_major_formatter(yearFMT)
#
#     date1='2015-8-15'
#     date2='2016-1-15'
#     ov1,v1=plotdensec(date1,ax3)
#     ov2,v2=plotdensec(date2,ax4)
#     cax1=f.add_axes([0.85,0.125,0.02,0.2])
#     colorbar(v2,ax=ax4,cax=cax1,ticks=arange(-0.6,0,0.2),label='velocity [m/s]')
#     ax4.set_yticklabels('')
#     ax3.set_title('August 2015',fontsize=15)
#     ax4.set_title('January 2016',fontsize=15)
#     ax3.set_ylabel('depth [m]')
#     f.text(0.4, 0.05, 'distance [km]', va='center', rotation='horizontal',fontsize=14)
#     for ii in [ax1,ax2]:
#         for vv in [date1,date2]:
#             [ii.axvline(vv,color='b',linewidth=4,alpha=0.4)]
#     savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans.pdf',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/paperfigs/EddyVar_plustrans.png',bbox_inches='tight')
#
#
# plot_eddyvar_plustrans()

#
#
# def plotdensec(dl,axx):
#     dmin=datetime.datetime.strptime(dl, '%Y-%m-%d')-datetime.timedelta(days=20)
#     dmax=datetime.datetime.strptime(dl, '%Y-%m-%d')+datetime.timedelta(days=20)
#     velslice=osnap.VELO.sel(TIME=slice(dmin,dmax))
#     vel=axx.contourf(osnap.LONGITUDE,osnap.DEPTH,velslice.mean(dim='TIME'),31,vmin=-0.5,vmax=0.5,cmap=cm.RdBu_r,extend='both')#arange(-0.6,-0.1,0.1),colors='C0',linewidths=3)
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.sel(TIME=slice(dmin,dmax)).mean(dim='TIME'),levels=dbnds,colors='k',linewidths=4)
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.sel(TIME=slice(dmin,dmax)).mean(dim='TIME'),levels=arange(27.61,d3,0.015),colors='k')
#     axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),2500*ones(len(osnap_bathy['bathy'].flatten())),color='k',zorder=22)
#     axx.axvline(CFlon[4],color=cf5col,linewidth=4)
#     axx.axvline(CFlon[5],color=cf6col,linewidth=4)
#     axx.axvline(CFlon[7],color=m1col,linewidth=4)
#     axx.axvline(ooi_lon['fla'],color=ooicol,linewidth=4)
#     axx.set_ylim(1500,0)
#     axx.set_xlim(-42.5,-41)
#     axx.set_xticks(arange(-42.5,-41,0.5))
#     axx.set_yticks(arange(0,1600,500))
#     return velslice.TIME.values,vel
