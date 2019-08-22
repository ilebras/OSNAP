from aux_funcs import *

#mooring resolution gridded CF data
dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1808lpfilt.pickle','rb'))

ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb'))

##OOI merged density profiles (sfc mooring and flanking mooring A)
oom=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_denmerged_xray.pickle','rb'))

[dendat,ooi_dendat,oom_dendat]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_dengrid/dendats_cf-ooi-oom.pickle','rb'))

grden=pickle.load(open(datadir+'OSNAP2016recovery/pickles/convection_dengrid/grden_cf3-ooi.pickle','rb'))


sal_uIIW=(dendat['sal'][:,p1:p2,:]*dendat['thickness'][:,p1:p2,:]).sum(dim='den')/dendat['thickness'][:,p1:p2,:].sum(dim='den')
sal_dIIW=(dendat['sal'][:,p2:p3,:]*dendat['thickness'][:,p2:p3,:]).sum(dim='den')/dendat['thickness'][:,p2:p3,:].sum(dim='den')

def plot_Fig2():
    f,axx=subplots(1,2,sharex=True,sharey=True,figsize=(14,5))
    grden_3d=grden.resample(date='3D').mean(dim='date')
    print(grden_3d.pden[:,1:-1,:])
    axx[0,0].contourf(grden_3d.date.values,range(5),(grden_3d.pden.where(grden.pden>=d1).where(grden.pden<d2).sum(dim='depth'))[:,1:-1,:],31,cmap=cm.viridis)
    axx[0,0].set_title('upper IIW')
    axx[0,1].set_title('deep uIIW'),
    axx[0,1].contourf(grden_3d.date.values,range(5),(grden_3d.pden.where(grden.pden>=d2).where(grden.pden<d3).sum(dim='depth'))[:,1:-1,:],31,cmap=cm.viridis)
    colorbar(thick,ax=axx[0,1],label='Thickness [m]')
    # axx[1,0].contourf(sal_uIIW.date.values,sal_uIIW.distance[4:],sal_uIIW[4:,:],31,cmap=cm.YlGnBu_r)
    # axx[1,0].set_ylabel('distance [km]')
    # axx[0,0].set_ylabel('distance [km]')
    # sal=axx[1,1].contourf(sal_dIIW.date.values,sal_uIIW.distance[4:],sal_dIIW[4:,:],31,cmap=cm.YlGnBu_r)
    # colorbar(sal,ax=axx[1,1],label='Salinity')
    savefig(figdir+'MixedLayer/paperfigs/Fig2.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/Fig2.png',dpi=300)

plot_Fig2()



def Thick_Tseries(m1,c1,d1,d2,ax1,var,labit):
        hfa=0.3
        N  = 2    # Filter order
        Wn = 0.0333333333 # Cutoff frequency (50 days)
        B, A = sig.butter(N, Wn, output='ba')
        if m1=='ooi':
            oomdate=range(5,len(oom.date)-5)
        if var=='sal':
            if m1=='ooi':
                thick1=(oom_dendat[var][d1:d2,oomdate]*oom_dendat['thickness'][d1:d2,oomdate]).sum(dim='den')/oom_dendat['thickness'][d1:d2,oomdate].sum(dim='den')
            else:
                thick1=(dendat[var][m1,d1:d2,:]*dendat['thickness'][m1,d1:d2,:]).sum(dim='den')/dendat['thickness'][m1,d1:d2,:].sum(dim='den')
        elif var=='thickness':
            if m1=='ooi':
                thick1=oom_dendat['thickness'][d1:d2,oomdate].sum(dim='den')
            else:
                thick1=dendat['thickness'][m1,d1:d2,:].sum(dim='den')
        nonan=argwhere(~isnan(thick1.values)).flatten()
        f=interp1d(nonan,thick1[nonan])
        thickinterp=f(range(len(thick1)))
        tfilt=tfilt(B,A,thickinterp)
        tfilt[isnan(thick1.values)]=NaN
        tresamp=thick1.resample(date='2D').mean(dim='date')
        ax1.plot(tresamp.date,tresamp,label='',color=c1,alpha=hfa)
        if m1=='ooi':
            ax1.plot(oom_dendat.date[oomdate],tfilt,color=ooicol,linewidth=2,label=labit)
        else:
            ax1.plot(dendat.date,tfilt,color=c1,linewidth=2,label=labit)

def plotallthick(ax1,ax2,var):
        Thick_Tseries(4,cf5col,p1,p2,ax1,var,'CF5: Boundary current maximum')
        Thick_Tseries(5,cf6col,p1,p2,ax1,var,'CF6: Deep boundary current mooring')
        # Thick_Tseries(6,cf7col,p1,p2,ax11)
        Thick_Tseries(7,m1col,p1,p2,ax1,var,'M1: Offshore of the boundary current')

        Thick_Tseries(4,cf5col,p2,p3,ax2,var,'CF5: Boundary current maximum')
        Thick_Tseries(5,cf6col,p2,p3,ax2,var,'CF6: Deep boundary current mooring')
        # Thick_Tseries(6,cf7col,p2,p3,ax22)
        Thick_Tseries(7,m1col,p2,p3,ax2,var,'M1: Offshore of the boundary current')

        Thick_Tseries('ooi',ooicol,p1,p2,ax1,var,'OOI: Irminger gyre interior')
        Thick_Tseries('ooi',ooicol,p2,p3,ax2,var,'OOI: Irminger gyre interior')



lonvec_less=[CFlon[4],CFlon[5],CFlon[-1],ooi_lon['fla']]
colorvec=[cf5col,cf6col,m1col,ooicol]

moorlab=['CF5: Boundary current maximum','CF6: Deep boundary current mooring','M1: Offshore of the boundary current','OOI: Irminger gyre interior']


def SalEvol():

    gs1 = gridspec.GridSpec(2, 2,wspace=0.1)

    f=figure(figsize=(20,10))

    ax1 = plt.subplot(gs1[0,0])
    ax2 = plt.subplot(gs1[0,1])

    plotallthick(ax1,ax2,'thickness')

    ax11 = plt.subplot(gs1[1,0])
    ax22 = plt.subplot(gs1[1,1])

    plotallthick(ax11,ax22,'sal')

    ax11.legend(loc=(-0.1,-0.5),ncol=4,fontsize=14,markerscale=3)

    ax1.set_ylim(0,1200)
    ax1.set_yticks(range(0,1400,300))
    ax2.set_ylim(0,2400)
    ax2.set_yticks(range(0,2500,600))
    ax1.yaxis.grid('on')
    ax2.yaxis.grid('on')

    ax11.set_ylim(34.88,34.97)
    ax22.set_ylim(34.88,34.97)
    ax11.yaxis.grid('on')
    ax22.yaxis.grid('on')


    ax1.set_ylabel('Layer thickness [m]')
    ax11.set_ylabel('Salinity')

    # ax111.set_ylabel('Transport [Sv]')
    ax1.set_title('upper IIW',fontsize=18)
    ax2.set_title('deep IIW',fontsize=18)

    for axx in [ax1,ax2,ax11,ax22]:
        # axx.xaxis.tick_top()
        axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        if (axx==ax11) | (axx==ax22):
            axx.xaxis.set_minor_formatter(monthFMT)
            axx.xaxis.set_major_formatter(yearFMT)

    ax1.set_xticklabels('')
    ax2.set_xticklabels('')


    savefig(figdir+'MixedLayer/paperfigs/ThickSalEvol_linesonly.png',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/ThickSalEvol_linesonly.pdf',bbox_inches='tight')


SalEvol()






# def TS_tvar(axx):
#     cmapo = matplotlib.cm.get_cmap('GnBu')
#     for ii in range(len(date_list)-1):
#         tsm,den=TSmrange(ooi['salinity'],ooi['temperature'],axx,date_list[ii],date_list[ii+1],cmapo(20*(ii+5)),labvec[ii])
#         # tsm,den=TSmrange(oom['sal'],oom['ptmp'],axx,date_list[ii],date_list[ii+1],cmapo(20*(ii+5)),'')
#     axx.set_xlim(34.87,34.95)
#     axx.set_ylim(3,4)
#     axx.set_xlabel('Salinity')
#     axx.set_ylabel('Theta [$^\circ$C]')
#     legend(loc=(1.1,0.2))


#
# def veldensec(axx):
#     vel=axx.contourf(osnap.LONGITUDE,osnap.DEPTH,osnap.VELO.mean(dim='TIME'),31,vmin=-0.4,vmax=0.4,cmap=cm.RdBu_r,extend='both')
#     colorbar(vel,ax=axx,ticks=arange(-0.4,0.5,0.2),label='velocity [m/s]')
#     axx.contour(osnap.LONGITUDE,osnap.DEPTH,osnap.PDEN.mean(dim='TIME'),levels=dbnds,colors='k',linewidths=3)
#     axx.fill_between(osnap_bathy['lon'].flatten(),-osnap_bathy['bathy'].flatten(),[4000]*len(osnap_bathy['lon'].flatten()),color='k',zorder=100)
#     xlim(-42.5,-39.5)
#     ylim(1500,0)
#     xlabel('Longitude [$^\circ$W]')
#     ylabel('depth [m]')
#     [plot([mm,mm],[2000,0],color=colorvec[ii],linewidth=3,label=moorlab[ii]) for ii,mm in enumerate(lonvec_less)]
#     legend(loc=(1.7,0.2),fontsize=14)
#     return vel
#
# def TSmrange(salvar,tmpvar,axx,da1,da2,col,labit):
#     tsm=axx.plot(salvar.sel(date=slice(da1,da2)).mean(dim='date').values.flatten(),tmpvar.sel(date=slice(da1,da2)).mean(dim='date').values.flatten(),
#                  color=col,linewidth=3,label=labit)
#     den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[d1,d2,d3])
#     return tsm,den


# date_list = [datetime.datetime(2014,9,15) + datetime.timedelta(days=120*x) for x in range(0, 6)]
# date_mid = [datetime.datetime(2014,11,15) + datetime.timedelta(days=120*x) for x in range(0, 5)]
#
# labvec=['November 2014','March 2015','July 2015','November 2015','March 2016']
#
#
# def TS_tvar(axx):
#     cmapo = matplotlib.cm.get_cmap('GnBu')
#     for ii in range(len(date_list)-1):
#         tsm,den=TSmrange(ooi['salinity'],ooi['temperature'],axx,date_list[ii],date_list[ii+1],cmapo(20*(ii+5)),labvec[ii])
#         # tsm,den=TSmrange(oom['sal'],oom['ptmp'],axx,date_list[ii],date_list[ii+1],cmapo(20*(ii+5)),'')
#     axx.set_xlim(34.87,34.95)
#     axx.set_ylim(3,4)
#     axx.set_xlabel('Salinity')
#     axx.set_ylabel('Theta [$^\circ$C]')
#     legend(loc=(1.1,0.2))
#
# def SalEvol():
#     #make outer gridspec
#     outer = gridspec.GridSpec(2, 1, height_ratios = [2, 1.1],hspace=0.5)
#     #make nested gridspecs
#     gs1 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec = outer[0],wspace=0.1)
#     gs2 = gridspec.GridSpecFromSubplotSpec(1, 10, subplot_spec = outer[1], hspace = 1)
#
#     f=figure(figsize=(20,8))
#
#     ax1 = plt.subplot(gs1[0,0])
#     ax2 = plt.subplot(gs1[0,1])
#
#     plotallthick(ax1,ax2,'thickness')
#
#     ax11 = plt.subplot(gs1[1,0])
#     ax22 = plt.subplot(gs1[1,1])
#
#     plotallthick(ax11,ax22,'sal')
#
#     ax3 = plt.subplot(gs2[:2])
#     vel=veldensec(ax3)
#     ax2.set_title('Mean velocity and density boundaries')
#     # cax1=f.add_axes([0.95,0.55,0.008,0.32])
#     # colorbar(vel,ax=ax3,cax=cax1,ticks=arange(-0.4,0.5,0.2))
#
#     ax33 = plt.subplot(gs2[6:8])
#     TS_tvar(ax33)
#     ax33.set_title('Evolution of dIIW salinity minimum')
#
#     ax1.set_ylim(0,1200)
#     ax1.set_yticks(range(0,1400,300))
#     ax2.set_ylim(0,2400)
#     ax2.set_yticks(range(0,2500,600))
#     ax1.yaxis.grid('on')
#     ax2.yaxis.grid('on')
#
#     ax11.set_ylim(34.88,34.97)
#     ax22.set_ylim(34.88,34.97)
#     ax11.yaxis.grid('on')
#     ax22.yaxis.grid('on')
#
#
#     ax1.set_ylabel('Layer thickness [m]')
#     ax11.set_ylabel('Salinity')
#
#     # ax111.set_ylabel('Transport [Sv]')
#     ax1.set_title('upper IIW',fontsize=18)
#     ax2.set_title('deep IIW',fontsize=18)
#
#     for axx in [ax1,ax2,ax11,ax22]:
#         # axx.xaxis.tick_top()
#         axx.set_xlim([datetime.datetime(2014,9,1),datetime.datetime(2016,7,15)])
#         axx.xaxis.set_major_locator(years)
#         axx.xaxis.set_minor_locator(threemonth)
#         if (axx==ax11) | (axx==ax22):
#             axx.xaxis.set_minor_formatter(monthFMT)
#             axx.xaxis.set_major_formatter(yearFMT)
#
#     ax1.set_xticklabels('')
#     ax2.set_xticklabels('')
#
#
#     savefig(figdir+'MixedLayer/paperfigs/ThickSalEvol_lines.png',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/paperfigs/ThickSalEvol_lines.pdf',bbox_inches='tight')
#
#
# SalEvol()


#####ARROWS
#
# import matplotlib.dates as md
# x0 = md.date2num(datetime.datetime(2015,3,1))
# y0 = 800
# xw = md.date2num(datetime.datetime(2015,5,10)) - x0
# yw = 0
# ax1.arrow(x0, y0, xw, yw, color='black',head_width=30,head_length=15,linewidth=2)
# x0 = md.date2num(datetime.datetime(2016,1,1))
# xw = md.date2num(datetime.datetime(2016,3,15)) - x0
# ax1.arrow(x0, y0, xw, yw, color='black',head_width=30,head_length=15,linewidth=2)
#
# x0 = md.date2num(datetime.datetime(2015,2,1))
# y0 = 600
# xw = md.date2num(datetime.datetime(2015,4,1)) - x0
# yw = 400
# ax2.arrow(x0, y0, xw, yw, color='black',head_width=0,linewidth=2)
# x0 = md.date2num(datetime.datetime(2016,1,1))
# xw = md.date2num(datetime.datetime(2016,3,15)) - x0
# ax2.arrow(x0, y0, xw, yw, color='black',head_width=30,head_length=15,linewidth=2)

# ax1.add_patch(arr)



#
#
# ###############################################################################################################################
# ###############################################################################################################################
# ############################################## TS before and after   ####################################################
# ###############################################################################################################################
# ###############################################################################################################################
#
#
# def TSmrange(salvar,tmpvar,axx,da1,da2,col,labit,shade='no'):
#
#     tsm=axx.plot(salvar.sel(date=slice(da1,da2)).mean(dim='date').values.flatten(),tmpvar.sel(date=slice(da1,da2)).mean(dim='date').values.flatten(),
#                  color=col,linewidth=3,label=labit)
#     den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[d1,d2,d3])
#
#     return tsm,den
#
#     # if shade=='yes':
#     #     smin=salvar.sel(date=slice(da1,da2)).min(dim='date').values.flatten()
#     #     smax=salvar.sel(date=slice(da1,da2)).max(dim='date').values.flatten()
#     #     # smin=(salvar.sel(date=slice(d1,d2)).mean(dim='date')-2*salvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
#     #     # smax=(salvar.sel(date=slice(d1,d2)).mean(dim='date')+2*salvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
#     #     sall=hstack((smin[::-1],smax))
#     #     tmin=tmpvar.sel(date=slice(da1,da2)).min(dim='date').values.flatten()
#     #     tmax=tmpvar.sel(date=slice(da1,da2)).max(dim='date').values.flatten()
#     #     # tmin=(tmpvar.sel(date=slice(d1,d2)).mean(dim='date')-2*tmpvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
#     #     # tmax=(tmpvar.sel(date=slice(d1,d2)).mean(dim='date')+2*tmpvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
#     #     tall=hstack((tmin[::-1],tmax))
#     #     p = plt.Polygon(np.column_stack((sall,tall)), facecolor=col, alpha=.25, edgecolor=col)
#     #     axx.add_artist(p)
#     #     axx.plot(smin,tmin,color=col)
#     #     axx.plot(smax,tmax,color=col)
#
# #
# # def TS_BAshaded():
# #     f,axx=subplots(1,2,sharex=True,sharey=True,figsize=(8,3))
# #     #Before
# #     b1='2014-9-10'
# #     b2='2014-9-30'
# #     l1,na=TSmrange(ooi_dendat['sal'],ooi_dendat['tmp'],axx[0],b1,b2,ooicol,'Irminger gyre interior, OOI')
# #     l2,na=TSmrange(dendat['sal'][-1,:,:],dendat['tmp'][-1,:,:],axx[0],b1,b2,m1col,'Offshore of the boundary current, M1',shade='yes')
# #     l3,na=TSmrange(dendat['sal'][4,:,:],dendat['tmp'][4,:,:],axx[0],b1,b2,cf5col,'Boundary current maximum, CF5',shade='yes')
# #
# #     a1='2015-4-10'
# #     a2='2015-4-30'
# #     TSmrange(ooi_dendat['sal'],ooi_dendat['tmp'],axx[1],a1,a2,ooicol,'')
# #     TSmrange(dendat['sal'][-1,:,:],dendat['tmp'][-1,:,:],axx[1],a1,a2,m1col,'',shade='yes')
# #     na,d3=TSmrange(dendat['sal'][4,:,:],dendat['tmp'][4,:,:],axx[1],a1,a2,cf5col,'',shade='yes')
# #
# #
# #     axx[0].text(34.955,4.45,'uLSW',fontsize=14)
# #     axx[0].text(34.95,3.85,'dLSW',fontsize=14)
# #
# #     t=axx[1].text(34.95,4.7,'27.68')
# #     t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
# #     t=axx[1].text(34.95,4.1,'27.74')
# #     t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
# #     t=axx[1].text(34.95,3.55,'27.8')
# #     t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
# #
# #     axx[0].set_ylabel('Potential temperature [$^\circ$C]',fontsize=12)
# #     f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=12)
# #     titf=14
# #     axx[0].set_title('September 10-30, 2014\n (before convection)',fontsize=titf)
# #     axx[1].set_title('April 10-30, 2015\n (after convection)',fontsize=titf)
# #     axx[0].legend(loc=(-0.5,-0.4),ncol=3,fontsize=11)
# #     savefig(figdir+'MixedLayer/paperfigs/TS_BeforeAfter.pdf',bbox_inches='tight')
# #     savefig(figdir+'MixedLayer/paperfigs/TS_BeforeAfter.png',bbox_inches='tight')
# #
# # TS_BAshaded()
#
# ###############################################################################################################################
# ###############################################################################################################################
# ############################################## Layer thickness evolution   ####################################################
# ###############################################################################################################################
# ###############################################################################################################################
#
#
# # want to plot full density, not just these bounds, also want to merge properly for geostrophy calc etc
# # #make an xarray which charters the depths of the bounding isopycnals at all cf + ooi
# # cf_bnds=xr.concat([dendat.depth[4:,p1,:],dendat.depth[4:,p2,:],dendat.depth[4:,p3,:]],dim='den_bnd')
# # ooi_dpth_daily=ooi_dendat.depth.resample(date='1D').mean(dim='date')
# # ooi_bnds=xr.concat([ooi_dpth_daily[p1,:],ooi_dpth_daily[p2,:],ooi_dpth_daily[p3,:]],dim='den_bnd')
# # all_bnds=xr.concat([cf_bnds,ooi_bnds],dim='distance')
#
#
# dlist=['2014-10-15','2015-1-15','2015-4-15','2015-10-1','2016-1-1','2016-4-1']
# dtitlist=['Oct 15 2014','Jan 15 2015','Apr 15 2015','Oct 1 2015','Jan 1 2016','Apr 1 2016']
# clist=['#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b']
#
# def TS_triptich(axx,b1,b2):
#     l1,na=TSmrange(ooi_dendat['sal'],ooi_dendat['tmp'],axx,b1,b2,ooicol,'Irminger gyre interior, OOI')
#     # l1,na=TSmrange(oom['sal'],oom['ptmp'],axx,b1,b2,oomcol,'Irminger gyre interior, OOI merged')
#     l2,na=TSmrange(dendat['sal'][-1,:,:],dendat['tmp'][-1,:,:],axx,b1,b2,m1col,'Offshore of the boundary current, M1')
#     # l3,na=TSmrange(dendat['sal'][6,:,:],dendat['tmp'][6,:,:],axx,b1,b2,'r','CF7')
#     l3,na=TSmrange(dendat['sal'][5,:,:],dendat['tmp'][5,:,:],axx,b1,b2,cf6col,'CF6')
#     l3,na=TSmrange(dendat['sal'][4,:,:],dendat['tmp'][4,:,:],axx,b1,b2,cf5col,'Boundary current maximum, CF5')
#     return l1,l2,l3
#
# def Thick_Tseries(m1,c1,d1,d2,ax1,var='thickness'):
#         N  = 2    # Filter order
#         Wn = 0.02 # Cutoff frequency (50 days)
#         B, A = sig.butter(N, Wn, output='ba')
#         hfa=0.3
#
#         if m1=='ooi':
#             # thick1=nansum(ooi_dendat[var][d1:d2,:],axis=0)
#             # ax1.plot(ooi_dendat.date,thick1,label='',color=c1,alpha=hfa)
#             # ax1.plot(ooi_dendat.date[::2],tfilt(B,A,thick1[::2]),color=c1,linewidth=3)
#             oomdate=range(5,len(oom.date)-5)
#             thick2=nansum(oom_dendat[var][d1:d2,:],axis=0)[oomdate]
#             ax1.plot(oom_dendat.date[oomdate],thick2,label='',color=ooicol,alpha=hfa)
#             ax1.plot(oom_dendat.date[oomdate][::2],tfilt(B,A,thick2[::2]),color=ooicol,linewidth=3)
#         else:
#             thick1=nansum(dendat[var][m1,d1:d2,:],axis=0)
#             ax1.plot(dendat.date,thick1,alpha=hfa,color=c1)
#             ax1.plot(dendat.date,sig.filtfilt(B,A,thick1),color=c1,linewidth=3)
#
#
#
# def ThickEvol_wTS():
#
#     f=figure(figsize=(20,12))
#     gs = gridspec.GridSpec(4, 6)
#
#     gs.update(hspace=0.5,wspace=0.3)
#
#     ax1 = plt.subplot(gs[0, :3])
#     ax2 = plt.subplot(gs[0, 3:])
#
#     Thick_Tseries('ooi',ooicol,p1,p2,ax1)
#     Thick_Tseries('ooi',ooicol,p2,p3,ax2)
#
#     ax11 = plt.subplot(gs[1, :3])
#     ax22 = plt.subplot(gs[1, 3:])
#
#
#     Thick_Tseries(4,cf5col,p1,p2,ax11)
#     Thick_Tseries(5,cf6col,p1,p2,ax11)
#     # Thick_Tseries(6,cf7col,p1,p2,ax11)
#     Thick_Tseries(7,m1col,p1,p2,ax11)
#
#     Thick_Tseries(4,cf5col,p2,p3,ax22)
#     Thick_Tseries(5,cf6col,p2,p3,ax22)
#     # Thick_Tseries(6,cf7col,p2,p3,ax22)
#     Thick_Tseries(7,m1col,p2,p3,ax22)
#
#     ax1.set_ylabel('Layer thickness [m]')
#     ax11.set_ylabel('Layer thickness [m]')
#     ax1.set_title('upper LSW \n',fontsize=18)
#     ax2.set_title('deep LSW \n',fontsize=18)
#
#     ax1.set_ylim(0,1200)
#     ax2.set_ylim(300,2000)
#     ax11.set_ylim(0,1100)
#     ax22.set_ylim(0,1100)
#
#     for axx in [ax1,ax2,ax11,ax22]:
#         axx.xaxis.tick_top()
#         axx.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
#         axx.xaxis.set_major_locator(years)
#         axx.xaxis.set_minor_locator(threemonth)
#         if (axx==ax1) | (axx==ax2):
#             axx.xaxis.set_minor_formatter(monthFMT)
#             axx.xaxis.set_major_formatter(yearFMT)
#         [axx.axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]
#
#     ax11.set_xticklabels('')
#     ax22.set_xticklabels('')
#
#     for ii,dl in enumerate(dlist):
#         axx=plt.subplot(gs[2,ii])
#         dmin=datetime.datetime.strptime(dl, '%Y-%m-%d')-datetime.timedelta(days=5)
#         dmax=datetime.datetime.strptime(dl, '%Y-%m-%d')+datetime.timedelta(days=5)
#         # axx.contour(grden.distance[:-1],grden.depth,grden.pden.sel(date=slice(dmin,dmax)).mean(dim='date')[:-1,:].T,levels=arange(27,29,0.02),colors='grey',cmap=cm.YlGnBu)
#         # axx.contour(gvel.distance,gvel.depth,gvel.u.sel(date=slice(dmin,dmax)).mean(dim='date').T,colors='k',levels=[-0.3,-0.2,-0.1])
#         axx.contour(grden.distance[2:-1],grden.depth,grden.pden.sel(date=slice(dmin,dmax)).mean(dim='date')[2:-1,:].T,colors=clist[ii],linewidths=3,levels=[d1,d2,d3])
#         # axx.plot(all_bnds.distance,all_bnds.sel(date=slice(dmin,dmax)).mean(dim='date').T,color=clist[ii],linewidth=3)
#
#         axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
#         axx.axvline(distvec[4],color=cf5col,linewidth=4)
#         axx.axvline(distvec[5],color=cf6col,linewidth=4)
#         # axx.axvline(distvec[6],color=cf7col,linewidth=4)
#         axx.axvline(distvec[7],color=m1col,linewidth=4)
#         axx.axvline(oom_dist,color=ooicol,linewidth=4)
#         if ii!=0:
#             gca().set_yticklabels('')
#         else:
#             ylabel('depth [m]',fontsize=12)
#             text(120,600,'uIIW',fontsize=16)
#             text(120,1000,'dIIW',fontsize=16)
#         axx.set_title(dtitlist[ii])
#         axx.set_ylim(1500,0)
#         axx.set_xlim(30,175)
#
#     for ii,dl in enumerate(dlist):
#         axx=plt.subplot(gs[3,ii])
#         dmin=datetime.datetime.strptime(dl, '%Y-%m-%d')-datetime.timedelta(days=5)
#         dmax=datetime.datetime.strptime(dl, '%Y-%m-%d')+datetime.timedelta(days=5)
#         l1,l2,l3=TS_triptich(axx,dmin,dmax)
#         axx.set_ylim(3,5)
#         axx.set_yticks(arange(3,5,0.5))
#         axx.set_xticks(arange(34.86,35,0.04))
#         axx.set_xlim(34.87,34.985)
#         if ii!=0:
#             gca().set_yticklabels('')
#         else:
#             ylabel('pot. temperature [$^\circ$C]',fontsize=12)
#             axx.text(34.955,4.5,'uIIW',fontsize=14)
#             axx.text(34.955,4.08,'dIIW',fontsize=14)
#
#
#     f.text(0.5,0.29,'distance [km]', ha='center',fontsize=12)
#     f.text(0.5,0.08,'salinity', ha='center',fontsize=12)
#     savefig(figdir+'MixedLayer/paperfigs/ThicknessEvol_wTS.pdf',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/paperfigs/ThicknessEvol_wTS.png',bbox_inches='tight')
#
# ThickEvol_wTS()
#
#
# XXXXXXXXXXXXX
#
#
#
# # for mm in range(4,8):
# #     for ii in range(14):
# #         figure()
# #         plot(dat['potential density'][mm,:,ii*50:(ii+1)*50],dat.depth)
# #         title('CF'+str(mm+1)+': '+str(dat.date[ii*50].values)[:10]+' -- '+str(dat.date[(ii+1)*50].values)[:10])
# #         [axvline(dd,color='grey') for dd in arange(d1,d3,0.01)]
# #         [axvline(dd,color='k') for dd in [d1,d2,d3]]
# #         ylim(1500,0)
# #         xlim(27.4,27.85)
#
#
# ###############################################################################################################################
# ###############################################################################################################################
# ############################################## Mean vel and density contours   #################################################
# ###############################################################################################################################
# ###############################################################################################################################
# ######## This may move to another script with a map!
