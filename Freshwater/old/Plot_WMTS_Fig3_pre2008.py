from firstfuncs_1618 import *

figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'

# WM_nor=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc')
WM_obs=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_2004.nc')
WM_obs_new=xr.open_dataset(datadir+'FW_WM/OSNAP2014-18_Tsub2020_WM_2008.nc')
# WM_mb=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_mb_2004.nc')


WM_obs.PTMP.sel(WM='AWN').plot()
WM_obs_new.CT.sel(WM='AWN').plot()
axhline(0)

#copied over from GetandPlotSec...
sigmax=27.6
sigmax_obs=27.54

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



coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}

WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month')

def plot_TS_splitboth(WM_obs,namtit,xlims,ylims,xlim2,ylim2):
    f,[ax1,ax2]=subplots(1,2,figsize=(10,3.75))
    tscale=100
    for wm in ['AWS','PWS','DWS']:
        # ax1.scatter(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').sel(WM=wm).values,
                # s=WM_obs['TRANS'].groupby('TIME.month').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax1.scatter(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                s=abs(WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values)*tscale,linewidth=3,alpha=0.8,label=wm,color=coldic[wm],zorder=100)
        if ('AWS' in wm):
            ax1.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for wm in ['AWN','PWN']:
        # ax2.scatter(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').sel(WM=wm).values,
                # s=WM_obs['TRANS'].groupby('TIME.month').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax2.scatter(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                s=abs(WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values)*tscale,linewidth=3,label=wm,alpha=0.8,color=coldic[wm],zorder=100)
        if ('PWN' in wm):
            ax2.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax2.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax2.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    ax2.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax1.set_ylabel('pot. temperature [$^\circ$C]',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    ax1.set_xlim(xlims)
    ax2.set_xlim(xlim2)
    ax1.set_ylim(ylims)
    ax2.set_ylim(ylim2)
    ax2.scatter(0,0,s=1*tscale,color='k',label='1 Sv')
    ax2.scatter(0,0,s=5*tscale,color='k',label='5 Sv')
    lgnd=f.legend(loc=(0.9,0.25),fontsize=14)
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [150]
    # f.suptitle('Transport-weighted water mass properties\n\n',fontsize=16)
    ax1.set_title('Southern boundary',fontsize=16)
    ax2.set_title('Northern boundary',fontsize=16)
    if 'obs' in namtit:
        ax2.set_yticklabels('')
    savefig(figdir_paper+'TS_split_'+str(namtit)+'.png',bbox_inches='tight')
    savefig(figdir_paper+'TS_split_'+str(namtit)+'.pdf',bbox_inches='tight')


plot_TS_splitboth(WM_obs,'obs',[33.5,35.5],[-1.5,10.5],[33.25,35.5],[-2,12])

# plot_TS_splitboth(WM_mb,'mb',[33.25,35.5],[-2,12],[33.25,35.5],[-2,17.5])
#
#
# plot_TS_splitboth(WM_nor,'nor',[33.5,35.8],[-2,12],[33.5,35.8],[-2,12])

############################################################################################
################  NORESM TS PLOT OF ANNUAL VARIATIONS  ############
############################################################################################
def plot_TS_yearly(WM_obs,namtit,xlims,ylims,xlim2,ylim2):
    f,[ax1,ax2]=subplots(1,2,figsize=(10,3.75))
    for wm in ['AWS','PWS','DWS']:
        ax1.scatter(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=WM_obs['TRANS'].groupby('TIME.year').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax1.scatter(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,
                s=WM_obs['TRANS'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values**2*4,linewidth=3,label='',color='k',zorder=100)
        if ('AWS' in wm):
            ax1.plot(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for wm in ['AWN','PWN']:
        ax2.scatter(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=WM_obs['TRANS'].groupby('TIME.year').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax2.scatter(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,
                s=WM_obs['TRANS'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values**2*4,linewidth=3,label='',color='k',zorder=100)
        if ('PWN' in wm):
            ax2.plot(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax2.plot(WM_obs['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)

    ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax2.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    ax2.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax1.set_ylabel('pot.temperature [$^\circ$C]',fontsize=14)
    f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
    ax1.set_xlim(xlims)
    ax2.set_xlim(xlim2)
    ax1.set_ylim(ylims)
    ax2.set_ylim(ylim2)
    lgnd=f.legend(loc='center right')
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [60]
    # f.suptitle('Transport-weighted water mass properties\n\n',fontsize=16)
    ax1.set_title('Southern boundary',fontsize=16)
    ax2.set_title('Northern boundary',fontsize=16)
    if 'obs' in namtit:
        ax2.set_yticklabels('')
    savefig(figdir_paper+'TS_split_yearly'+str(namtit)+'.png',bbox_inches='tight')
    savefig(figdir_paper+'TS_split_yearly'+str(namtit)+'.pdf',bbox_inches='tight')

plot_TS_yearly(WM_nor,'nor',[34.8,35.8],[4,12],[33.5,35.8],[-2,12])
############################################################################################
################ TS for PWS and DWS  ############
############################################################################################

def plot_TS_PWS(WM_obs,namtit,xlims,ylims):
    f,ax1=subplots(1,1,figsize=(5,3.75))
    for wm in ['AWS','PWS','PWN']:
        ax1.scatter(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                s=WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100)
        if ('AWS' in wm) | ('PWN' in wm):
            ax1.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)

    ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax1.set_ylabel('pot.temperature [$^\circ$C]',fontsize=14)
    ax1.set_xlabel('salinity',fontsize=14)
    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    lgnd=ax1.legend(loc=(1.05,0.3))
    for ii in range(3):
        lgnd.legendHandles[ii]._sizes = [60]
    savefig(figdir_paper+'TS_PWS.png',bbox_inches='tight')
    savefig(figdir_paper+'TS_PWS.pdf',bbox_inches='tight')


plot_TS_PWS(WM_obs,'obs',[33.25,35.5],[-2,12])

def plot_TS_DWS(WM_obs,namtit,xlims,ylims):
    f,ax1=subplots(1,1,figsize=(5,3.75))
    for wm in ['AWS','DWS']:
        ax1.scatter(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                s=WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100)
        if ('AWS' in wm) | ('PWN' in wm):
            ax1.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM_obs['PSAL'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['PTMP'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)

    ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
    ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax1.set_ylabel('pot.temperature [$^\circ$C]',fontsize=14)
    ax1.set_xlabel('salinity',fontsize=14)
    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    lgnd=ax1.legend(loc=(1.05,0.3))
    for ii in range(2):
        lgnd.legendHandles[ii]._sizes = [60]
    savefig(figdir_paper+'TS_DWS.png',bbox_inches='tight')
    savefig(figdir_paper+'TS_DWS.pdf',bbox_inches='tight')

plot_TS_DWS(WM_obs,'obs',[33.25,35.5],[-2,12])

############################################################################################
################  MAURITZEN TS PLOT + GRAVEYARD (see _pre2004 version for more GRAVEYARD)  ############
############################################################################################
#
# M96={}
# M96['WMS']=['AW_IS','AW_BS','AW_FS','AE_PW','AE_AAW','AW_IW','NE_RAW','ED_PW','ED_AAW','ED_RAW','AW_IR','AW_EI','DW']
# M96['PTMP']=[7,4,2,-1,0.5,0,1.5,-1,0.5,1,4,2.5,-1]
# M96['PSAL']=[35.3,35.05,35,34.3,34.85,34.89,34.95,34.3,34.85,34.92,34,34.9,34.9]
# M96['TRANS']=[6.8,-1.6,-3.4,1.5,1.5,0.17,1.1,-1.6,-2,-0.8,0.9,-0.7,-2]
#
# sign=array(M96['TRANS'])/array([abs(tt) for tt in M96['TRANS']])
#
# sign
# sum(M96['TRANS'])
#
# def plot_M96_TS():
#     f,axx=subplots(1,1,figsize=(6,5),sharex=True,sharey=True)
#     for ii in range(len(M96['WMS'])):
#         axx.scatter(M96['PSAL'][ii],M96['PTMP'][ii],s=M96['TRANS'][ii]**2*50,zorder=50,linewidth=3,label=M96['WMS'][ii])
#         if sign[ii]>0:
#             axx.plot(M96['PSAL'][ii],M96['PTMP'][ii],'k+',zorder=51)
#         else:
#             axx.plot(M96['PSAL'][ii],M96['PTMP'][ii],'k_',zorder=51)
#     axx.legend(loc=(1.05,0))
#     axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
#     xlim(34,36)
#     ylim(-2,12)
#     f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
#     f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
#
#     f.suptitle('Mauritzen 1996 (Pt.II)\nTransport-weighted water mass properties\n',fontsize=16)
#     savefig(figdir+'TS_Mauritzen.png',bbox_inches='tight')
#     savefig(figdir+'TS_Mauritzen.pdf',bbox_inches='tight')
#
# plot_M96_TS()


# def plot_TS_obs():
#     f,ax1=subplots(1,1,figsize=(5,4),sharex=True,sharey=True)
#     for wm in WM_obs.WM:
#         ax1.scatter(WM_obs['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,WM_obs['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,
#         s=WM_obs['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').values**2*4,zorder=50,linewidth=3,label=str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
#         ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,
#                 s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label='',color='k',zorder=100)
#         if ('AWS' in wm) | ('PWN' in wm):
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
#         else:
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
#     lgnd=ax1.legend(loc=(1.05,0.3))
#     for ii in range(5):
#         lgnd.legendHandles[ii]._sizes = [40]
#
#     ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
#     ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
#     f.text(0, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
#     f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
#     xlim(34,35.5)
#     ylim(-2,12)
#     f.suptitle('Transport-weighted water mass properties\n',fontsize=16)
#     # ax1.set_title('Observations',fontsize=14)
#     # ax2.set_title('NorESM',fontsize=14)
#     # savefig(figdir_paper+'TS_obs_only.png',bbox_inches='tight')
#     # savefig(figdir_paper+'TS_obs_only.pdf',bbox_inches='tight')
#
# plot_TS_obs()

# def plot_TS_obs_split1():
#     f,ax1=subplots(1,1,figsize=(5,4),sharex=True,sharey=True)
#     for wm in ['AWS','PWS','DWS']:
#         # ax1.scatter(WM_obs['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,WM_obs['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,
#         # s=WM_obs['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').values**2*4,zorder=50,linewidth=3,label=wm,color=coldic[wm],alpha=0.5)
#         ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,
#                 s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100)
#         if ('AWS' in wm) | ('PWN' in wm):
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
#         else:
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
#     lgnd=ax1.legend(loc=(1.05,0.3))
#     for ii in range(3):
#         lgnd.legendHandles[ii]._sizes = [40]
#     ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
#     ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
#     f.text(0, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
#     f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
#     xlim(34,35.5)
#     ylim(-2,12)
#     f.suptitle('Southern Boundary: OSNAP East\n',fontsize=16)
#     # ax1.set_title('Observations',fontsize=14)
#     # ax2.set_title('NorESM',fontsize=14)
#     savefig(figdir_paper+'TS_obs_OSNAP_only.png',bbox_inches='tight')
#     savefig(figdir_paper+'TS_obs_OSNAP_only.pdf',bbox_inches='tight')
#
# plot_TS_obs_split1()
#
# def plot_TS_obs_split2():
#     f,ax1=subplots(1,1,figsize=(5,4),sharex=True,sharey=True)
#     for wm in ['AWN','PWN']:
#         # ax1.scatter(WM_obs['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,WM_obs['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,
#         # s=WM_obs['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').values**2*4,zorder=50,linewidth=3,label=wm,color=coldic[wm],alpha=0.5)
#         ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,
#                 s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label='',color='k',zorder=100,label=wm,color=coldic[wm])
#         if ('AWS' in wm) | ('PWN' in wm):
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
#         else:
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
#     lgnd=ax1.legend(loc=(1.05,0.3))
#     for ii in range(2):
#         lgnd.legendHandles[ii]._sizes = [40]
#     ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
#     ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
#     f.text(0, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
#     f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
#     xlim(34,35.5)
#     ylim(-2,12)
#     f.suptitle('Northern Boundary\n',fontsize=16)
#     # ax1.set_title('Observations',fontsize=14)
#     # ax2.set_title('NorESM',fontsize=14)
#     savefig(figdir_paper+'TS_obs_FSBSO_only.png',bbox_inches='tight')
#     savefig(figdir_paper+'TS_obs_FSBSO_only.pdf',bbox_inches='tight')
#
# plot_TS_obs_split2()
#
#
# def plot_TS_bylayer():
#     f,[ax1,ax2]=subplots(1,2,figsize=(11,5),sharex=True,sharey=True)
#
#     for wm in WM_nor.WM:
#         ax2.scatter(WM_nor['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,WM_nor['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,
#         s=WM_nor['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME')**2*4,zorder=50,linewidth=3,label=str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
#         lgnd=ax2.legend(loc=(1.05,0.2))
#     for wm in WM_nor.WM:
#         ax2.scatter(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,
#                 s=WM_nor['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label='NorESM : '+str(wm.values),color='k',zorder=100)
#         if ('AWS' in wm) | ('PWN' in wm):
#             ax2.plot(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
#         else:
#             ax2.plot(WM_nor['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_nor['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
#     for wm in WM_obs.WM:
#         ax1.scatter(WM_obs['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,WM_obs['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').values,
#         s=WM_obs['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').values**2*4,zorder=50,linewidth=3,label='NorESM : '+str(wm.values),color=coldic[str(wm.values)],alpha=0.5)
#         ax1.scatter(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,
#                 s=WM_obs['TRANS'].sel(WM=wm).mean(dim='TIME').values**2*4,linewidth=3,label=str(wm.values),color='k',zorder=100)
#         if ('AWS' in wm) | ('PWN' in wm):
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'+',markersize=12,color='w',zorder=101)
#         else:
#             ax1.plot(WM_obs['PSAL'].sel(WM=wm).mean(dim='TIME').values,WM_obs['PTMP'].sel(WM=wm).mean(dim='TIME').values,'_',markersize=12,color='w',zorder=101)
#     for ii in range(5):
#         lgnd.legendHandles[ii]._sizes = [40]
#     for axx in [ax1,ax2]:
#         axx.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax-2,sigmax+2,0.2),zorder=5,alpha=0.5)
#         axx.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
#     f.text(0.05, 0.5, 'pot.temperature [$^\circ$C]', va='center', rotation='vertical',fontsize=14)
#     f.text(0.5, 0, 'salinity', ha='center',fontsize=14)
#     xlim(34,36)
#     ylim(-2,12)
#     f.suptitle('Transport-weighted water mass properties\n',fontsize=16)
#     ax1.set_title('Observations',fontsize=14)
#     ax2.set_title('NorESM',fontsize=14)
#     savefig(figdir_paper+'TS_comp.png',bbox_inches='tight')
#     savefig(figdir_paper+'TS_comp.pdf',bbox_inches='tight')
#
# plot_TS_bylayer()
