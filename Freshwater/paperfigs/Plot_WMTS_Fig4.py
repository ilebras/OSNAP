from firstfuncs_1618 import *

figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'

# WM_nor=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc')
WM_obs=xr.open_dataset(datadir+'FW_WM/OSNAP2014-18_Tsub2020_WM_2008.nc')
# WM_mb=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_mb_2004.nc')

#copied over from CalcWMs...
sigmax_obs=27.56

salvec=linspace(31,36,103)
tmpvec=linspace(-3,40,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)
SAvec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CTvec=gsw.CT_from_pt(SAvec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
sigma1mat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SAvec[ii],CTvec[jj])

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange'}

def plot_TS_splitboth(WM_obs,namtit,xlims,ylims,xlim2,ylim2):
    f,[ax1,ax2]=subplots(1,2,figsize=(10,3.75))#,sharex=True,sharey=True)
    # ax2=ax1
    tscale=40
    for wm in ['AWS','PWS','DWS']:
        for ii,axx in enumerate([ax1,ax2]):
            if ii==0:
                axx.scatter(WM_obs['SA'].sel(TIME=slice('2015-1-1','2018-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM_obs['CT'].sel(TIME=slice('2015-1-1','2018-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                    s=abs(WM_obs['TRANS'].sel(TIME=slice('2015-1-1','2018-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values)*tscale,linewidth=3,alpha=0.8,label=wm,color=coldic[wm],zorder=100)
                axx.scatter(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                    s=abs(WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values)*tscale,linewidth=3,alpha=0.8,label='',color='k',zorder=100)
            else:
                axx.scatter(WM_obs['SA'].sel(TIME=slice('2015-1-1','2018-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM_obs['CT'].sel(TIME=slice('2015-1-1','2018-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                    s=abs(WM_obs['TRANS'].sel(TIME=slice('2015-1-1','2018-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values)*tscale*2,linewidth=3,alpha=0.8,label='',color=coldic[wm],zorder=100)
                axx.scatter(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                    s=abs(WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values)*tscale*2,linewidth=3,alpha=0.8,label='',color='k',zorder=100)
            if ('AWS' in wm):
                axx.plot(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
            else:
                axx.plot(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for wm in ['AWN','PWN']:
        ax1.scatter(WM_obs['SA'].sel(TIME=slice('2006-1-1','2011-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM_obs['CT'].sel(TIME=slice('2006-1-1','2011-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=abs(WM_obs['TRANS'].sel(TIME=slice('2006-1-1','2011-1-1')).groupby('TIME.year').mean('TIME').sel(WM=wm).values)*tscale,linewidth=3,label=wm,alpha=0.8,color=coldic[wm],zorder=100)
        ax1.scatter(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,
                        s=abs(WM_obs['TRANS'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values)*tscale,linewidth=3,label='',alpha=0.8,color='k',zorder=100)
        if ('PWN' in wm):
            ax1.plot(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM_obs['SA'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,WM_obs['CT'].groupby('TIME.month').mean('TIME').mean('month').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for axx in [ax1,ax2]:
        axx.contour(SAvec,CTvec,pdenmat,colors='grey',levels=arange(sigmax_obs-8,sigmax_obs+2,0.4),zorder=5,alpha=0.3)
        axx.contour(SAvec,CTvec,pdenmat,colors='k',levels=[sigmax_obs],zorder=500)
    ax1.set_ylabel('Conservative Temperature [$^\circ$C]',fontsize=14)
    f.text(0.5, 0, 'Absolute Salinity [g/kg]', ha='center',fontsize=14)
    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    ax2.set_xlim(xlim2)
    ax2.set_ylim(ylim2)
    # ax2.scatter(0,0,s=1*tscale,color='k',label='1 Sv')
    # ax2.scatter(0,0,s=5*tscale,color='k',label='5 Sv')
    lgnd=f.legend(loc=(0.9,0.25),fontsize=14)
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [150]
    # f.suptitle('Transport-weighted water mass properties\n\n',fontsize=16)
    ax2.set_title('Southern boundary',fontsize=16)
    ax1.set_title('All water masses',fontsize=16)
    ax1.plot([xlim2[0],xlim2[0],xlim2[1],xlim2[1],xlim2[0]],[ylim2[0],ylim2[1],ylim2[1],ylim2[0],ylim2[0]],'k-')
    for wm in ['AWN','PWN']:
        ax1.text(xlab[wm],ylab[wm],str(np.round(WM_obs.TRANS.sel(WM=wm).groupby('TIME.month').mean('TIME').mean('month').values,1))+' Sv',color=coldic[wm],fontsize=16,fontweight='bold')
    for wm in WM_obs.WM.values:
        ax1.text(xtit[wm],ytit[wm],wm,color=coldic[wm],fontsize=16,fontweight='bold')
    for wm in ['AWS','PWS','DWS']:
        ax2.text(xlab[wm],ylab[wm],str(np.round(WM_obs.TRANS.sel(WM=wm).groupby('TIME.month').mean('TIME').mean('month').values,1))+' Sv',color=coldic[wm],fontsize=16,fontweight='bold')
        # ax2.text(xtit[wm],ytit[wm],wm,color=coldic[wm],fontsize=16,fontweight='bold')
    savefig(figdir_paper+'TS_split_'+str(namtit)+'.png',bbox_inches='tight')
    savefig(figdir_paper+'TS_split_'+str(namtit)+'.pdf',bbox_inches='tight')


xlab={}
ylab={}
xlab['AWN']=34.9
ylab['AWN']=18
xlab['PWN']=33.6
ylab['PWN']=0.5
xlab['AWS']=35.22
ylab['AWS']=10.5
xlab['PWS']=34.45
ylab['PWS']=5.5
xlab['DWS']=35.19
ylab['DWS']=1.5

xtit={}
ytit={}
xtit['AWN']=35
ytit['AWN']=20.5
xtit['AWS']=34.85
ytit['AWS']=10
xtit['DWS']=35.3
ytit['DWS']=2.5
xtit['PWS']=34.45
ytit['PWS']=6
xtit['PWN']=33.6
ytit['PWN']=3

plot_TS_splitboth(WM_obs,'obs',[33.5,35.75],[-2,24],[34.4,35.7],[1,12])
