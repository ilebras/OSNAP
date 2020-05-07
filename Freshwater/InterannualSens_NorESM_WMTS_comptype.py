from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM_interannual/'

WM={}
WM['dpth']=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc')
WM['sig']=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004_sigma.nc')
WM['tcorr']=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004_sigma_transcorr.nc')
# WM['tandhcorr']=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004_sigma_transandheatcorr.nc')

WM['dpth']['PSAL'].sel(WM='AWN').plot()
WM['dpth']['PTMP'].sel(WM='AWN').plot()


#copied over from GetandPlotSec...
sigmax={}
sigmax['dpth']=27.6
sigmax['sig']=27.58
sigmax['tcorr']=27.58
# sigmax['tandhcorr']=27.58

startyear=2000
startmonth=1
endyear=2018
endmonth=12

def load_allyears(varlab):
    tmp1=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_200001-200912.nc')[0])
    tmp2=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_201001-201812.nc')[0])
    varout=xr.concat([tmp1,tmp2],dim='time')
    varout=varout.rename({'time':'TIME'})
    varout['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return varout

hf=load_allyears('heatloss')
so=xr.open_dataset(datadir+'NorESM/NorESM_source_storage_xray_18yrs_2004.nc')

salvec=linspace(31,36.5,103)
tmpvec=linspace(-3,20,103)
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


coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'cyan','Q':'limegreen'}

############################################################################################
################  NORESM TS PLOT OF ANNUAL VARIATIONS  ############
############################################################################################
def plot_TS_yearly(WM,namtit,xlims,ylims,xlim2,ylim2):
    f,[ax1,ax2]=subplots(1,2,figsize=(10,3.75))
    for wm in ['AWS','PWS','DWS']:
        ax1.scatter(WM['PSAL'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax1.scatter(WM['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values**2*4,linewidth=3,label='',color='k',zorder=100)
        if ('AWS' in wm):
            ax1.plot(WM['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for wm in ['AWN','PWN']:
        ax2.scatter(WM['PSAL'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax2.scatter(WM['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values**2*4,linewidth=3,label='',color='k',zorder=100)
        if ('PWN' in wm):
            ax2.plot(WM['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax2.plot(WM['PSAL'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['PTMP'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)

    ax1.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax[namtit]-6,sigmax[namtit]+2,0.2),zorder=5,alpha=0.5)
    ax1.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax[namtit]],zorder=500)
    ax2.contour(salvec,tmpvec,pdenmat,colors='grey',levels=arange(sigmax[namtit]-6,sigmax[namtit]+2,0.2),zorder=5,alpha=0.5)
    ax2.contour(salvec,tmpvec,pdenmat,colors='k',levels=[sigmax[namtit]],zorder=500)
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
    savefig(figdir+'TS_split_yearly_'+str(namtit)+'.png',bbox_inches='tight')
    savefig(figdir+'TS_split_yearly_'+str(namtit)+'.pdf',bbox_inches='tight')

typevec=['dpth','sig','tcorr']#,'tandhcorr']

# future me: maybe add a box in right panel showing range of left panel
for kk in typevec:
    plot_TS_yearly(WM[kk],kk,[34.8,35.8],[2,12],[33.5,36.3],[-3,20])

############################################################################################
################### Solve inverse model with NorESM inputs #################################
############################################################################################


# First, solve for 19 year average
# Then, solve for each year
# Finally, solve for each  combination of years (FOR NORTH AND SOUTH ONLY)

#######################################################################
################### Full 19 year mean #################################
#######################################################################

cp=3850
rhow=1000

Snorm=35
Tnorm=5
Winv=diag([1,1/Snorm,1/Tnorm])

def get_U_S_T_from_WM(kk):
    U={}
    S={}
    T={}
    for wm in WM[kk].WM:
        U[str(wm.values)]=float(WM[kk]['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        S[str(wm.values)]=float(WM[kk]['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        T[str(wm.values)]=float(WM[kk]['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)

    U['FW']=so['FW+SI'].groupby('TIME.month').mean('TIME').mean(dim='month').values
    U['Q']=hf['NORDIC_hflx'].groupby('TIME.month').mean('TIME').mean(dim='month').values*1e6/rhow/cp
    S['FW']=0
    T['FW']=0
    T['Q']=1

    return U,S,T


def get_U_from_x(x):
    U={}
    U['PWS']=x[0]
    U['AWS']=x[1]
    U['DWS']=x[2]
    U['PWN']=x[3]
    U['AWN']=x[4]
    U['FW']=x[5]
    U['Q']=x[6]
    return U

def plot_base_case_simple(kk):
    f,axx=subplots(1,4,figsize=(12,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,1,1]))
    alf=0.85
    capi=7
    #U
    axx[0].bar(range(2),[Ubase[kk] for kk in ['AWS','DWS']],color=[coldic[kk] for kk in ['AWS','DWS']],yerr=[Ue[kk] for kk in ['AWS','DWS']],capsize=capi,alpha=alf)
    axx[0].plot(range(2),[U[kk] for kk in ['AWS','DWS']],'o',color='k')

    ylimi=20
    axx[0].set_ylim(-ylimi,ylimi)
    ylimi=4.25
    axx[1].set_ylim(-ylimi,ylimi)
    axx[1].bar(range(3),[Ubase[kk] for kk in ['PWS','PWN','AWN']],color=[coldic[kk] for kk in ['PWS','PWN','AWN']],yerr=[Ue[kk] for kk in ['PWS','PWN','AWN']],capsize=capi,alpha=alf)
    axx[1].plot(range(3),[U[kk] for kk in ['PWS','PWN','AWN']],'o',color='k')

    axx[2].bar(0,[Ubase[kk] for kk in ['FW']],color=[coldic[kk] for kk in ['FW']],yerr=[Ue[kk] for kk in ['FW']],capsize=capi,alpha=alf)
    axx[2].plot(0,[U[kk] for kk in ['FW']],'o',color='k')
    fwlim=0.2
    axx[2].set_ylim(-fwlim,fwlim)

    fsz=14
    axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    axx[3].set_ylabel('Heat flux [TW]',fontsize=fsz)
    axx[3].bar(0,cp*rhow*(Ubase['Q'])/1e6,color=coldic['Q'],yerr=cp*rhow*Ue['Q']/1e6,capsize=capi,alpha=alf)
    axx[3].plot(0,cp*rhow*(U['Q'])/1e6,'o',color='k')

    for ii in range(3):
        axx[ii].axhline(0,color='k')
    axx[0].set_xticks(range(2))
    axx[0].set_xticklabels(['AWS','DWS'])
    axx[1].set_xticks(range(3))
    axx[1].set_xticklabels(['PWS','PWN','AWN'])
    axx[2].set_xticks(range(1))
    axx[2].set_xticklabels(['FW'])
    axx[3].set_xticks([0])
    axx[3].set_xticklabels('Q')

    savefig(figdir+'InvBudSol_NorESM_fullmean_'+kk+'.png',bbox_inches='tight')
    savefig(figdir+'InvBudSol_NorESM_fullmean_'+kk+'.pdf',bbox_inches='tight')


for kk in typevec:
    U,S,T=get_U_S_T_from_WM(kk)

    AM=array([[1,1,1,1,1,1,0],\
    [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],0],\
    [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],1]])

    x0=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['Q']]

    dv=-AM.dot(x0)

    Qvar=10
    Evec=hstack((abs(array([U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN']])/10),0.02,Qvar))
    # Evec=hstack((5*[1],0.02,Qvar))
    E=diag(Evec)

    Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

    Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
    Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
    xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
    xsol_Ad=E.dot(xsol_prime)
    xbase=x0+xsol_Ad
    Ubase=get_U_from_x(xbase)
    Ue=get_U_from_x(Evec)
    plot_base_case_simple(kk)

Uinit_mean=U
Ubase_mean=Ubase

#######################################################################
################### Each year individually ############################
#######################################################################

def get_U_S_T_from_WM_eachyear(kk):
    U={}
    S={}
    T={}
    for wm in WM[kk].WM:
        U[str(wm.values)]=WM[kk]['TRANS'].sel(WM=wm).groupby('TIME.year').mean(dim='TIME').values
        S[str(wm.values)]=WM[kk]['PSAL'].sel(WM=wm).groupby('TIME.year').mean(dim='TIME').values
        T[str(wm.values)]=WM[kk]['PTMP'].sel(WM=wm).groupby('TIME.year').mean(dim='TIME').values

    U['FW']=so['FW+SI'].groupby('TIME.year').mean(dim='TIME')
    U['Q']=hf['NORDIC_hflx'].groupby('TIME.year').mean(dim='TIME')*1e6/rhow/cp
    S['FW']=zeros(len(U['FW']))
    T['FW']=zeros(len(U['FW']))
    T['Q']=ones(len(U['FW']))

    return U,S,T

kk='tcorr'
U,S,T=get_U_S_T_from_WM_eachyear(kk)

Qvar=10
Evec=hstack((abs(array([mean(U['PWS']),mean(U['AWS']),mean(U['DWS']),mean(U['PWN']),mean(U['AWN'])])/10),0.02,Qvar))
E=diag(Evec)
yrs=19

Usol={}
for ii in range(yrs):
        AM=array([[1,1,1,1,1,1,0],\
        [S['PWS'][ii],S['AWS'][ii],S['DWS'][ii],S['PWN'][ii],S['AWN'][ii],S['FW'][ii],0],\
        [T['PWS'][ii],T['AWS'][ii],T['DWS'][ii],T['PWN'][ii],T['AWN'][ii],T['FW'][ii],1]])

        x0=[U['PWS'][ii],U['AWS'][ii],U['DWS'][ii],U['PWN'][ii],U['AWN'][ii],U['FW'][ii],U['Q'][ii]]

        dv=-AM.dot(x0)

        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        xbase=x0+xsol_Ad
        Usol[ii]=get_U_from_x(xbase)

#######################################################################
################### Now try combining years ############################
#######################################################################

Ucomb={}
for ii in range(yrs):
    Ucomb[ii]={}
    for jj in range(yrs):
        if ii==jj:
            Ucomb[ii][jj]=get_U_from_x(NaN*ones(7))
        else:
            AM=array([[1,1,1,1,1,1,0],\
            [S['PWS'][ii],S['AWS'][ii],S['DWS'][ii],S['PWN'][jj],S['AWN'][jj],0,0],\
            [T['PWS'][ii],T['AWS'][ii],T['DWS'][ii],T['PWN'][jj],T['AWN'][jj],0,1]])

            x0=[U['PWS'][ii],U['AWS'][ii],U['DWS'][ii],U['PWN'][jj],U['AWN'][jj],mean(U['FW']),mean(U['Q'])]

            dv=-AM.dot(x0)

            Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

            Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
            Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
            xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
            xsol_Ad=E.dot(xsol_prime)
            xbase=x0+xsol_Ad
            Ucomb[ii][jj]=get_U_from_x(xbase)




def plot_erranom(kk):
    f,axx=subplots(1,1,figsize=(12,3),constrained_layout=True)
    alf=0.85
    capi=7
    xshi1=0.8
    xshi2=1
    xshi3=1.2
    xwi=0.2
    Ucomb_all={}
    for ii in range(yrs):
        for jj in range(yrs):
            for ll in ['AWS','DWS','PWS','AWN','PWN','FW','Q']:
                if ii!=jj:
                    if ii==0:
                        if jj==1:
                            Ucomb_all[ll]=Ucomb[ii][jj][ll]
                    else:
                        Ucomb_all[ll]=hstack((Ucomb_all[ll], Ucomb[ii][jj][ll]))
    boxo1=axx.boxplot([(Ucomb_all[ll]-Ubase_mean[ll])/Ubase_mean[ll] for ll in ['AWS','DWS','PWS','AWN','PWN','FW','Q']],patch_artist=True,notch=True,positions=(arange(7)+xshi3), showfliers=False,widths=xwi)
    boxi1=axx.boxplot([[(Usol[kk][ll]-Ubase_mean[ll])/Ubase_mean[ll] for kk in Usol] for ll in ['AWS','DWS','PWS','AWN','PWN','FW','Q']],patch_artist=True,notch=True,positions=(arange(7)+xshi2), showfliers=False,widths=xwi)
    box1=axx.boxplot([[(Usol[kk][ll]-U[ll][kk])/Ubase_mean[ll] for kk in Usol] for ll in ['AWS','DWS','PWS','AWN','PWN','FW','Q']],patch_artist=True,notch=True,positions=(arange(7)+xshi1), showfliers=False,widths=xwi)
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
            plt.setp(boxo1[item], color='#02818a')
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
            plt.setp(boxi1[item], color='magenta')
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
            plt.setp(box1[item], color='k')
    axhline(0,color='k')
    axx.set_xticks(range(1,8))
    axx.set_xticklabels(['AWS','DWS','PWS','AWN','PWN','FW','Q'])
    axx.set_ylim(-1,1)
    axx.set_ylabel('Error ratio')
    savefig(figdir+'InvBudSol_NorESM_ErrAnom.png',bbox_inches='tight')
    savefig(figdir+'InvBudSol_NorESM_ErrAnom.pdf',bbox_inches='tight')

plot_erranom(kk)



# def plot_base_case_eachyear(kk):
#     f,axx=subplots(1,4,figsize=(12,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,1,1]))
#     alf=0.85
#     capi=7
#     xshi1=0.8
#     xshi2=1
#     xshi3=1.2
#     xwi=0.2
#     boxo1=axx[0].boxplot([array([[Ucomb[ii][jj][ll] for ii in range(yrs)] for jj in range(yrs)]).flatten() for ll in ['AWS','DWS']],patch_artist=True,notch=True,positions=(arange(2)+xshi3), showfliers=False,widths=xwi)
#     boxi1=axx[0].boxplot([[Usol[kk][ll] for kk in Usol] for ll in ['AWS','DWS']],patch_artist=True,notch=True,positions=(arange(2)+xshi2), showfliers=False,widths=xwi)
#     box1=axx[0].boxplot([U[kk] for kk in ['AWS','DWS']],patch_artist=True,notch=True,positions=(arange(2)+xshi1), showfliers=False,widths=xwi)
#     axx[0].plot(arange(2)+xshi1,[Uinit_mean[ll] for ll in ['AWS','DWS']],'wo',zorder=100,mec='k')
#     axx[0].plot(arange(2)+xshi2,[Ubase_mean[ll] for ll in ['AWS','DWS']],'wo',zorder=100,mec='k')
#     boxo2=axx[1].boxplot([array([[Ucomb[ii][jj][ll] for ii in range(yrs)] for jj in range(yrs)]).flatten() for ll in ['PWS','PWN','AWN']],patch_artist=True,notch=True,positions=(arange(3)+xshi3), showfliers=False,widths=xwi)
#     boxi2=axx[1].boxplot([[Usol[kk][ll] for kk in Usol] for ll in ['PWS','PWN','AWN']],patch_artist=True,notch=True,positions=(arange(3)+xshi2), showfliers=False,widths=xwi)
#     box2=axx[1].boxplot([U[kk] for kk in ['PWS','PWN','AWN']],patch_artist=True,notch=True,positions=(arange(3)+xshi1), showfliers=False,widths=xwi)
#     axx[1].plot(arange(3)+xshi1,[Uinit_mean[ll] for ll in ['PWS','PWN','AWN']],'wo',zorder=100,mec='k')
#     axx[1].plot(arange(3)+xshi2,[Ubase_mean[ll] for ll in ['PWS','PWN','AWN']],'wo',zorder=100,mec='k')
#     boxo3=axx[2].boxplot(array([[Ucomb[ii][jj]['FW'] for ii in range(yrs)] for jj in range(yrs)]).flatten(),patch_artist=True,notch=True,positions=(arange(1)+xshi3), showfliers=False,widths=xwi)
#     boxi3=axx[2].boxplot([Usol[kk]['FW'] for kk in Usol],patch_artist=True,notch=True,positions=(arange(1)+xshi2), showfliers=False,widths=xwi)
#     box3=axx[2].boxplot(U['FW'].values,patch_artist=True,notch=True,positions=(arange(1)+xshi1), showfliers=False,widths=xwi)
#     axx[2].plot(arange(1)+xshi1,Uinit_mean['FW'],'wo',zorder=100,mec='k')
#     axx[2].plot(arange(1)+xshi2,Ubase_mean['FW'],'wo',zorder=100,mec='k')
#     boxo4=axx[3].boxplot(array([[cp*rhow*Ucomb[ii][jj]['Q']/1e6 for ii in range(yrs)] for jj in range(yrs)]).flatten(),patch_artist=True,notch=True,positions=(arange(1)+xshi3), showfliers=False,widths=xwi)
#     boxi4=axx[3].boxplot([cp*rhow*Usol[kk]['Q']/1e6 for kk in Usol],patch_artist=True,notch=True,positions=(arange(1)+xshi2), showfliers=False,widths=xwi)
#     box4=axx[3].boxplot([cp*rhow*(U['Q'])/1e6],patch_artist=True,notch=True,positions=(arange(1)+xshi1), showfliers=False,widths=xwi)
#     axx[3].plot(arange(1)+xshi1,cp*rhow*Uinit_mean['Q']/1e6,'wo',zorder=100,mec='k')
#     axx[3].plot(arange(1)+xshi2,cp*rhow*Ubase_mean['Q']/1e6,'wo',zorder=100,mec='k')
#     ylimi=20
#     axx[0].set_ylim(-ylimi,ylimi)
#     ylimi=3.5
#     axx[1].set_ylim(-ylimi,ylimi)
#     fwlim=0.2
#     axx[2].set_ylim(-fwlim,fwlim)
#     fsz=14
#     axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
#     axx[3].set_ylabel('Heat flux [TW]',fontsize=fsz)
#     for boxi in [boxo1,boxo2,boxo3,boxo4]:
#         for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
#                 plt.setp(boxi[item], color='#02818a')
#     for boxi in [boxi1,boxi2,boxi3,boxi4]:
#         for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
#                 plt.setp(boxi[item], color='magenta')
#     for boxi in [box1,box2,box3,box4]:
#         for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
#                 plt.setp(boxi[item], color='k')
#     for ii in range(3):
#         axx[ii].axhline(0,color='k')
#     axx[0].set_xticks(range(1,3))
#     axx[0].set_xticklabels(['AWS','DWS'])
#     axx[1].set_xticks(range(1,4))
#     axx[1].set_xticklabels(['PWS','PWN','AWN'])
#     axx[2].set_xticks(range(1,2))
#     axx[2].set_xticklabels(['FW'])
#     axx[3].set_xticks([1])
#     axx[3].set_xticklabels('Q')
#
#     savefig(figdir+'InvBudSol_NorESM_eachyear.png',bbox_inches='tight')
#     savefig(+'InvBudSol_NorESM_eachyear.pdf',bbox_inches='tight')
#
# plot_base_case_eachyear(kk)



epsilon=arange(0,1.1,0.1)
a_pwmat=zeros((len(epsilon),yrs,yrs))
a_mean=zeros(len(epsilon))
for ii,ee in enumerate(1-epsilon):
    a_mean[ii]=(ee*Ubase_mean['PWN']*(Smean['PWN']/Smean['AWS']-1)+Ubase_mean['PWS']*(Smean['PWS']/Smean['AWS']-1))/Ubase_mean['FW']
    for jj in range(yrs):
        for kk in range(yrs):
                a_pwmat[ii,jj,kk]=(ee*Ucomb[jj][kk]['PWN']*(S['PWN'][kk]/S['AWS'][jj]-1)+Ucomb[jj][kk]['PWS']*(S['PWS'][jj]/S['AWS'][jj]-1))/Ucomb[jj][kk]['FW']


def plot_epsilon_dep():
    [plot(1-epsilon,a_pwmat[:,:,ii],'r',alpha=0.2) for ii in range(19)];
    plot(1-epsilon,a_mean,'k',linewidth=3)
    axhline(0,color='k')
    xlim(0,1)

plot_epsilon_dep()
