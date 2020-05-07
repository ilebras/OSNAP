from firstfuncs_1618 import *
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'


WM=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004_sigma_transcorr.nc')

#copied over from GetandPlotSec...
sigmax=27.58

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


coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'k','Q':'limegreen'}

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

plot_TS_yearly(WM,'nor',[34.8,35.8],[4,12.25],[33.5,36.2],[-2.5,17])

############################################################################################
################### Solve inverse model with NorESM inputs #################################
############################################################################################

# First, solve for 19 year average
# Then, solve for each year
# Finally, solve for each combination of years

cp=3850
rhow=1000
tera=10**12
#Noresm (taking sea ice into account)
Q=-252*tera/rhow/cp/1e6 #for the Sverdrups

def get_U_S_T_from_WM():
    U={}
    S={}
    T={}
    for wm in WM.WM:
        U[str(wm.values)]=float(WM['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        S[str(wm.values)]=float(WM['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        T[str(wm.values)]=float(WM['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)

    U['SI']=0.073 # NorESM fresh water input v. similar to Kwok et al. 2004 70mSv
    U['FW']=0.028 # mean E-P from JRA55
    U['Q']=Q
    S['SI']=0
    S['FW']=0
    T['SI']=0
    T['FW']=0
    T['Q']=1

    return U,S,T

U,S,T=get_U_S_T_from_WM()
#Save out the initial conditions
Ug=U.copy()
Sg=S.copy()
Tg=T.copy()

Ug


def get_U_from_x(x):
    U={}
    U['PWS']=x[0]
    U['AWS']=x[1]
    U['DWS']=x[2]
    U['PWN']=x[3]
    U['AWN']=x[4]
    U['FW']=x[5]
    U['SI']=x[6]
    U['Q']=x[7]
    return U

AM={}

AM['base']=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

x0=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]


zz='base'
dv=-AM[zz].dot(x0)

Snorm=35
Tnorm=5
Winv=diag([1,1/Snorm,1/Tnorm])

Qvar=10
Evec=hstack((abs(array([U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN']])/5),0.02,0.02,Qvar))
E=diag(Evec)

Umat,D,VmatT=linalg.svd(Winv.dot(AM[zz].dot(E)))

Lambda_inv = zeros((AM[zz].shape[0], AM[zz].shape[1])).T
Lambda_inv[:AM[zz].shape[0], :AM[zz].shape[0]] = diag(1/D)
xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
xsol_Ad=E.dot(xsol_prime)
xbase=x0+xsol_Ad
P=diag(E-E.dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E.dot(AM[zz].T))+linalg.inv(Winv)).dot(AM[zz].dot(E)))))
Ubase=get_U_from_x(xbase)
Ue=get_U_from_x(P)

Ubase


def plot_base_case_simple():
    f,axx=subplots(1,4,figsize=(13,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,2,1]))
    alf=0.75
    capi=7
    #U
    axx[0].bar(range(2),[Ubase[kk] for kk in ['AWS','DWS']],color=[coldic[kk] for kk in ['AWS','DWS']],yerr=[Ue[kk] for kk in ['AWS','DWS']],capsize=capi,alpha=alf)
    axx[0].plot(range(2),[Ug[kk] for kk in ['AWS','DWS']],'o',color='k')

    ylimi=20
    axx[0].set_ylim(-ylimi,ylimi)
    ylimi=4.25
    axx[1].set_ylim(-ylimi,ylimi)
    axx[1].bar(range(3),[Ubase[kk] for kk in ['PWS','PWN','AWN']],color=[coldic[kk] for kk in ['PWS','PWN','AWN']],yerr=[Ue[kk] for kk in ['PWS','PWN','AWN']],capsize=capi,alpha=alf)
    axx[1].plot(range(3),[Ug[kk] for kk in ['PWS','PWN','AWN']],'o',color='k')

    axx[2].bar(range(2),[Ubase[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']],yerr=[Ue[kk] for kk in ['SI','FW']],capsize=capi,alpha=alf)
    axx[2].plot(range(2),[Ug[kk] for kk in ['SI','FW']],'o',color='k')
    fwlim=0.1
    axx[2].set_ylim(-fwlim,fwlim)

    fsz=14
    axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    axx[3].set_ylabel('Heat flux [TW]',fontsize=fsz)
    axx[3].bar(0,cp*rhow*(Ubase['Q'])/1e6,color=coldic['Q'],yerr=cp*rhow*Ue['Q']/1e6,capsize=capi,alpha=alf)
    axx[3].plot(0,cp*rhow*(Ug['Q'])/1e6,'o',color='k')

    for ii in range(3):
        axx[ii].axhline(0,color='k')
        axx[ii].set_xticks(range(2))
    axx[0].set_xticklabels(['AWS','DWS'])
    axx[1].set_xticks(range(3))
    axx[1].set_xticklabels(['PWS','PWN','AWN'])
    axx[2].set_xticklabels(['SI','FW'])
    axx[3].set_xticks([0])
    axx[3].set_xticklabels('Q')

    savefig(figdir_paper+'InvBudSol_NorESM_fullmean_sigma_transcorr.png',bbox_inches='tight')
    savefig(figdir_paper+'InvBudSol_NorESM_fullmean_sigma_transcorr.pdf',bbox_inches='tight')

WM

WM.TRANS.sum(dim='WM').mean(dim='TIME')

plot_base_case_simple()

basediff=[(kk,Ubase[kk]-U[kk]) for kk in Ubase]

basediff
