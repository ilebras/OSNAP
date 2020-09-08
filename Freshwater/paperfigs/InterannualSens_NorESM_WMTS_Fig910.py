from firstfuncs_1618 import *
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs'

WM=xr.open_dataset(datadir+'NorESM/NorESM_WMs_18yrs_2004.nc')

#copied over from CalcWMs...
sigmax=27.59

salvec=linspace(31,36,103)
tmpvec=linspace(-8,28,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)
SAvec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CTvec=gsw.CT_from_pt(SAvec,tmpvec)
pdenmat=zeros((shape(salmat)))
pdenmat2=zeros((shape(salmat)))
sigma1mat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SAvec[ii],CTvec[jj])

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'cyan','Q':'limegreen'}

############################################################################################
################  NORESM TS PLOT OF ANNUAL VARIATIONS  ############
############################################################################################
def plot_TS_yearly(WM,namtit,xlims,ylims,xlim2,ylim2):
    f,[ax2,ax1]=subplots(1,2,figsize=(10,3.75))
    for wm in ['AWS','PWS','DWS']:
        ax1.scatter(WM['SA'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax1.scatter(WM['SA'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values**2*4,linewidth=3,label='',color='k',zorder=100)
        if ('AWS' in wm):
            ax1.plot(WM['SA'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax1.plot(WM['SA'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for wm in ['AWN','PWN']:
        ax2.scatter(WM['SA'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').sel(WM=wm).values**2*4,linewidth=3,label=wm,color=coldic[wm],zorder=100,alpha=0.8)
        ax2.scatter(WM['SA'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,
                s=WM['TRANS'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values**2*4,linewidth=3,label='',color='k',zorder=100)
        if ('PWN' in wm):
            ax2.plot(WM['SA'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'+',markersize=12,color='w',zorder=101)
        else:
            ax2.plot(WM['SA'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,WM['CT'].groupby('TIME.year').mean('TIME').mean('year').sel(WM=wm).values,'_',markersize=12,color='w',zorder=101)
    for axi in [ax1,ax2]:
        axi.contour(SAvec,CTvec,pdenmat,colors='grey',levels=arange(sigmax-4,sigmax+2,0.4),zorder=5,alpha=0.5)
        axi.contour(SAvec,CTvec,pdenmat,colors='k',levels=[sigmax],zorder=500)
    ax2.set_ylabel('Conservative Temperature [$^\circ$C]',fontsize=14)
    f.text(0.5, 0, 'Absolute Salinity [g/kg]', ha='center',fontsize=14)
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
    ax2.plot([xlims[0],xlims[0],xlims[1],xlims[1],xlims[0]],[ylims[0],ylims[1],ylims[1],ylims[0],ylims[0]],'k-')
    if 'obs' in namtit:
        ax2.set_yticklabels('')
    for wm in ['AWN','PWN']:
        ax2.text(xlab[wm],ylab[wm],str(np.round(WM.TRANS.sel(WM=wm).groupby('TIME.month').mean('TIME').mean('month').values,1))+' Sv',color=coldic[wm],fontsize=16,fontweight='bold')
    for wm in ['AWS','PWS','DWS']:
        ax1.text(xlab[wm],ylab[wm],str(np.round(WM.TRANS.sel(WM=wm).groupby('TIME.month').mean('TIME').mean('month').values,1))+' Sv',color=coldic[wm],fontsize=16,fontweight='bold')
    savefig(figdir_paper+'/TS_split_yearly'+str(namtit)+'.png',bbox_inches='tight')
    savefig(figdir_paper+'/TS_split_yearly'+str(namtit)+'.pdf',bbox_inches='tight')

xlab={}
ylab={}
xlab['AWN']=35.25
ylab['AWN']=13
xlab['PWN']=34
ylab['PWN']=0
xlab['AWS']=35.35
ylab['AWS']=11
xlab['PWS']=35.02
ylab['PWS']=8
xlab['DWS']=35.5
ylab['DWS']=4

plot_TS_yearly(WM,'nor',[35,36],[3,12],[33.5,36.2],[-3,19])



############################################################################################
################### Solve inverse model with NorESM inputs #################################
############################################################################################

# First, solve for 19 year average
# Then, solve for each year
# Finally, solve for each  combination of years (FOR NORTH AND SOUTH ONLY)

#######################################################################
################### Full 19 year mean #################################
#######################################################################

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

cp=3850
rhow=1025

Snorm=35
Tnorm=5
Winv=diag([1,1/Snorm,1/Tnorm])

def get_U_S_T_from_WM():
    U={}
    U_std={}
    S={}
    T={}
    for wm in WM.WM:
        U[str(wm.values)]=float(WM['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        U_std[str(wm.values)]=float(WM['TRANS'].sel(WM=wm).std('TIME').values)
        S[str(wm.values)]=float(WM['SA'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        T[str(wm.values)]=float(WM['CT'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)

    U['FW']=so['FW+SI'].groupby('TIME.month').mean('TIME').mean(dim='month').values
    U['Q']=hf['NORDIC_hflx'].groupby('TIME.month').mean('TIME').mean(dim='month').values*1e6/rhow/cp
    U_std['FW']=so['FW+SI'].std(dim='TIME').values
    U_std['Q']=hf['NORDIC_hflx'].std(dim='TIME').values*1e6/rhow/cp #1e6 for the Sverdrups
    S['FW']=0
    S['Q']=0
    T['FW']=0
    T['Q']=1

    return U,U_std,S,T


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

def plot_base_case_simple():
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

    savefig(figdir_paper+'_extra_2004/InvBudSol_NorESM_fullmean.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/InvBudSol_NorESM_fullmean.pdf',bbox_inches='tight')


kk='UVcorr'
U,U_std,S,T=get_U_S_T_from_WM()

AM=array([[1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],1]])

x0=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['Q']]

dv=-AM.dot(x0)

wmvec=['PWS','AWS','DWS','PWN','AWN','FW','Q']
Evec=[abs(U_std[wm]/sqrt(19*4)) for wm in wmvec]
E=diag(Evec)

Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
xsol_Ad=E.dot(xsol_prime)
xbase=x0+xsol_Ad
Ubase=get_U_from_x(xbase)
Ue=get_U_from_x(Evec)
plot_base_case_simple()

Uinit_mean

Uinit_mean=U
Ubase_mean=Ubase
Smean=S.copy()
Tmean=T.copy()

Udiff_mean={}
for wm in Ubase_mean:
    Udiff_mean[wm]=Ubase_mean[wm]-Uinit_mean[wm]

Udiff_mean['FW']*1e3

#######################################################################
################### Each year individually ############################
#######################################################################

def get_U_S_T_from_WM_eachyear():
    U={}
    U_std={}
    S={}
    T={}
    for wm in WM.WM:
        U[str(wm.values)]=WM['TRANS'].sel(WM=wm).groupby('TIME.year').mean(dim='TIME').values
        U_std[str(wm.values)]=WM['TRANS'].sel(WM=wm).groupby('TIME.year').std(dim='TIME').values
        S[str(wm.values)]=WM['SA'].sel(WM=wm).groupby('TIME.year').mean(dim='TIME').values
        T[str(wm.values)]=WM['CT'].sel(WM=wm).groupby('TIME.year').mean(dim='TIME').values

    U['FW']=so['FW+SI'].groupby('TIME.year').mean(dim='TIME')
    U['Q']=hf['NORDIC_hflx'].groupby('TIME.year').mean(dim='TIME')*1e6/rhow/cp
    U_std['FW']=so['FW+SI'].groupby('TIME.year').std(dim='TIME')
    U_std['Q']=hf['NORDIC_hflx'].groupby('TIME.year').std(dim='TIME')*1e6/rhow/cp
    S['FW']=zeros(len(U['FW']))
    T['FW']=zeros(len(U['FW']))
    T['Q']=ones(len(U['FW']))
    S['Q']=zeros(len(U['FW']))

    return U,S,T,U_std

U,S,T,U_std=get_U_S_T_from_WM_eachyear()

yrs=19

Usol={}
for ii in range(yrs):
        AM=array([[1,1,1,1,1,1,0],\
        [S['PWS'][ii],S['AWS'][ii],S['DWS'][ii],S['PWN'][ii],S['AWN'][ii],S['FW'][ii],0],\
        [T['PWS'][ii],T['AWS'][ii],T['DWS'][ii],T['PWN'][ii],T['AWN'][ii],T['FW'][ii],1]])

        x0=[U['PWS'][ii],U['AWS'][ii],U['DWS'][ii],U['PWN'][ii],U['AWN'][ii],U['FW'][ii],U['Q'][ii]]

        dv=-AM.dot(x0)
        Evec=[abs(U_std[wm][ii]/sqrt(4)) for wm in wmvec]
        E=diag(Evec)

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
if 'AWN' in ['AWN','PWN']:
    print('w')
T

Ucomb={}
Tcomb={}
Scomb={}
for ii in range(yrs):
    Ucomb[ii]={}
    Tcomb[ii]={}
    Scomb[ii]={}
    for jj in range(yrs):
        Tcomb[ii][jj]={}
        Scomb[ii][jj]={}
        for wm in S:
            if wm in ['AWS','DWS','PWS']:
                Tcomb[ii][jj][wm]=T[wm][ii]
                Scomb[ii][jj][wm]=S[wm][ii]
            elif wm in ['AWN','PWN']:
                Tcomb[ii][jj][wm]=T[wm][jj]
                Scomb[ii][jj][wm]=S[wm][jj]
            else:
                Tcomb[ii][jj][wm]=T[wm][0]
                Scomb[ii][jj][wm]=S[wm][0]
        if ii==jj:
            Ucomb[ii][jj]=get_U_from_x(NaN*ones(7))
        else:
            AM=array([[1,1,1,1,1,1,0],\
            [S['PWS'][ii],S['AWS'][ii],S['DWS'][ii],S['PWN'][jj],S['AWN'][jj],0,0],\
            [T['PWS'][ii],T['AWS'][ii],T['DWS'][ii],T['PWN'][jj],T['AWN'][jj],0,1]])

            x0=[U['PWS'][ii],U['AWS'][ii],U['DWS'][ii],U['PWN'][jj],U['AWN'][jj],mean(U['FW']),mean(U['Q'])]

            dv=-AM.dot(x0)

            Evec=[abs(U_std['PWS'][ii]/2),abs(U_std['AWS'][ii]/2),abs(U_std['DWS'][ii]/2),abs(U_std['PWN'][jj]/2),abs(U_std['AWN'][jj]/2),abs(mean(U_std['FW'])/2),abs(mean(U_std['Q'])/2)]
            E=diag(Evec)

            Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

            Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
            Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
            xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
            xsol_Ad=E.dot(xsol_prime)
            xbase=x0+xsol_Ad
            Ucomb[ii][jj]=get_U_from_x(xbase)


#######################################################################
############## Each SET OF 4 YEARS individually #######################
#######################################################################

num4=16
def get_U_S_T_from_WM_4yrcombos(Uin,Sin,Tin):
    U={}
    U_std={}
    S={}
    T={}
    for wm in Uin:
        U[wm]=zeros(num4)
        U_std[wm]=zeros(num4)
        S[wm]=zeros(num4)
        T[wm]=zeros(num4)
        for ii in range(num4):
            U[wm][ii]=mean(Uin[wm][ii:ii+4])
            U_std[wm][ii]=std(Uin[wm][ii:ii+4])
            S[wm][ii]=mean(Sin[wm][ii:ii+4])
            T[wm][ii]=mean(Tin[wm][ii:ii+4])

    return U,S,T,U_std

U_4yr,S_4yr,T_4yr,Ustd_4yr=get_U_S_T_from_WM_4yrcombos(U,S,T)

Usol_4yr={}
for ii in range(num4):
        AM=array([[1,1,1,1,1,1,0],\
        [S_4yr['PWS'][ii],S_4yr['AWS'][ii],S_4yr['DWS'][ii],S_4yr['PWN'][ii],S_4yr['AWN'][ii],S_4yr['FW'][ii],0],\
        [T_4yr['PWS'][ii],T_4yr['AWS'][ii],T_4yr['DWS'][ii],T_4yr['PWN'][ii],T_4yr['AWN'][ii],T_4yr['FW'][ii],1]])

        x0=[U_4yr['PWS'][ii],U_4yr['AWS'][ii],U_4yr['DWS'][ii],U_4yr['PWN'][ii],U_4yr['AWN'][ii],U_4yr['FW'][ii],U_4yr['Q'][ii]]

        dv=-AM.dot(x0)
        Evec=[abs(Ustd_4yr[wm][ii]/sqrt(4*4)) for wm in wmvec]
        E=diag(Evec)

        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        xbase=x0+xsol_Ad
        Usol_4yr[ii]=get_U_from_x(xbase)



#######################################################################
########### Now try combining 4 YEAR CHUNKS ###########################
#######################################################################

Ucomb_4yr={}
Tcomb_4yr={}
Scomb_4yr={}
for ii in range(num4):
    Ucomb_4yr[ii]={}
    Tcomb_4yr[ii]={}
    Scomb_4yr[ii]={}
    for jj in range(num4):
        Tcomb_4yr[ii][jj]={}
        Scomb_4yr[ii][jj]={}
        for wm in S:
            if wm in ['AWS','DWS','PWS']:
                Tcomb_4yr[ii][jj][wm]=T[wm][ii]
                Scomb_4yr[ii][jj][wm]=S[wm][ii]
            elif wm in ['AWN','PWN']:
                Tcomb_4yr[ii][jj][wm]=T[wm][jj]
                Scomb_4yr[ii][jj][wm]=S[wm][jj]
            else:
                Tcomb_4yr[ii][jj][wm]=T[wm][0]
                Scomb_4yr[ii][jj][wm]=S[wm][0]
        if ii==jj:
            Ucomb_4yr[ii][jj]=get_U_from_x(NaN*ones(7))
        else:
            AM=array([[1,1,1,1,1,1,0],\
            [S_4yr['PWS'][ii],S_4yr['AWS'][ii],S_4yr['DWS'][ii],S_4yr['PWN'][jj],S_4yr['AWN'][jj],0,0],\
            [T_4yr['PWS'][ii],T_4yr['AWS'][ii],T_4yr['DWS'][ii],T_4yr['PWN'][jj],T_4yr['AWN'][jj],0,1]])

            x0=[U_4yr['PWS'][ii],U_4yr['AWS'][ii],U_4yr['DWS'][ii],U_4yr['PWN'][jj],U_4yr['AWN'][jj],mean(U['FW']),mean(U['Q'])] #As before, full time series mean for FW and Q

            dv=-AM.dot(x0)

            Evec=[abs(Ustd_4yr['PWS'][ii]/4),abs(Ustd_4yr['AWS'][ii]/4),abs(Ustd_4yr['DWS'][ii]/4),abs(Ustd_4yr['PWN'][jj]/4),abs(Ustd_4yr['AWN'][jj]/4),abs(mean(Ustd_4yr['FW'])/4),abs(mean(Ustd_4yr['Q'])/4)]
            E=diag(Evec)

            Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

            Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
            Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
            xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
            xsol_Ad=E.dot(xsol_prime)
            xbase=x0+xsol_Ad
            Ucomb_4yr[ii][jj]=get_U_from_x(xbase)

#######################################################################
###### Define an error metric and plot it for different cases #########
#######################################################################
def get_err_from_19yrmean(U_choose,T_choose,S_choose):
    err={}
    err_vec_S=[abs(abs(Ubase_mean[wm]*Smean[wm]/35-U_choose[wm]*S_choose[wm]/35)/Ubase_mean[wm]*Smean[wm]/35) for wm in Ubase_mean]
    err['S']=sum(err_vec_S)/7*100
    err_vec_T=[abs(abs(Ubase_mean[wm]*Tmean[wm]/10-U_choose[wm]*T_choose[wm]/10)/Ubase_mean[wm]*Tmean[wm]/10) for wm in Ubase_mean]
    err['T']=sum(err_vec_T)/7*100
    err_vec_U=[abs(abs(Ubase_mean[wm]-U_choose[wm])/Ubase_mean[wm]) for wm in Ubase_mean]
    err['U']=sum(err_vec_U)/7*100
    return err


def vec2dic(test):
    dic={}
    for ii in range(len(test['AWS'])):
        dic[ii]={}
        for dd in test:
            dic[ii][dd]=test[dd][ii]
    return dic


err_1yr_overlap=[get_err_from_19yrmean(Usol[ii],vec2dic(T)[ii],vec2dic(S)[ii]) for ii in range(19)]
err_4yr_overlap=[get_err_from_19yrmean(Usol_4yr[ii],vec2dic(T_4yr)[ii],vec2dic(S_4yr)[ii]) for ii in range(num4)]

err_1yr_diff={}
for ii in range(19):
    for jj in range(19):
        if jj!=ii:
            err_1yr_diff[ii,jj]=get_err_from_19yrmean(Ucomb[ii][jj],Tcomb[ii][jj],Scomb[ii][jj])

err_4yr_diff={}
for ii in range(num4):
    for jj in range(num4):
        if abs(ii-jj)>3:
            err_4yr_diff[ii,jj]=get_err_from_19yrmean(Ucomb_4yr[ii][jj],Tcomb_4yr[ii][jj],Scomb_4yr[ii][jj])



def noresm_err_boxplot():
    f,ax=subplots(3,1,figsize=(5,6))
    diffcol='limegreen'
    overcol='teal'
    medcol='k'
    plvec=['U','S','T']
    for kk,axx in enumerate(ax):
        diffs=axx.boxplot([[err_1yr_diff[ii][plvec[kk]] for ii in err_1yr_diff],[err_4yr_diff[ii][plvec[kk]] for ii in err_4yr_diff]],positions=[0,2],widths=[0.5,0.5],patch_artist=True);
        overs=axx.boxplot([[ii[plvec[kk]] for ii in err_1yr_overlap],[ii[plvec[kk]] for ii in err_4yr_overlap]],positions=[1,3],widths=[0.5,0.5],patch_artist=True);
        axx.set_xticklabels('')
        axx.set_xticks([])
        axx.axvline(1.5,color='k')
        plt.setp(diffs['boxes'], color=diffcol)
        plt.setp(overs['boxes'], color=overcol)
        for kk in [overs,diffs]:
            plt.setp(kk['medians'], color=medcol,linewidth=2)
        axx.set_ylim(0,35)
    ax[2].plot([0,1],2*[1e3],color=diffcol,linewidth=5,label='Offset')
    ax[2].plot([0,1],2*[1e3],color=overcol,linewidth=5,label='Contemporaneous')
    ax[2].legend(loc=(0.13,-0.42),title='North and South initial conditions',ncol=2)
    ax[0].text(-0.2,40,'One-year means',fontsize=12)
    ax[0].text(1.75,40,'Four-year means',fontsize=12)
    ax[0].set_ylabel('Inverse model\nsolutions',fontsize=12)
    ax[1].set_ylabel('Salt\ntransport',fontsize=12)
    # ax[1].text(-1.75,-60,'Mean absolute percentage error from NorESM 19 year mean solution',fontsize=12,rotation='vertical')
    ax[2].set_ylabel('Temperature\ntransport',fontsize=12)
    savefig(figdir_paper+'/NorESM_subsamp.png',bbox_inches='tight')
    savefig(figdir_paper+'/NorESM_subsamp.pdf',bbox_inches='tight')

noresm_err_boxplot()

help(boxplot)

def plot_model_timeseries():
    f,ax=subplots(3,1,figsize=(20,10))
    [ax[0].plot(WM.resample(TIME='1Y').mean(dim='TIME').TIME.values,U[wm],'o-') for wm in wmvec[:5]]
    [ax[0].plot(WM.resample(TIME='6M').mean(dim='TIME').TIME.values[1::2][2:-1],U_4yr[wm],'*-') for wm in wmvec[:5]]
    [ax[1].plot(WM.resample(TIME='1Y').mean(dim='TIME').TIME.values,U[wm]*S[wm],'o-') for wm in wmvec[:5]]
    [ax[1].plot(WM.resample(TIME='6M').mean(dim='TIME').TIME.values[1::2][2:-1],U_4yr[wm]*S_4yr[wm],'*-') for wm in wmvec[:5]]
    [ax[2].plot(WM.resample(TIME='1Y').mean(dim='TIME').TIME.values,U[wm]*T[wm],'o-') for wm in wmvec[:]]
    [ax[2].plot(WM.resample(TIME='6M').mean(dim='TIME').TIME.values[1::2][2:-1],U_4yr[wm]*T_4yr[wm],'*-') for wm in wmvec[:]]


plot_model_timeseries()

#######################################################################
########## Plot the distributions #####################
#######################################################################

def get_Uc_wm(wm):
    Uc_wm=array([[Ucomb[ii][jj][wm] for ii in range(yrs)] for jj in range(yrs)]).flatten()
    Uc_wm=Uc_wm[~isnan(Uc_wm)]
    return Uc_wm


def plot_base_case_eachyear():
    f,axx=subplots(1,4,figsize=(9,2.5),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,1,1]))
    alf=0.85
    capi=7
    xshi1=0.68
    xshi2=1
    xshi3=1.32
    xwi=0.25
    whisnum=2
    boxo1=axx[0].boxplot([get_Uc_wm(wm) for wm in ['AWS','DWS']],patch_artist=True,notch=False,positions=(arange(2)+xshi3), showfliers=True,widths=xwi,whis=whisnum)
    boxi1=axx[0].boxplot([[Usol[kk][ll] for kk in Usol] for ll in ['AWS','DWS']],patch_artist=True,notch=False,positions=(arange(2)+xshi2), showfliers=True,widths=xwi,whis=whisnum)
    box1=axx[0].boxplot([U[kk] for kk in ['AWS','DWS']],patch_artist=True,notch=False,positions=(arange(2)+xshi1), showfliers=True,widths=xwi,whis=whisnum)
    axx[0].plot(arange(2)+xshi1,[Uinit_mean[ll] for ll in ['AWS','DWS']],'wo',zorder=100,mec='k')
    axx[0].plot(arange(2)+xshi2,[Ubase_mean[ll] for ll in ['AWS','DWS']],'w^',zorder=100,mec='k')
    boxo2=axx[1].boxplot([get_Uc_wm(ll) for ll in ['PWS','PWN','AWN']],patch_artist=True,notch=False,positions=(arange(3)+xshi3), showfliers=True,widths=xwi,whis=whisnum)
    boxi2=axx[1].boxplot([[Usol[kk][ll] for kk in Usol] for ll in ['PWS','PWN','AWN']],patch_artist=True,notch=False,positions=(arange(3)+xshi2), showfliers=True,widths=xwi,whis=whisnum)
    box2=axx[1].boxplot([U[kk] for kk in ['PWS','PWN','AWN']],patch_artist=True,notch=False,positions=(arange(3)+xshi1), showfliers=True,widths=xwi,whis=whisnum)
    axx[1].plot(arange(3)+xshi1,[Uinit_mean[ll] for ll in ['PWS','PWN','AWN']],'wo',zorder=100,mec='k')
    axx[1].plot(arange(3)+xshi2,[Ubase_mean[ll] for ll in ['PWS','PWN','AWN']],'w^',zorder=100,mec='k')
    boxo3=axx[2].boxplot(get_Uc_wm('FW'),patch_artist=True,notch=False,positions=(arange(1)+xshi3), showfliers=True,widths=xwi,whis=whisnum)
    boxi3=axx[2].boxplot([Usol[kk]['FW'] for kk in Usol],patch_artist=True,notch=False,positions=(arange(1)+xshi2), showfliers=True,widths=xwi,whis=whisnum)
    box3=axx[2].boxplot(U['FW'].values,patch_artist=True,notch=False,positions=(arange(1)+xshi1), showfliers=True,widths=xwi,whis=whisnum)
    axx[2].plot(arange(1)+xshi1,Uinit_mean['FW'],'wo',zorder=100,mec='k')
    axx[2].plot(arange(1)+xshi2,Ubase_mean['FW'],'w^',zorder=100,mec='k')
    boxo4=axx[3].boxplot(cp*rhow*get_Uc_wm('Q')/1e6,patch_artist=True,notch=False,positions=(arange(1)+xshi3), showfliers=True,widths=xwi,whis=whisnum)
    boxi4=axx[3].boxplot([cp*rhow*Usol[kk]['Q']/1e6 for kk in Usol],patch_artist=True,notch=False,positions=(arange(1)+xshi2), showfliers=True,widths=xwi,whis=whisnum)
    box4=axx[3].boxplot([cp*rhow*(U['Q'])/1e6],patch_artist=True,notch=False,positions=(arange(1)+xshi1), showfliers=True,widths=xwi,whis=whisnum)
    axx[3].plot(arange(1)+xshi1,cp*rhow*Uinit_mean['Q']/1e6,'wo',zorder=100,mec='k')
    axx[3].plot(arange(1)+xshi2,cp*rhow*Ubase_mean['Q']/1e6,'w^',zorder=100,mec='k')
    ylimi=20
    axx[0].set_ylim(-ylimi,ylimi)
    ylimi=4
    axx[1].set_ylim(-ylimi,ylimi)
    fwlim=0.2
    axx[2].set_ylim(-fwlim,fwlim)
    axx[3].set_ylim(-300,0)
    fsz=14
    axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    axx[3].set_ylabel('Heat flux [TW]',fontsize=fsz)
    for boxi in [boxo1,boxo2,boxo3,boxo4]:
        for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
                plt.setp(boxi[item], color='#02818a')
        plt.setp(boxi['fliers'],markeredgecolor='#02818a')
    for boxi in [boxi1,boxi2,boxi3,boxi4]:
        for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
                plt.setp(boxi[item], color='magenta')
        plt.setp(boxi['fliers'],markeredgecolor='m')
    for boxi in [box1,box2,box3,box4]:
        for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
                plt.setp(boxi[item], color='k')
        plt.setp(boxi['fliers'],markeredgecolor='k')
    for ii in range(3):
        axx[ii].axhline(0,color='k')
    axx[0].set_xticks(range(1,3))
    axx[0].set_xticklabels(['AWS','DWS'])
    axx[1].set_xticks(range(1,4))
    axx[1].set_xticklabels(['PWS','PWN','AWN'])
    axx[2].set_xticks(range(1,2))
    axx[2].set_xticklabels(['FW'])
    axx[3].set_xticks([1])
    axx[3].set_xticklabels('Q')

    savefig(figdir+'InvBudSol_NorESM_eachyear.png',bbox_inches='tight')
    savefig(figdir+'InvBudSol_NorESM_eachyear.pdf',bbox_inches='tight')

plot_base_case_eachyear()


#######################################################################
########## Get the impact on fresh water partitioning #####################
#######################################################################

epsilon=arange(0,1.1,0.1)
a_pwmat=zeros((len(epsilon),yrs,yrs))
b_pwmat=zeros((len(epsilon),yrs,yrs))
a_mean=zeros(len(epsilon))
b_mean=zeros(len(epsilon))
for ii,ee in enumerate(1-epsilon):
    a_mean[ii]=(ee*Ubase_mean['PWN']*(Smean['PWN']/Smean['AWS']-1)+Ubase_mean['PWS']*(Smean['PWS']/Smean['AWS']-1))/Ubase_mean['FW']
    b_mean[ii]=((1-ee)*Ubase_mean['PWN']*(Smean['PWN']/Smean['AWS']-1)+Ubase_mean['DWS']*(Smean['DWS']/Smean['AWS']-1))/Ubase_mean['FW']
    for jj in range(yrs):
        for kk in range(yrs):
                a_pwmat[ii,jj,kk]=(ee*Ucomb[jj][kk]['PWN']*(S['PWN'][kk]/S['AWS'][jj]-1)+Ucomb[jj][kk]['PWS']*(S['PWS'][jj]/S['AWS'][jj]-1))/Ucomb[jj][kk]['FW']
                b_pwmat[ii,jj,kk]=((1-ee)*Ucomb[jj][kk]['PWN']*(S['PWN'][kk]/S['AWS'][jj]-1)+Ucomb[jj][kk]['DWS']*(S['DWS'][jj]/S['AWS'][jj]-1))/Ucomb[jj][kk]['FW']

plot((1-a_pwmat-b_pwmat).flatten())


def plot_epsilon_dep():
    f,axx=subplots(1,2,figsize=(13,4.5),sharex=True)
    f.subplots_adjust(wspace=0.3)
    axx[0].set_xlim(0,1)
    for ii in range(19):
        for jj in range(19):
            if ii==jj:
                axx[0].plot(1-epsilon,a_pwmat[:,ii,ii],'red',zorder=3)
            else:
                axx[0].plot(1-epsilon,a_pwmat[:,ii,jj],'grey',linewidth=3,alpha=0.3,zorder=2)
    axx[0].plot(1-epsilon,a_mean,'k',linewidth=3,zorder=4)
    for ii in range(19):
        for jj in range(19):
            if ii==jj:
                axx[1].plot(1-epsilon,b_pwmat[:,ii,ii],'red',zorder=3)
            else:
                axx[1].plot(1-epsilon,b_pwmat[:,ii,jj],'grey',linewidth=3,alpha=0.3,zorder=2)
    axx[1].plot(1-epsilon,b_mean,'k',linewidth=3,zorder=4)
    axx[0].set_ylabel('$\mathbf{a}$, Fraction of FW in Polar Water South')
    axx[1].set_ylabel('$\mathbf{b}$, Fraction of FW in Deep Water South')
    for axi in axx:
        axi.set_yticks(arange(-0.5,1.5,0.25))
        axi.axhline(0,color='k',zorder=5)
        axi.set_xlabel('$\mathbf{1- \epsilon}$\n Fraction of PWN in PWS')
        axi.set_ylim(-0.5,1.4)
    savefig(figdir+'FWfrac_NorESM_yrdep.pdf',bbox_inches='tight')
    savefig(figdir+'FWfrac_NorESM_yrdep.png',bbox_inches='tight')
#

plot_epsilon_dep()
