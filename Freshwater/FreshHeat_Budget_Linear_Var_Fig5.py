from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Linear/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################

WM=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')

cp=3850
rhow=1000
tera=10**12
#Noresm (taking sea ice into account)
Q=-251*tera/rhow/cp/1e6 #for the Sverdrups

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

dv=-AM['base'].dot(x0)

Snorm=35
Tnorm=5
Winv=diag([1,1/Snorm,1/Tnorm])

abs(array([U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN']])/5)

Qvar=10
Evec=hstack((abs(array([U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN']])/5),0.02,0.02,Qvar))
E=diag(Evec)

zz='base'
Umat,D,VmatT=linalg.svd(Winv.dot(AM[zz].dot(E)))

Lambda_inv = zeros((AM[zz].shape[0], AM[zz].shape[1])).T
Lambda_inv[:AM[zz].shape[0], :AM[zz].shape[0]] = diag(1/D)
xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
xsol_Ad=E.dot(xsol_prime)
res=x0+xsol_Ad
P=diag(E-E.dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E.dot(AM[zz].T))+linalg.inv(Winv)).dot(AM[zz].dot(E)))))

Ue=get_U_from_x(Evec)

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'k','Q':'limegreen'}


def plot_base_case_simple():
    f,axx=subplots(1,4,figsize=(13,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,2,1]))
    Ue=get_U_from_x(P)
    Ubase=get_U_from_x(res)
    alf=0.75
    capi=7
    #U
    axx[0].bar(range(2),[Ubase[kk] for kk in ['AWS','DWS']],color=[coldic[kk] for kk in ['AWS','DWS']],yerr=[Ue[kk] for kk in ['AWS','DWS']],capsize=capi,alpha=alf)
    axx[0].plot(range(2),[Ug[kk] for kk in ['AWS','DWS']],'o',color='k')

    ylimi=20
    axx[0].set_ylim(-ylimi,ylimi)
    ylimi=4.5
    axx[1].set_ylim(-ylimi,ylimi)
    axx[1].bar(range(3),[Ubase[kk] for kk in ['PWS','PWN','AWN']],color=[coldic[kk] for kk in ['PWS','PWN','AWN']],yerr=[Ue[kk] for kk in ['PWS','PWN','AWN']],capsize=capi,alpha=alf)
    axx[1].plot(range(3),[Ug[kk] for kk in ['PWS','PWN','AWN']],'o',color='k')

    axx[2].bar(range(2),[Ubase[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']],yerr=[Ue[kk] for kk in ['SI','FW']],capsize=capi,alpha=alf)
    axx[2].plot(range(2),[Ug[kk] for kk in ['SI','FW']],'o',color='k')
    fwlim=0.125
    axx[2].set_ylim(-fwlim,fwlim)

    fsz=14
    axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    axx[3].set_ylabel('Heat flux [TW]',fontsize=fsz)
    axx[3].bar(0,cp*rhow*(Ubase['Q'])/1e6,color=coldic['Q'],yerr=cp*rhow*Ue['Q']/1e6,capsize=capi,alpha=alf)
    axx[3].plot(0,cp*rhow*(Ug['Q'])/1e6,'o',color='k')

    # axx[0].set_xticks(range(varlen+1))
    # axx[0].set_xticklabels(Sg.keys())
    # axx[0].set_xlim(-0.5,6.5)
    for ii in range(3):
        axx[ii].axhline(0,color='k')
        axx[ii].set_xticks(range(2))
    axx[0].set_xticklabels(['AWS','DWS'])
    axx[1].set_xticks(range(3))
    axx[1].set_xticklabels(['PWS','PWN','AWN'])
    axx[2].set_xticklabels(['SI','FW'])
    axx[3].set_xticks([0])
    axx[3].set_xticklabels('Q')


    savefig(figdir_paper+'InvBudSol_SI0.png',bbox_inches='tight')
    savefig(figdir_paper+'InvBudSol_SI0.pdf',bbox_inches='tight')

plot_base_case_simple()
Ubase=get_U_from_x(res)


basediff=[(kk,Ubase[kk]-U[kk]) for kk in Ubase]
basediff


U['PWN']+U['AWN']

U['AWS']+U['DWS']+U['PWS']

0.095+0.12+0.065+0.007+0.006

S['PWS']
S['PWN']


##################################################################################
# Calculate fraction of fresh water vs. other water masses that goes into each limb
#################################################################################

#fraction of PWN in DWS limb
epsilon=arange(0,1.1,0.1)

#fraction of FW in PWS, as a function of epsilon
a=((1-epsilon)*Ubase['PWN']*(S['PWN']/S['AWS']-1)+Ubase['PWS']*(S['PWS']/S['AWS']-1))/(Ubase['FW']+Ubase['SI'])
#fraction of FW in DWS, as a function of epsilon
b=(epsilon*Ubase['PWN']*(S['PWN']/S['AWS']-1)+Ubase['DWS']*(S['DWS']/S['AWS']-1))/(Ubase['FW']+Ubase['SI'])

a

plot(epsilon,a)
plot(epsilon,b)

# def plot_base_case_simple():
#     f,axx=subplots(1,2,figsize=(10,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[8,1]))
#     varlen=len(S)
#     Ue=get_U_from_x(P)
#     Ubase=get_U_from_x(res)
#     alf=0.75
#     capi=7
#     #U
#     axx[0].bar(range(varlen)[:-2],[Ubase[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2],yerr=[Ue[kk] for kk in Sg][:-2],capsize=capi,alpha=alf)
#     axx[0].plot(range(varlen)[:-2],[Ug[kk] for kk in Sg][:-2],'o',color='k')
#     # for ii,kk in enumerate(['AWS','PWS','DWS','AWN','PWN']):
#     ylimi=20
#     axx[0].set_ylim(-ylimi,ylimi)
#     ax1=axx[0].twinx()
#     ax1.bar([5,6],[Ubase[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']],yerr=[Ue[kk] for kk in ['SI','FW']],capsize=capi,alpha=alf)
#     ax1.plot([5,6],[Ug[kk] for kk in ['SI','FW']],'o',color='k')
#     fwlim=0.125
#     ax1.set_ylim(-fwlim,fwlim)
#     axx[0].axvline(4.5,color='k')
#     fsz=14
#     axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
#     axx[1].set_ylabel('Heat flux [TW]',fontsize=fsz)
#     axx[1].bar(0,cp*rhow*(Ubase['Q'])/1e6,color=coldic['Q'],yerr=cp*rhow*Ue['Q']/1e6,capsize=capi,alpha=alf)
#     axx[1].plot(0,cp*rhow*(Ug['Q'])/1e6,'o',color='k')
#
#     axx[0].set_xticks(range(varlen+1))
#     axx[0].set_xticklabels(Sg.keys())
#     axx[0].set_xlim(-0.5,6.5)
#     axx[0].axhline(0,color='k')
#     ax1.axhline(0,color='k')
#     axx[1].set_xticks([0])
#     axx[1].set_xticklabels('Q')
#     axx[1].set_ylim(-400,0)
#
#     savefig(figdir_paper+'InvBudSol.png',bbox_inches='tight')
#     savefig(figdir_paper+'InvBudSol.pdf',bbox_inches='tight')



##### Look into how much Sea ice properties matter
sivar={}

for S_SI in range(0,10,2):
    sivar[S_SI]={}
    for T_SI in range(-90,0,10):

        AM=array([[1,1,1,1,1,1,1,0],\
        [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S_SI,0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T_SI,1]])

        Evec=hstack((5*[1],0.01,0.01,10))
        E=diag(Evec)

        x0=res
        dv=-AM.dot(x0)

        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        sivar[S_SI][T_SI]=x0+xsol_Ad


def plot_varTS():
    f,axx=subplots(3,2,figsize=(8,6),constrained_layout=True,gridspec_kw=dict(width_ratios=[8,1],height_ratios=[2,1,1]))
    varlen=len(S)
    Ue=get_U_from_x(P)
    Ubase=get_U_from_x(res)
    alf=0.75
    #U
    axx[0,0].bar(range(varlen)[:-2],[Ubase[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2])
    axx[0,0].plot(range(varlen)[:-2],[Ug[kk] for kk in Sg][:-2],'o',color='k')

    ylimi=25
    axx[0,0].set_ylim(-ylimi,ylimi)
    ax1=axx[0,0].twinx()
    ax1.bar([5,6],[Ubase[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']])
    ax1.plot([5,6],[Ug[kk] for kk in ['SI','FW']],'o',color='k')
    # for aa in varsol:
    #     for bb in varsol[aa]:
    #         Utmp=get_U_from_x(varsol[aa][bb])
    #         axx[0,0].plot(arange(varlen)[:-2],[Utmp[kk] for kk in Sg][:-2],'o',color='y',mec='k')
    #         ax1.plot(array([5,6]),[Utmp[kk] for kk in ['SI','FW']],'o',color='y',mec='k')
    #         axx[0,1].plot(0,cp*rhow*(Utmp['Q'])/1e6,'o',color='y',mec='k')
    #         if aa==S['PWS']:
    #             axx[0,0].plot(arange(varlen)[:-2]+0.4,[Utmp[kk] for kk in Sg][:-2],'o',color=coldic['PWN'],mec='k')
    #             ax1.plot(array([5,6])+0.4,[Utmp[kk] for kk in ['SI','FW']],'o',color=coldic['PWN'],mec='k')
    #             axx[0,1].plot(+0.4,cp*rhow*(Utmp['Q'])/1e6,'o',color=coldic['PWN'],mec='k')
    #         if bb==S['PWN']:
    #             axx[0,0].plot(arange(varlen)[:-2]-0.4,[Utmp[kk] for kk in Sg][:-2],'o',color=coldic['PWS'],mec='k')
    #             ax1.plot(array([5,6])-0.4,[Utmp[kk] for kk in ['SI','FW']],'o',color=coldic['PWS'],mec='k')
    #             axx[0,1].plot(-0.4,cp*rhow*(Utmp['Q'])/1e6,'o',color=coldic['PWS'],mec='k')
    for aa in sivar:
        for bb in sivar[aa]:
            Utmp=get_U_from_x(sivar[aa][bb])
            if aa==S['SI']:
                        axx[0,0].plot(arange(varlen)[:-2]+0.4,[Utmp[kk] for kk in Sg][:-2],'o',color='g',mec='k')
                        ax1.plot(array([5,6])+0.4,[Utmp[kk] for kk in ['SI','FW']],'o',color='g',mec='k')
                        axx[0,1].plot(+0.4,cp*rhow*(Utmp['Q'])/1e6,'o',color='g',mec='k')
            if bb==T['SI']:
                axx[0,0].plot(arange(varlen)[:-2]-0.4,[Utmp[kk] for kk in Sg][:-2],'o',color='pink',mec='k')
                ax1.plot(array([5,6])-0.4,[Utmp[kk] for kk in ['SI','FW']],'o',color='pink',mec='k')
                axx[0,1].plot(-0.4,cp*rhow*(Utmp['Q'])/1e6,'o',color='pink',mec='k')
    fwlim=0.2
    ax1.set_ylim(-fwlim,fwlim)
    axx[0,0].axvline(4.5,color='k')
    fsz=14
    axx[0,0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    #Q
    axx[0,1].set_ylabel('Heat flux [TW]',fontsize=fsz)
    axx[0,1].bar(0,cp*rhow*(Ubase['Q'])/1e6,color=coldic['Q'])
    axx[0,1].plot(0,cp*rhow*(Ug['Q'])/1e6,'o',color='k')

    axx[0,0].set_xticks(range(varlen+1))
    axx[0,0].set_xticklabels(Sg.keys())
    axx[0,0].set_xlim(-0.5,6.5)
    axx[0,0].axhline(0,color='k')

    axx[0,1].set_xticks([0])
    axx[0,1].set_xticklabels('Q')
    axx[0,1].set_ylim(-400,0)

    #S
    axx[1,0].bar(range(varlen)[:-2],[Sg[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2])
    axx[1,0].set_ylim(33,36)
    ax2=axx[1,0].twinx()
    ax2.bar([5,6],[Sg[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']],yerr=[0,0])
    axx[1,0].axvline(4.5,color='k')
    axx[1,0].set_xticks(range(varlen+1))
    axx[1,0].set_xticklabels(Sg.keys())
    axx[1,1].axis('off')
    axx[1,0].set_xlim(-0.5,6.5)

    #T
    axx[2,0].bar(range(varlen)[:-2],[Tg[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2])
    ax3=axx[2,0].twinx()
    ax3.bar([5,6],[Tg[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']])
    axx[2,0].axvline(4.5,color='k')
    axx[2,0].set_xticks(range(varlen+1))
    axx[2,0].set_xticklabels(Sg.keys())
    axx[2,0].set_xlim(-0.5,6.5)
    axx[2,1].axis('off')


plot_varTS()
