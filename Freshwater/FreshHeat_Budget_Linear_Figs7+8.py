from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Linear/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs/'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################

WM=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')
Q=-65

def get_U_S_T_from_WM():
    U={}
    S={}
    T={}
    for wm in WM.WM:
        U[str(wm.values)]=float(WM['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        S[str(wm.values)]=float(WM['PSAL'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        T[str(wm.values)]=float(WM['PTMP'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)

    U['SI']=0.05
    U['FW']=0.05
    U['Q']=Q
    S['SI']=6
    S['FW']=0
    T['SI']=-90
    T['FW']=0
    T['Q']=1

    return U,S,T

U,S,T=get_U_S_T_from_WM()

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


#############################################################
# Testing sensitivities in Us, Ss, i.e. A and d
#############################################################
AM={}
# set up a few different cases,
AM['base']=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

AM['no_FW']=AM['base'].copy()
AM['no_SI']=AM['base'].copy()
AM['no_FW_or_SI']=AM['base'].copy()

AM['2xFW']=AM['base'].copy()
AM['2xSI']=AM['base'].copy()
# AM['2xPWN']=AM['base'].copy()

AM['NS_massbal']=array([[1,1,1,0,0,0.5,0.5,0],\
[0,0,0,1,1,0.5,0.5,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

AM['NS_massbal_adj']=AM['NS_massbal'].copy()
AM['NS_thruflo']=AM['NS_massbal'].copy()

AM['no_FW_thruflo']=AM['NS_massbal'].copy()
AM['no_SI_thruflo']=AM['NS_massbal'].copy()
AM['no_FW_or_SI_thruflo']=AM['NS_massbal'].copy()

AM['2xFW_thruflo']=AM['NS_massbal'].copy()
AM['2xSI_thruflo']=AM['NS_massbal'].copy()


AM['fresher_PWS_2']=array([[1,1,1,1,1,1,1,0],\
[S['PWS']-0.2,S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
AM['fresher_PW_2']=array([[1,1,1,1,1,1,1,0],\
[S['PWS']-0.2,S['AWS'],S['DWS'],S['PWN']-0.2,S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
AM['fresher_PWN_2']=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN']-0.2,S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

AM['fresher_PWS_2_thruflo']=array([[1,1,1,0,0,0.5,0.5,0],\
[0,0,0,1,1,0.5,0.5,0],\
[S['PWS']-0.2,S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
AM['fresher_PW_2_thruflo']=array([[1,1,1,0,0,0.5,0.5,0],\
[0,0,0,1,1,0.5,0.5,0],\
[S['PWS']-0.2,S['AWS'],S['DWS'],S['PWN']-0.2,S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
AM['fresher_PWN_2_thruflo']=array([[1,1,1,0,0,0.5,0.5,0],\
[0,0,0,1,1,0.5,0.5,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN']-0.2,S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

AM['fresher_PWS_5']=array([[1,1,1,1,1,1,1,0],\
[S['PWS']-0.5,S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
AM['fresher_PW_5']=array([[1,1,1,1,1,1,1,0],\
[S['PWS']-0.5,S['AWS'],S['DWS'],S['PWN']-0.5,S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
AM['fresher_PWN_5']=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN']-0.5,S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])


x0={}
x0['base']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]
x0['2xFW']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],2*U['FW'],U['SI'],Q]
x0['2xFW_thruflo']=x0['2xFW']
x0['2xSI']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],2*U['SI'],Q]
x0['2xSI_thruflo']=x0['2xSI']
x0['2xPWN']=[U['PWS'],U['AWS'],U['DWS'],2*U['PWN'],U['AWN'],U['FW'],U['SI'],Q]
x0['no_FW']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],0,U['SI'],Q]
x0['no_FW_thruflo']=x0['no_FW'].copy()
x0['no_SI']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],0,Q]
x0['no_SI_thruflo']=x0['no_SI'].copy()
x0['no_FW_or_SI']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],0,0,Q]
x0['no_FW_or_SI_thruflo']=x0['no_FW_or_SI'].copy()
x0['NS_massbal']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]
thruflo=1.6
x0['NS_massbal_adj']=[U['PWS'],U['AWS']-thruflo/2,U['DWS']-thruflo/2,U['PWN']+thruflo/2,U['AWN']+thruflo/2,U['FW'],U['SI'],Q]
x0['NS_thruflo']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]


dv={}
for kk in AM:
    if kk in x0:
        dv[kk]=-AM[kk].dot(x0[kk])
    else:
        dv[kk]=-AM[kk].dot(x0['base'])
    if 'thruflo' in kk:
        dv[kk][0]=dv[kk][0]+thruflo
        dv[kk][1]=dv[kk][1]-thruflo


Winv={}
Snorm=35
Tnorm=5
Winv['base']=diag([1,1/Snorm,1/Tnorm])
for kk in AM:
    if ('NS' in kk) | ('thruflo' in kk):
        Winv[kk]=diag([1,1,1/Snorm,1/Tnorm])

E={}
Qvar=5
Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,0.025,0.025,Qvar))
E['base']=diag(Evec)
E['no_FW']=diag(hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,1e-10,0.025,Qvar)))
E['no_FW_thruflo']=E['no_FW'].copy()
E['no_SI']=diag(hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,0.025,1e-10,Qvar)))
E['no_SI_thruflo']=E['no_SI'].copy()
E['no_FW_or_SI']=diag(hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,1e-10,1e-10,Qvar)))
E['no_FW_or_SI_thruflo']=E['no_FW_or_SI'].copy()



### this is set up so you only have to define an E or Winv if it is different than the base case:
xsol_Ad={}
P={}
for zz in AM:
    if zz in Winv:
        if zz in E:
            Umat,D,VmatT=linalg.svd(Winv[zz].dot(AM[zz].dot(E[zz])))
            P[zz]=diag(E[zz]-E[zz].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E[zz].dot(AM[zz].T))+linalg.inv(Winv[zz])).dot(AM[zz].dot(E[zz])))))
        else:
            Umat,D,VmatT=linalg.svd(Winv[zz].dot(AM[zz].dot(E['base'])))
            P[zz]=diag(E['base']-E['base'].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E['base'].dot(AM[zz].T))+linalg.inv(Winv[zz])).dot(AM[zz].dot(E['base'])))))
    else:
        if zz in E:
            Umat,D,VmatT=linalg.svd(Winv['base'].dot(AM[zz].dot(E[zz])))
            P[zz]=diag(E[zz]-E[zz].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E[zz].dot(AM[zz].T))+linalg.inv(Winv['base'])).dot(AM[zz].dot(E[zz])))))
        else:
            Umat,D,VmatT=linalg.svd(Winv['base'].dot(AM[zz].dot(E['base'])))
            P[zz]=diag(E['base']-E['base'].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E['base'].dot(AM[zz].T))+linalg.inv(Winv['base'])).dot(AM[zz].dot(E['base'])))))
    Lambda_inv = zeros((AM[zz].shape[0], AM[zz].shape[1])).T
    Lambda_inv[:AM[zz].shape[0], :AM[zz].shape[0]] = diag(1/D)
    if zz in Winv:
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv[zz].dot(dv[zz]))))
    else:
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv['base'].dot(dv[zz]))))
    if zz in E:
        xsol_Ad[zz]=E[zz].dot(xsol_prime)
    else:
        xsol_Ad[zz]=E['base'].dot(xsol_prime)


#ditto for x0
res={}
for xx in AM:
    if xx in x0:
        res[xx]=x0[xx]+xsol_Ad[xx]
    else:
        res[xx]=x0['base']+xsol_Ad[xx]


Ug=U.copy()
Sg=S.copy()
Tg=T.copy()
Ue=get_U_from_x(Evec)

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'k','Q':'limegreen'}

Ug['AWS']+Ug['PWS']+Ug['DWS']

def base_plot():
    hr=2
    f,axx=subplots(6,1,figsize=(9,10),constrained_layout=True,gridspec_kw=dict(height_ratios=[hr,1,hr,1,hr,1]))
    varlen=len(S)
    alf=0.75
    #U
    axx[0].bar(range(varlen)[:-2],[Ug[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2],yerr=[Ue[kk] for kk in Sg][:-2],capsize=10,alpha=alf)
    axx[0].set_ylim(-20,20)
    axx[1].set_ylim(-2.2,2.2)
    ax1=axx[0].twinx()
    ax1.bar([5,6],[Ug[kk] for kk in ['SI','FW',]],color=[coldic[kk] for kk in ['SI','FW',]],yerr=[Ue[kk] for kk in ['SI','FW']],capsize=10,alpha=alf)
    ax1.set_ylim(-0.15,0.15)
    ax1_1=axx[1].twinx()
    ax1_1.set_ylim(-0.01,0.01)
    # #UxS
    axx[2].bar(range(varlen),[Ug[kk]*Sg[kk] for kk in Sg],color=[coldic[kk] for kk in Sg],yerr=[Ue[kk]*Sg[kk] for kk in Sg],capsize=10,alpha=alf)
    axx[2].set_ylim(-750,750)
    axx[3].set_ylim(-75,75)
    # #UxT
    axx[4].bar(range(varlen+1),[Ug[kk]*Tg[kk] for kk in Ug],color=[coldic[kk] for kk in Ug],yerr=[Ue[kk]*Tg[kk] for kk in Ug],capsize=10,alpha=alf)
    axx[5].set_ylim(-7.5,7.5)
    colvec=['k','k','k']
    labvec=['No volume divergence','Enforce volume throughflow','Enforce volume balance at the boundaries']
    lsty=['-','dotted','--']
    for ii,case in enumerate(['base','NS_thruflo','NS_massbal']):
        Ux=get_U_from_x(res[case])
        axx[0].plot(range(varlen)[:5],[Ux[kk] for kk in S][:5],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        ax1.plot(range(varlen)[5:],[Ux[kk] for kk in S][5:],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        axx[2].plot(range(varlen),[(Ux[kk]*S[kk]) for kk in S],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        axx[4].plot(range(varlen+1),[(Ux[kk]*T[kk]) for kk in T],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        axx[1].plot(range(varlen)[:5],[Ux[kk]-Ug[kk] for kk in S][:5],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        ax1_1.plot(range(varlen)[5:],[Ux[kk]-Ug[kk] for kk in S][5:],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        axx[3].plot(range(varlen),[(Ux[kk]*S[kk])-(Ug[kk]*Sg[kk]) for kk in S],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])
        axx[5].plot(range(varlen+1),[(Ux[kk]*T[kk])-(Ug[kk]*Tg[kk]) for kk in T],'o-',label=labvec[ii],color=colvec[ii],linestyle=lsty[ii])


    for axi in [axx[0],axx[1]]:
        axi.axvline(4.5,color='k')

    for axxi in axx:
            axxi.set_xticks(range(varlen+1))
            axxi.set_xlim(-0.5,7.5)
            axxi.axhline(0,color='k')

    for axi in [axx[0],axx[2]]:
        axi.set_xticklabels(Sg.keys())
    axx[4].set_xticklabels(Tg.keys())
    for axi in [axx[1],axx[3],axx[5]]:
        axi.set_xticklabels('')

    fsz=14
    xp=-0.02
    f.text(xp,0.7,'Volume transport [Sv]',rotation='vertical',fontsize=fsz)
    f.text(xp,0.35,'Salinity transport, S x [Sv]',rotation='vertical',fontsize=fsz)
    f.text(xp,0,'Temperature transport [$^\circ$C x Sv]',rotation='vertical',fontsize=fsz)

    axx[0].legend(loc=(-0.04,1.075),ncol=3)
    savefig(figdir_paper+'InvBud_base.png',bbox_inches='tight')
    savefig(figdir_paper+'InvBud_base.pdf',bbox_inches='tight')

base_plot()

cp=3850
rhow=1000

def plot_a_case(casey,titi,ylimi):
    f,axx=subplots(1,2,figsize=(10,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[8,1]))
    varlen=len(S)
    Ux=get_U_from_x(res[casey])
    Ue=get_U_from_x(P[casey])
    Ubase=get_U_from_x(res['base'])
    alf=0.75
    #U
    axx[0].bar(range(varlen)[:-2],[Ux[kk]-Ubase[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2]),#yerr=[Ue[kk] for kk in Sg][:-2],capsize=10,alpha=alf)
    axx[0].set_ylim(-ylimi,ylimi)
    ax1=axx[0].twinx()
    ax1.bar([5,6],[Ux[kk]-Ubase[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']])#,yerr=[Ue[kk] for kk in ['SI','FW']],capsize=10,alpha=alf)
    ax1.set_ylim(-0.1,0.1)
    axx[0].axvline(4.5,color='k')
    fsz=14
    axx[0].set_ylabel('Anomaly [Sv]',fontsize=fsz)
    axx[1].set_ylabel('Anomaly [TW]',fontsize=fsz)
    axx[1].bar(0,cp*rhow*(Ux['Q']-Ubase['Q'])/1e6,color=coldic['Q'])#,yerr=Ue['Q'],capsize=10,alpha=alf)
    # colvec=['k']
    # labvec=['']
    # ii=0
    # Ux0=get_U_from_x(x0[casey])
    # Ubase0=get_U_from_x(x0['base'])
    # axx[0].plot(range(varlen)[:5],[Ux0[kk]-Ubase0[kk] for kk in S][:5],'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
    # ax1.plot([5,6],[Ux0[kk]-Ubase0[kk] for kk in ['SI','FW']],'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
    # axx[1].plot(0,cp*rhow*(Ux0['Q']-Ubase0['Q'])/1e6,'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')

    axx[0].set_xticks(range(varlen+1))
    axx[0].set_xticklabels(Sg.keys())
    axx[0].set_xlim(-0.5,6.5)
    axx[0].axhline(0,color='k')
    axx[1].set_xticks([0])
    axx[1].set_xticklabels('Q')
    axx[1].set_ylim(-15,15)
    axx[0].set_title(titi,fontsize=16)

    savefig(figdir_paper+'InvBud_case_'+casey+'_anom.png',bbox_inches='tight')
    savefig(figdir_paper+'InvBud_case_'+casey+'_anom.pdf',bbox_inches='tight')

plot_a_case('NS_massbal','Enforce mass balance at north and south',2)

plot_a_case('2xFW','Double fresh water input',2)
plot_a_case('2xSI','Double sea ice input',5)
plot_a_case('fresher_PWN_2','Freshen PWN input by 0.2',5)
plot_a_case('no_FW_or_SI','Remove all fresh water and sea ice',7.5)


def plot_base_case_simple(titi,ylimi):
    f,axx=subplots(1,2,figsize=(10,3),constrained_layout=True,gridspec_kw=dict(width_ratios=[8,1]))
    varlen=len(S)
    Ue=get_U_from_x(P['base'])
    Ubase=get_U_from_x(res['base'])
    alf=0.75
    #U
    axx[0].bar(range(varlen)[:-2],[Ubase[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2],yerr=[Ue[kk] for kk in Sg][:-2],capsize=10,alpha=alf)
    axx[0].set_ylim(-ylimi,ylimi)
    ax1=axx[0].twinx()
    ax1.bar([5,6],[Ubase[kk] for kk in ['SI','FW']],color=[coldic[kk] for kk in ['SI','FW']],yerr=[Ue[kk] for kk in ['SI','FW']],capsize=10,alpha=alf)
    ax1.set_ylim(-0.1,0.1)
    axx[0].axvline(4.5,color='k')
    fsz=14
    axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    axx[1].set_ylabel('Heat flux [TW]',fontsize=fsz)
    axx[1].bar(0,cp*rhow*(Ubase['Q'])/1e6,color=coldic['Q'],yerr=Ue['Q'],capsize=10,alpha=alf)
    # colvec=['k']
    # labvec=['']
    # ii=0
    # Ux0=get_U_from_x(x0[casey])
    # Ubase0=get_U_from_x(x0['base'])
    # axx[0].plot(range(varlen)[:5],[Ux0[kk]-Ubase0[kk] for kk in S][:5],'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
    # ax1.plot([5,6],[Ux0[kk]-Ubase0[kk] for kk in ['SI','FW']],'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
    # axx[1].plot(0,cp*rhow*(Ux0['Q']-Ubase0['Q'])/1e6,'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')

    axx[0].set_xticks(range(varlen+1))
    axx[0].set_xticklabels(Sg.keys())
    axx[0].set_xlim(-0.5,6.5)
    axx[0].axhline(0,color='k')
    axx[1].set_xticks([0])
    axx[1].set_xticklabels('Q')
    axx[1].set_ylim(-275,0)
    axx[0].set_title(titi,fontsize=16)

    savefig(figdir_paper+'InvBud_base_simple.png',bbox_inches='tight')
    savefig(figdir_paper+'InvBud_base_simple.pdf',bbox_inches='tight')


plot_base_case_simple('Base case solution',20)

# def plot_no_FWSI(casey):
#         hr=2
#         f,axx=subplots(1,2,figsize=(8,4),constrained_layout=True,gridspec_kw=dict(width_ratios=[7,1]))
#         varlen=len(S)
#         Ux=get_U_from_x(res[casey])
#         Ue=get_U_from_x(P[casey])
#         alf=0.75
#         #U
#         axx[0].bar(range(varlen)[:-2],[Ux[kk] for kk in Sg][:-2],color=[coldic[kk] for kk in Sg][:-2],yerr=[Ue[kk] for kk in Sg][:-2],capsize=10,alpha=alf)
#         axx[0].set_ylim(-20,20)
#         axx[0].set_ylabel('Volume transport [Sv]')
#         axx[1].set_ylabel('Heat flux [TW]')
#         axx[1].bar(0,cp*rhow*Ux['Q']/1e6,color=coldic['Q'],yerr=Ue['Q'],capsize=10,alpha=alf)
#         colvec=['k']
#         labvec=['']
#         for ii,case in enumerate(['base']):
#             Ux=get_U_from_x(res[case])
#             axx[0].plot(range(varlen)[:5],[Ux[kk] for kk in S][:5],'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
#             # ax1.plot([5,6],[Ux[kk] for kk in ['FW','SI']],'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
#             axx[1].plot(0,cp*rhow*Ux['Q']/1e6,'*',label=labvec[ii],color=colvec[ii],markersize=16,mfc='w')
#
#         axx[0].set_xticks(range(varlen+1))
#         axx[0].set_xticklabels(Sg.keys())
#         axx[0].set_xlim(-0.5,4.5)
#         axx[0].axhline(0,color='k')
#         axx[1].set_xticks([0])
#         axx[1].set_xticklabels('Q')
#         axx[1].set_ylim(-275,0)
#
#         savefig(figdir_paper+'InvBud_noFWSI.png',bbox_inches='tight')
#         savefig(figdir_paper+'InvBud_noFWSI.pdf',bbox_inches='tight')
#
# plot_no_FWSI('no_FW_or_SI')

# def complot_sens(casename):
#     # Ue=get_U_from_x(Evec[:-1])
#     Ug=U.copy()
#     Sg=S.copy()
#     Tg=T.copy()
#     f,axx=subplots(3,2,figsize=(16,10),sharex='col')
#     varlen=len(S)
#     # alf=0.4
#     # #U
#     # axx[0,0].bar(range(varlen)[:-2],[Ug[kk] for kk in Ug][:-2],color=colvec,yerr=[Ue[kk] for kk in Ug][:-2],capsize=10,alpha=alf)
#     axx[0,0].set_ylim(-25,25)
#     # axx[0,0].axvline(4.5,color='k')
#     ax1=axx[0,0].twinx()
#     # ax1.bar([5,6],[Ug[kk] for kk in ['FW','SI']],color=colvec,yerr=[Ue[kk] for kk in ['FW','SI']],capsize=10,alpha=alf)
#     ax1.set_ylim(-0.2,0.2)
#     axx[0,1].set_ylim(-5,5)
#     ax1_1=axx[0,1].twinx()
#     ax1_1.set_ylim(-0.1,0.1)
#     # #UxS
#     # axx[1,0].bar(range(varlen),[Ug[kk]*Sg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Sg[kk] for kk in Ug],capsize=10,alpha=alf)
#     axx[1,0].set_ylim(-750,750)
#     # #UxT
#     # axx[2,0].bar(range(varlen),[Ug[kk]*Tg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Tg[kk] for kk in Ug],capsize=10,alpha=alf)
#     # axx[2,0].bar(varlen,xbase[-1],color=colvec,yerr=Evec[-1],capsize=10,alpha=alf)
#     ii=0
#     for case in res:
#         if (casename in case) | ('base' in case):
#             Ux=get_U_from_x(res[case])
#             if case in x0:
#                 Ux0=get_U_from_x(x0[case])
#             else:
#                 Ux0=get_U_from_x(x0['base'])
#
#             axx[0,0].plot(range(varlen)[:5],[Ux[kk] for kk in S][:5],'o-',label=case,color='C'+str(ii))
#             axx[0,0].plot(range(varlen)[:5],[Ux0[kk] for kk in S][:5],'o--',color='C'+str(ii))
#
#             ax1.plot(range(varlen)[5:],[Ux[kk] for kk in S][5:],'o-',label=case,color='C'+str(ii))
#             ax1.plot(range(varlen)[5:],[Ux0[kk] for kk in S][5:],'o--',color='C'+str(ii))
#
#             axx[1,0].plot(range(varlen),[Ux[kk]*S[kk] for kk in S],'o-',label=case,color='C'+str(ii))
#             axx[2,0].plot(range(varlen+1),[Ux[kk]*T[kk] for kk in T],'o-',label=case,color='C'+str(ii))
#             axx[2,0].plot(varlen,res[case][-1],'o-',label=case,color='C'+str(ii))
#
#             axx[0,1].plot(range(varlen)[:5],[Ux[kk]-Ug[kk] for kk in S][:5],'o-',label=case,color='C'+str(ii))
#             axx[0,1].plot(range(varlen)[:5],[Ux[kk]-Ux0[kk] for kk in S][:5],'o--',color='C'+str(ii))
#             ax1_1.plot(range(varlen)[5:],[Ux[kk]-Ug[kk] for kk in S][5:],'o-',label=case,color='C'+str(ii))
#             ax1_1.plot(range(varlen)[5:],[Ux[kk]-Ux0[kk] for kk in S][5:],'o--',color='C'+str(ii))
#             axx[1,1].plot(range(varlen),[(Ux[kk]*S[kk])-(Ug[kk]*Sg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
#             axx[1,1].plot(range(varlen),[(Ux[kk]*S[kk])-(Ux0[kk]*Sg[kk]) for kk in S],'o--',color='C'+str(ii))
#
#             axx[2,1].plot(range(varlen+1),[(Ux[kk]*T[kk])-(Ug[kk]*Tg[kk]) for kk in T],'o-',label=case,color='C'+str(ii))
#             axx[2,1].plot(range(varlen+1),[(Ux[kk]*T[kk])-(Ux0[kk]*Tg[kk]) for kk in T],'o--',color='C'+str(ii))
#
#             ii+=1
#     for axi in [axx[0,0],axx[0,1]]:
#         axi.axvline(4.5,color='k')
#     axx[0,1].set_title('Anomalies from initial conditions')
#     # ax2.legend(loc='upper right',)
#     axx[1,1].legend(loc=(1.05,0))
#     for axesi in axx:
#         for axxi in axesi:
#             axxi.set_xticks(range(varlen+1))
#             axxi.set_xticklabels(Tg.keys())
#             axxi.axhline(0,color='k')
#     axx[0,0].set_ylabel('Transport, U [Sv]')
#     axx[1,0].set_ylabel('Salinity budget term\n UxS [Sv]')
#     axx[2,0].set_ylabel('Temperature budget term\n UxT [$^\circ$C Sv]')
#
#     savefig(figdir+'Budget_Linear_Ad_Var_'+casename+'.png',bbox_inches='tight')
#
# complot_sens('base')


##### Below was for testing sensitivity to weights, keeping it around, but don't want it getting in the way.

# #############################################################
# # Plotting infrastructure
# #############################################################
#
#
# def BudgetOut(x):
#     mall=x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]
#     sall=(x[0]*S['PWS']+x[1]*S['AWS']+x[2]*S['DWS']+x[3]*S['PWN']+x[4]*S['AWN']+x[6]*S['SI'])
#     tall=(x[0]*T['PWS']+x[1]*T['AWS']+x[2]*T['DWS']+x[3]*T['PWN']+x[4]*T['AWN']+x[5]*T['FW']+x[6]*T['SI']+x[7])
#     return mall,sall,tall
#
# mb,sb,tb=BudgetOut(xbase)
# mg,sg,tg=BudgetOut(xbase+xsol)
#
# mb,mg
# sb,sg
# tb,tg
#
# def plot_budget_out(casename):
#     f,axx=subplots(3,1,figsize=(8,9),sharex=True)
#     ii=0
#     for kk in res:
#         if casename in kk:
#             mall,sall,tall=BudgetOut(res[kk].x)
#             axx[0].plot(ii,mall,'o',label=kk)
#             axx[1].plot(ii,sall,'o',label=kk)
#             axx[2].plot(ii,tall,'o',label=kk)
#             ii+=1
#     for axi in axx:
#         axi.axhline(0,color='grey')
#     axx[1].legend(loc=(1.05,0))
#     axx[2].set_xticklabels('')
#     axx[0].set_ylabel('Volume balance [Sv]')
#     axx[1].set_ylabel('Salinity balance [S x Sv]')
#     axx[2].set_ylabel('Temperature balance [$^\circ$C x Sv]')
#     axx[0].set_title('Balance residuals for '+casename)
#     savefig(figdir+'Gamma_dependence_'+casename+'.png',bbox_inches='tight')
#
# colvec='grey'
#
# def barplot(casename):
#     Ue=get_U_from_x(Evec[:-1])
#     Ug=U.copy()
#     Sg=S.copy()
#     Tg=T.copy()
#     f,axx=subplots(3,2,figsize=(16,10),sharex='col')
#     varlen=len(Ug)
#     alf=0.4
#     #U
#     axx[0,0].bar(range(varlen)[:-2],[Ug[kk] for kk in Ug][:-2],color=colvec,yerr=[Ue[kk] for kk in Ug][:-2],capsize=10,alpha=alf)
#     axx[0,0].set_ylim(-25,25)
#     axx[0,0].axvline(4.5,color='k')
#     ax1=axx[0,0].twinx()
#     ax1.bar([5,6],[Ug[kk] for kk in ['FW','SI']],color=colvec,yerr=[Ue[kk] for kk in ['FW','SI']],capsize=10,alpha=alf)
#     ax1.set_ylim(-0.2,0.2)
#     #UxS
#     axx[1,0].bar(range(varlen),[Ug[kk]*Sg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Sg[kk] for kk in Ug],capsize=10,alpha=alf)
#     axx[1,0].set_ylim(-750,750)
#     #UxT
#     axx[2,0].bar(range(varlen),[Ug[kk]*Tg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Tg[kk] for kk in Ug],capsize=10,alpha=alf)
#     axx[2,0].bar(varlen,xbase[-1],color=colvec,yerr=Evec[-1],capsize=10,alpha=alf)
#     ii=0
#     ax1_1=axx[0,1].twinx()
#     for case in res:
#         if casename in case:
#             Ux=get_U_from_x(res[case])
#             axx[0,0].plot(range(varlen)[:5],[Ux[kk] for kk in U][:5],'o-',label=case,color='C'+str(ii))
#             ax1.plot(range(varlen)[5:],[Ux[kk] for kk in U][5:],'o-',label=case,color='C'+str(ii))
#             axx[1,0].plot(range(varlen),[Ux[kk]*S[kk] for kk in S],'o-',label=case,color='C'+str(ii))
#             axx[2,0].plot(range(varlen),[Ux[kk]*T[kk] for kk in S],'o-',label=case,color='C'+str(ii))
#             axx[2,0].plot(varlen,res[case][-1],'o-',label=case,color='C'+str(ii))
#             axx[0,1].plot(range(varlen)[:5],[Ux[kk]-Ug[kk] for kk in U][:5],'o-',label=case,color='C'+str(ii))
#
#             ax1_1.plot(range(varlen)[5:],[Ux[kk]-Ug[kk] for kk in U][5:],'o-',label=case,color='C'+str(ii))
#             axx[1,1].plot(range(varlen),[(Ux[kk]*S[kk])-(Ug[kk]*Sg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
#             axx[2,1].plot(range(varlen),[(Ux[kk]*T[kk])-(Ug[kk]*Tg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
#             axx[2,1].plot(varlen,res[case][-1]-xbase[-1],'o-',label=case,color='C'+str(ii))
#             ii+=1
#     for axi in [axx[0,0],axx[1,0],axx[0,1]]:
#         axi.axvline(4.5,color='k')
#     axx[0,1].set_title('Anomalies from initial conditions')
#     # ax2.legend(loc='upper right',)
#     axx[1,1].legend(loc=(1.05,0))
#     for axesi in axx:
#         for axxi in axesi:
#             axxi.set_xticks(range(varlen+1))
#             axxi.set_xticklabels(Tg.keys())
#             axxi.axhline(0,color='k')
#     axx[0,0].set_ylabel('Transport, U [Sv]')
#     axx[1,0].set_ylabel('Salinity budget term\n UxS [Sv]')
#     axx[2,0].set_ylabel('Temperature budget term\n UxT [$^\circ$C Sv]')
#
#     savefig(figdir+'FWBudget_Linear_'+casename+'.png',bbox_inches='tight')

#############################################################
# Set up A and d
#############################################################
#
# Amat=array([[1,1,1,1,1,1,1,0],\
# [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
# [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
#
# Q=-65
# xbase=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]
#
# Usum=0
# SUsum=0
# TUsum=Q
# for wm in U:
#     Usum=Usum+U[wm]
#     SUsum=SUsum+U[wm]*S[wm]
#     TUsum=TUsum+U[wm]*T[wm]
#
# dvec=[-Usum,-SUsum,-TUsum]

#############################################################
# Testing W, E
#############################################################
#
# # set up a few different cases,
#
# Winv={}
# Snorm=35
# Tnorm=5
# Winv['W=base']=diag([1,1/Snorm,1/Tnorm])
# Winv['W=1']=diag([1,1,1])
# Winv['W=1000']=diag([1000,1000,1000])
#
#
# E={}
# Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,0.025,0.025,10))
# E['E=base']=diag(Evec)
# Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values/100,0.025,0.025,10))
# E['E=U_WM/100']=diag(Evec)
# Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,0.025,0.025,10/100))
# E['E=q_T/100']=diag(Evec)
# Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,0.025,0.025,10))
# E['E=base/1e4']=diag(Evec/1e4)
#
# xsol={}
# for xx in E:
#     for yy in Winv:
#         Umat,D,VmatT=linalg.svd(Winv[yy].dot(Amat.dot(E[xx])))
#         Lambda = zeros((Amat.shape[0], Amat.shape[1]))
#         # populate Sigma with n x n diagonal matrix
#         Lambda[:Amat.shape[0], :Amat.shape[0]] = diag(D)
#         Lambda_inv=Lambda.T
#         Lambda_inv[:Amat.shape[0], :Amat.shape[0]] = diag(1/D)
#
#         xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv[yy].dot(dvec))))
#
#         xsol[xx+'; '+yy]=E[xx].dot(xsol_prime)
#
#
# res={}
# for xx in xsol:
#     res[xx]=xbase+xsol[xx]
#
# barplot('E=base; W=base')
#
# barplot('W=base')
