from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Linear/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################
#ocean water masses
WM=xr.open_dataset(datadir+'FW_WM/OSNAP2014-18_Tsub2020_WM_2008.nc')
WM_mb=xr.open_dataset(datadir+'FW_WM/OSNAP2014-18_WM_mb_2008.nc')

#surface fluxes
noresm=xr.open_dataset(datadir+'NorESM/NorESM_freshwater_fromspatial.nc')
hf1=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_200001-200912.nc')
hf2=xr.open_dataset(datadir+'NorESM/NorESM2-LM_omip2_NordicSeas_heatloss_201001-201812.nc')
noresm_hf=xr.concat([hf1,hf2],dim='time')
noresm['hf_int']=('time',-noresm_hf.NORDIC_hflx.values)


cp=3850
rhow=1025
tera=10**12

def get_U_S_T_from_WM(WM):
    U={}
    U_std={}
    S={}
    T={}
    for wm in WM.WM:
        U[str(wm.values)]=float(WM['TRANS'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        U_std[str(wm.values)]=float(WM['TRANS'].sel(WM=wm).std('TIME').values)
        S[str(wm.values)]=float(WM['SA'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)
        T[str(wm.values)]=float(WM['CT'].sel(WM=wm).groupby('TIME.month').mean('TIME').mean(dim='month').values)

    #note: minus sign is a carry-over from plotting atm fluxes
    U['FW']=(-noresm['runoff']-noresm['ep_int']-noresm['mltfrz']).mean(dim='time').values
    U['Q']=-noresm['hf_int'].mean(dim='time').values*tera/rhow/cp/1e6 #1e6 for the Sverdrups
    U_std['FW']=(-noresm['runoff']-noresm['ep_int']-noresm['mltfrz']).std(dim='time').values
    U_std['Q']=-noresm['hf_int'].std(dim='time').values*tera/rhow/cp/1e6 #1e6 for the Sverdrups
    S['FW']=0
    S['Q']=0
    T['FW']=0
    T['Q']=1

    return U,S,T,U_std

U,S,T,U_std=get_U_S_T_from_WM(WM)


(U['FW'])*1e3

U['Q']/(tera/rhow/cp/1e6)

U_mb,S_mb,T_mb,U_std_mb=get_U_S_T_from_WM(WM_mb)

wmvec=['PWS','AWS','DWS','PWN','AWN','FW','Q']

#get the full annual mean combos
Uyr={}
Ustdyr={}
Syr={}
Tyr={}
for yrS in ['2015','2016','2017']:
    for yrN in ['2006','2007','2008','2009','2010']:
        Uyr[yrS+'/'+yrN]={}
        Ustdyr[yrS+'/'+yrN]={}
        Syr[yrS+'/'+yrN]={}
        Tyr[yrS+'/'+yrN]={}
        for wm in ['AWN','PWN']:
            Uyr[yrS+'/'+yrN][wm]=WM['TRANS'].sel(WM=wm).sel(TIME=slice(yrN+'-1-1',yrN+'-12-31')).mean(dim='TIME').values
            Ustdyr[yrS+'/'+yrN][wm]=WM['TRANS'].sel(WM=wm).sel(TIME=slice(yrN+'-1-1',yrN+'-12-31')).std(dim='TIME').values
            Syr[yrS+'/'+yrN][wm]=WM['SA'].sel(WM=wm).sel(TIME=slice(yrN+'-1-1',yrN+'-12-31')).mean(dim='TIME').values
            Tyr[yrS+'/'+yrN][wm]=WM['CT'].sel(WM=wm).sel(TIME=slice(yrN+'-1-1',yrN+'-12-31')).mean(dim='TIME').values
        for wm in ['AWS','PWS','DWS']:
            Uyr[yrS+'/'+yrN][wm]=WM['TRANS'].sel(WM=wm).sel(TIME=slice(yrS+'-1-1',yrS+'-12-31')).mean(dim='TIME').values
            Ustdyr[yrS+'/'+yrN][wm]=WM['TRANS'].sel(WM=wm).sel(TIME=slice(yrS+'-1-1',yrS+'-12-31')).std(dim='TIME').values
            Syr[yrS+'/'+yrN][wm]=WM['SA'].sel(WM=wm).sel(TIME=slice(yrS+'-1-1',yrS+'-12-31')).mean(dim='TIME').values
            Tyr[yrS+'/'+yrN][wm]=WM['CT'].sel(WM=wm).sel(TIME=slice(yrS+'-1-1',yrS+'-12-31')).mean(dim='TIME').values
        for wm in ['FW','Q']:
            Uyr[yrS+'/'+yrN][wm]=U[wm]
            Ustdyr[yrS+'/'+yrN][wm]=U_std[wm]
            Syr[yrS+'/'+yrN][wm]=S[wm]
            Tyr[yrS+'/'+yrN][wm]=T[wm]


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

AM={}
x0={}

AM['base']=array([[1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],1]])

x0['base']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['Q']]

AM['massbal']=array([[1,1,1,0,0,0.5,0],\
[0,0,0,1,1,0.5,0],\
[S_mb['PWS'],S_mb['AWS'],S_mb['DWS'],S_mb['PWN'],S_mb['AWN'],S_mb['FW'],0],\
[T_mb['PWS'],T_mb['AWS'],T_mb['DWS'],T_mb['PWN'],T_mb['AWN'],T_mb['FW'],1]])


x0['massbal']=[U_mb['PWS'],U_mb['AWS'],U_mb['DWS'],U_mb['PWN'],U_mb['AWN'],U_mb['FW'],U_mb['Q']]

#vars that I want to be handy for later calcs
Snorm=35
Tnorm=5

####################################################################
#Calculate standard error for each unknown rather than just doing 20%
#Degrees of freedom.
#Using 4 degrees of freedom per year
pryr=4
ndegs={}
#OSNAP = 45 months
for kk in ['AWS','DWS','PWS']:
    ndegs[kk]=45/12*pryr
#NORTH = 68 months
for kk in ['AWN','PWN']:
    ndegs[kk]=68/12*pryr
# FW and Q = 19yrs
for kk in ['Q','FW']:
    ndegs[kk]=19*pryr

ndegs_yr=ndegs.copy()
for wm in wmvec[:5]:
    ndegs_yr[wm]=4


Winv=diag([1,1/Snorm,1/Tnorm])


Winv.dot(AM['base'].dot(x0['base']))

def run_inverse_model(zz,U,S,T,U_std,ndegs_var=ndegs):
    dv=-AM[zz].dot(x0[zz])

    Evec=[abs(U_std[wm])/sqrt(ndegs_var[wm]) for wm in wmvec]
    E=diag(Evec)

    if zz=='base':
        Winv=diag([1,1/Snorm,1/Tnorm])
    elif zz=='massbal':
        Winv=diag([1,1,1/Snorm,1/Tnorm])

    Umat,D,VmatT=linalg.svd(Winv.dot(AM[zz].dot(E)))

    Lambda_inv = zeros((AM[zz].shape[0], AM[zz].shape[1])).T
    Lambda_inv[:AM[zz].shape[0], :AM[zz].shape[0]] = diag(1/D)
    xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
    xsol_Ad=E.dot(xsol_prime)
    xbase=x0[zz]+xsol_Ad
    Pmat=E-E.dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E.dot(AM[zz].T))+linalg.inv(Winv)).dot(AM[zz].dot(E))))
    P=diag(Pmat)
    Ubase=get_U_from_x(xbase)
    Ue=get_U_from_x(Evec)
    Up=get_U_from_x(sqrt(P))
    return Ubase,Ue,xbase,Up

Ubase,Ue,xbase,Up=run_inverse_model('base',U,S,T,U_std)

Ue

Up

Ue['FW']*1e3
Up['FW']*1e3

Umb_sol,Umb_err,xmb,Upmb=run_inverse_model('massbal',U_mb,S_mb,T_mb,U_std_mb)

Usol_yr={}
Uerr_yr={}
xsol_yr={}
Up_yr={}
for yrS in ['2015','2016','2017']:
    for yrN in ['2006','2007','2008','2009','2010']:
        Usol_yr[yrS+'/'+yrN],Uerr_yr[yrS+'/'+yrN],xsol_yr[yrS+'/'+yrN],Up_yr[yrS+'/'+yrN]=run_inverse_model('base',Uyr[yrS+'/'+yrN],Syr[yrS+'/'+yrN],Tyr[yrS+'/'+yrN],Ustdyr[yrS+'/'+yrN])

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'#00F2F2','FW':'#00F2F2','Q':'limegreen'}

def plot_base_case_simple(Ubase,Ue,U,plt):
    f,axx=subplots(1,4,figsize=(9,2.5),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,1,1]))

    alf=0.75
    capi=7
    #U
    axx[0].bar(range(2),[U[kk] for kk in ['AWS','DWS']],color=[coldic[kk] for kk in ['AWS','DWS']],yerr=[Ue[kk] for kk in ['AWS','DWS']],capsize=capi,alpha=alf)
    axx[0].plot(range(2),[Ubase[kk] for kk in ['AWS','DWS']],'o',color='k')
    # for yrS in ['2015','2016','2017']:
    #     for yrN in ['2006','2007','2008','2009','2010']:
    #         axx[0].plot(arange(2)+0.2,[Usol_yr[yrS+'/'+yrN][kk] for kk in ['AWS','DWS']],'o',color='w',mec='k')


    ylimi=21
    axx[0].set_ylim(-ylimi,ylimi)
    ylimi=4.5
    axx[1].set_ylim(-ylimi,ylimi)
    axx[1].bar(range(3),[U[kk] for kk in ['PWS','PWN','AWN']],color=[coldic[kk] for kk in ['PWS','PWN','AWN']],yerr=[Ue[kk] for kk in ['PWS','PWN','AWN']],capsize=capi,alpha=alf)
    axx[1].plot(range(3),[Ubase[kk] for kk in ['PWS','PWN','AWN']],'o',color='k')
    # for yrS in ['2015','2016','2017']:
    #     for yrN in ['2006','2007','2008','2009','2010']:
    #         axx[1].plot(arange(3)+0.2,[Usol_yr[yrS+'/'+yrN][kk] for kk in ['PWS','PWN','AWN']],'o',color='w',mec='k')


    axx[2].bar(range(1),U['FW'],color=coldic['FW'],yerr=Ue['FW'],capsize=capi,alpha=alf)
    axx[2].plot(range(1),Ubase['FW'],'o',color='k')
    axx[2].set_ylim(-0.12,0.12)
    # for yrS in ['2015','2016','2017']:
    #     for yrN in ['2006','2007','2008','2009','2010']:
    #         axx[2].plot(0.2,Usol_yr[yrS+'/'+yrN]['FW'],'o',color='w',mec='k')


    fsz=14
    axx[0].set_ylabel('Volume transport [Sv]',fontsize=fsz)
    axx[3].set_ylabel('Heat flux [TW]',fontsize=fsz)
    axx[3].bar(0,cp*rhow*(U['Q'])/1e6,color=coldic['Q'],yerr=cp*rhow*Ue['Q']/1e6,capsize=capi,alpha=alf)
    axx[3].plot(0,cp*rhow*(Ubase['Q'])/1e6,'o',color='k')
    # for yrS in ['2015','2016','2017']:
    #     for yrN in ['2006','2007','2008','2009','2010']:
    #         axx[3].plot(0.2,cp*rhow*Usol_yr[yrS+'/'+yrN]['Q']/1e6,'o',color='w',mec='k')


    for ii in range(3):
        axx[ii].axhline(0,color='k')
    axx[0].set_xticks(range(2))
    axx[0].set_xticklabels(['AWS','DWS'])
    axx[1].set_xticks(range(3))
    axx[1].set_xticklabels(['PWS','PWN','AWN'])
    axx[2].set_xticks(range(1))
    axx[2].set_xticklabels(['SFW'])
    axx[3].set_xticks([0])
    axx[3].set_xticklabels('Q')

    savefig(figdir_paper+'_extra_2004/InvBudSol_'+plt+'.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/InvBudSol_'+plt+'.pdf',bbox_inches='tight')
    if plt=='base':
            savefig(figdir_paper+'/InvBudSol_'+plt+'.png',bbox_inches='tight')
            savefig(figdir_paper+'/InvBudSol_'+plt+'.pdf',bbox_inches='tight')

plot_base_case_simple(Ubase,Ue,U,'base')

T

U

Ue

Ue['Q']/(tera/rhow/cp/1e6)

Ubase
(Ubase['FW'])*1e3
Ubase['Q']/(tera/rhow/cp/1e6)

(Up['FW'])*1e3
Up['Q']/(tera/rhow/cp/1e6)

# plot_base_case_simple(Umb_sol,Umb_err,U_mb,'mb')
# [(  kk,Umb_sol[kk]-U_mb[kk]) for kk in Ubase]
##################################################################################
# Calculate fraction of fresh water vs. other water masses that goes into each limb
#################################################################################
#fraction of PWN in DWS limb
epsilon=arange(0,1.1,0.1)

def get_a_b_fracs(Ubase,S):
    #fraction of FW in PWS, as a function of epsilon
    a=((1-epsilon)*Ubase['PWN']*(S['PWN']/S['AWS']-1)+Ubase['PWS']*(S['PWS']/S['AWS']-1))/(Ubase['FW'])
    #fraction of FW in DWS, as a function of epsilon
    b=(epsilon*Ubase['PWN']*(S['PWN']/S['AWS']-1)+Ubase['DWS']*(S['DWS']/S['AWS']-1))/(Ubase['FW'])
    return a,b

a={}
b={}
a['base'],b['base']=get_a_b_fracs(Ubase,S)
a['mb'],b['mb']=get_a_b_fracs(Umb_sol,S_mb)

for yrS in ['2015','2016','2017']:
    for yrN in ['2006','2007','2008','2009','2010']:
        a[yrS+'/'+yrN],b[yrS+'/'+yrN]=get_a_b_fracs(Usol_yr[yrS+'/'+yrN],Syr[yrS+'/'+yrN])


def plot_adep():
    for ii,kk in enumerate(a):
        plot(1-epsilon,a[kk],linewidth=3,label=kk,color='C'+str(ii))

    xlabel('$\mathbf{1-\epsilon}$\nfraction of PWN in PWS')
    ylabel('$\mathbf{a}$\n fraction of SFW in PWS')
    xlim(0,1)
    axhline(0,color='k')
    legend(loc=(1.05,0))
    savefig(figdir_paper+'_extra_2004/FWfrac_mbdep.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/FWfrac_mbdep.pdf',bbox_inches='tight')


plot_adep()

#################################################################################
##### Look into how much Sea ice properties matter
#################################################################################
(-noresm['mltfrz']).mean(dim='time').values/U['FW']

x0_si=array([Ubase['PWS'],Ubase['AWS'],Ubase['DWS'],Ubase['PWN'],Ubase['AWN'],0.3*Ubase['FW'],0.7*Ubase['FW'],Ubase['Q']])

sivar={}
for S_SI in range(0,10,2):
    sivar[S_SI]={}
    for T_SI in range(-90,5,10):
        AM=array([[1,1,1,1,1,1,1,0],\
        [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S_SI,0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T_SI,1]])


        dv=-AM.dot(x0_si)

        Evec=array(hstack(([1]*5,x0_si[-3:]/5)))
        E=diag(Evec)
        Winv=diag([1,1/Snorm,1/Tnorm])
        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        sivar[S_SI][T_SI]=x0_si+xsol_Ad

def get_mats_from_dic(sivar):
    Svec=array([float(ff) for ff in sivar])
    Tvec=array([float(ff) for ff in sivar[Svec[0]]])
    simats={}
    for QQ,kk in enumerate(Ubase):
        simats[kk]=zeros((len(Svec),len(Tvec)))
        for ii,ss in enumerate(Svec):
            for jj,tt in enumerate(Tvec):
                simats[kk][ii,jj]=sivar[ss][tt][QQ]
    return Svec,Tvec,simats

Svec,Tvec,simats=get_mats_from_dic(sivar)

def plot_SIresponse():
    f,axx=subplots(2,4,figsize=(15,6),sharex=True,sharey=True)
    axivec=array([])
    for axirow in axx:
        for axi in axirow:
            axivec=hstack((axivec,axi))
    for axi,kk in zip(axivec,simats):
        if (kk=='FW') | (kk=='SI'):
            climi=10
            contit=axi.contourf(Svec,Tvec,(simats[kk].T-Ubase[kk])*1e3,vmin=-climi,vmax=climi,cmap=cm.RdBu)
            axi.set_title(kk+' [mSv]')
            cbar=colorbar(contit,ax=axi,format='%1.0f')
        elif kk=='Q':
            climi=30
            contit=axi.contourf(Svec,Tvec,cp*rhow*(simats['Q'].T-Ubase['Q'])/1e6,vmin=-climi,vmax=climi,cmap=cm.PiYG_r)
            axi.set_title(kk+' [TW]')
            cbar=colorbar(contit,ax=axi,format='%2.0f')
        else:
            climi=0.3
            contit=axi.contourf(Svec,Tvec,(simats[kk].T-Ubase[kk]),vmin=-climi,vmax=climi,cmap=cm.PuOr_r)
            axi.set_title(kk+' [Sv]')
            cbar=colorbar(contit,ax=axi,format='%0.2f')
        for label in cbar.ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)

    f.text(0.5, 0, 'sea ice salinity', ha='center',fontsize=14)
    f.text(0.05, 0.5, 'effective sea ice temperature [$^\circ$C]', va='center',rotation='vertical',fontsize=14)

    savefig(figdir_paper+'_extra_2004/SeaIce_paramdep.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/SeaIce_paramdep.pdf',bbox_inches='tight')

plot_SIresponse()


U_si=get_U_from_x(sivar[0][-30])
U_si_init=get_U_from_x(x0_si)

U_si['SI']=U_si['Q']
U_si['Q']=sivar[0][-30][-1]

U_si_init['SI']=U_si_init['Q']
U_si_init['Q']=x0_si[-1]

Udiff={}
for wm in U_si:
    Udiff[wm]=U_si[wm]-U_si_init[wm]
Udiff
Udiff['FW']*1e3
Udiff['SI']*1e3
Udiff['Q']/(tera/rhow/cp/1e6)

#################################################################################
##### Test dependence on PW salinity (both north and south)
#################################################################################

pwsvar={}
for S_PWNa in arange(-1,0.05,0.1):
    pwsvar[S_PWNa]={}
    for S_PWSa in arange(-1.0,0.05,0.1):
        AM=array([[1,1,1,1,1,1,0],\
        [S['PWS']+S_PWSa,S['AWS'],S['DWS'],S['PWN']+S_PWNa,S['AWN'],S['FW'],0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],1]])

        dv=-AM.dot(xbase)

        Evec=array(hstack(([1]*5,xbase[-2:]/5)))
        E=diag(Evec)
        Winv=diag([1,1/Snorm,1/Tnorm])
        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        pwsvar[S_PWNa][S_PWSa]=xbase+xsol_Ad

PWN_Svec,PWS_Svec,pwmats=get_mats_from_dic(pwsvar)
U_pw=get_U_from_x(pwsvar[-0.5000000000000001][-0.5000000000000001])

Udiff={}
for wm in Ubase:
    Udiff[wm]=U_pw[wm]-Ubase[wm]
Udiff
Udiff['FW']*1e3
Udiff['Q']/(tera/rhow/cp/1e6)
####################################################################################################
########  Response is pretty uniform: try to tease out a pattern (and look at other deps?)  #######
##################################################################################################
PWN_Smat,PWS_Smat=meshgrid(PWN_Svec,PWS_Svec)

sivar.keys()


U_pw=get_U_from_x(pwsvar[-0.5000000000000001][-0.5000000000000001])

def lineplot_PW_salinity():
    f,axx=subplots(1,3,figsize=(11,3),sharey=True)
    xind=-1
    yind=-1
    svr=len(PWS_Svec)
    xvar=[(S['PWN']+PWN_Smat)[xind,:],(S['PWS']+PWS_Smat)[:,yind],[(S['PWS']+PWS_Smat)[ii,ii] for ii in range(svr)]]
    ufw_tot=-Ubase['FW']
    yvar_fw=[pwmats['FW'].T[xind,:]+ufw_tot,pwmats['FW'].T[:,yind]+ufw_tot,array([pwmats['FW'].T[ii,ii]+ufw_tot for ii in range(svr)])]
    yvar_Q=[pwmats['Q'].T[xind,:]-Ubase['Q'],pwmats['Q'].T[:,yind]-Ubase['Q'],array([pwmats['Q'].T[ii,ii]-Ubase['Q'] for ii in range(svr)])]
    xlab=['PWN salinity','PWS salinity','PWS salinity']
    titvec=['a) Vary PWN salinity\n\nPWS = 34.4','b) Vary PWS salinity\n\nPWN = 33.7','c) Vary both PW salinities']
    lw=2
    for kk in ['AWS','PWS','DWS','AWN','PWN']:
            axx[0].plot(xvar[0],(pwmats[kk].T[xind,:]-Ubase[kk]),color=coldic[kk],label=kk,linewidth=lw)
            axx[1].plot(xvar[1],(pwmats[kk].T[:,yind]-Ubase[kk]),color=coldic[kk],label=kk,linewidth=lw)
            axx[2].plot(xvar[2],array([(pwmats[kk].T[ii,ii]-Ubase[kk])for ii in range(svr)]),color=coldic[kk],label=kk,linewidth=lw)
    for ii in range(3):
        ax1=axx[ii].twinx()
        for ll in ['']:
            ax1.plot(xvar[ii],(yvar_fw[ii])*1e3,color='c',linewidth=lw)
        ax2=axx[ii].twinx()
        ax2.plot(xvar[ii],cp*rhow*(yvar_Q[ii])/1e6,color='limegreen',linewidth=lw)
        axx[ii].set_xlabel(xlab[ii])
        ax1.set_ylim(-25,25)
        ax2.set_ylim(-60,60)
        axx[ii].set_title(titvec[ii],fontweight='bold')
        if ii!=2:
            ax1.set_yticklabels('')
            ax2.set_yticklabels('')
        axx[ii].set_xlim(xvar[ii][0],xvar[ii][-1])
    axx[0].set_ylim(-1.5,1.5)
    axx[0].set_yticks(arange(-1.5,1.6,0.5))
    ax2.spines["right"].set_position(("axes", 1.3))
    axx[0].set_ylabel('Transport anomaly [Sv]')
    ax1.set_ylabel('Fresh water flux anomaly [mSv]',color='c')
    ax2.set_ylabel('Heat flux anomaly [TW]',color='limegreen')
    ax1.tick_params(axis='y', colors='c')
    ax2.tick_params(axis='y', colors='limegreen')
    leg=axx[0].legend(loc=(0.5,-0.5),ncol=5,fontsize=13)
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    axi2=axx[2].twiny()
    axi2.set_xticks(arange(32.8,33.8,0.2))
    axi2.set_xlim(xvar[0][0],xvar[0][-1])
    axi2.set_xlabel('PWN salinity')
    # axx[2].axvline(34.4-0.5,color='k',zorder=0)
    # axx[0].set_title('a) Vary PWN salinities\n\n',fontweight='bold')
    # axx[1].set_title('b) Vary PWS salinities\n\n',fontweight='bold')
    # axx[2].set_title('c) Vary both PW salinities',fontweight='bold')
    savefig(figdir_paper+'/PWS_dep.png',bbox_inches='tight')
    savefig(figdir_paper+'/PWS_dep.pdf',bbox_inches='tight')

lineplot_PW_salinity()

#######################################################################################
##############  What happens if we add more FW? (Like 100mSv)  ###########################
#######################################################################################

fwvar={}
for U_FW in arange(0,0.11,0.01):
        AM=array([[1,1,1,1,1,1,0],\
        [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],1]])

        xinit=xbase.copy()
        xinit[5]=xinit[5]+U_FW
        dv=-AM.dot(xinit)

        Evec=[abs(U_std[wm]/sqrt(ndegs[wm])) for wm in wmvec]
        Evec[5]=1e-10
        E=diag(Evec)
        Winv=diag([1,1/Snorm,1/Tnorm])
        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        fwvar[U_FW]=xinit+xsol_Ad

Ubase['FW']*1e3

U_fwvar=get_U_from_x(fwvar[0.02])


Udiff={}
for wm in Ubase:
    Udiff[wm]=U_fwvar[wm]-Ubase[wm]
Udiff
Udiff['FW']*1e3
Udiff['Q']/(tera/rhow/cp/1e6)

a_fw,b_fw=get_a_b_fracs(U_fwvar,S)

U['FW']
Ubase['FW']
U_fwvar['FW']

U_fwvar['Q']*cp*rhow/1e6

#######################################################################################
##############  Now look at consequences  for FW dist  ###########################
#######################################################################################
a_pwmat=zeros((len(epsilon),shape(pwmats['Q'])[1],shape(pwmats['Q'])[0]))
b_pwmat=a_pwmat.copy()
for ii,ee in enumerate(1-epsilon):
    a_pwmat[ii,:,:]=(ee*pwmats['PWN'].T*((S['PWN']+PWN_Smat)/S['AWS']-1)+pwmats['PWS'].T*((S['PWS']+PWS_Smat)/S['AWS']-1))/(pwmats['FW'].T)
    b_pwmat[ii,:,:]=((1-ee)*pwmats['PWN'].T*((S['PWN']+PWN_Smat)/S['AWS']-1)+pwmats['DWS'].T*(S['DWS']/S['AWS']-1))/(pwmats['FW'].T)
c_pwmat=1-a_pwmat-b_pwmat

############################### Calculate out relevant transports/partitions

# From NorESM spatial patterns, 47% is in west, 53% in east
FWfrac_west=0.47
FWfrac_east=0.53

def get_split(U_choose,a_choose,b_choose):
    split={}
    split['FW']={}
    split['PWN']={}
    split['AWS']={}
    split['Q']={}
    # get amount of FW in each limb
    split['FW']['Est']=FWfrac_west*(U_choose['FW'])*1e3
    split['FW']['AW']=(1-a_choose[0]-b_choose[0])*(U_choose['FW'])*1e3 #doesn't matter which a and b you choose (as long as same index), will all have same sum.
    split['FW']['Ov']=FWfrac_east*(U_choose['FW'])*1e3-split['FW']['AW']
    # given this distribution, how much PWN flows into PWS?
    split['PWN']['PWS']=U_choose['PWN']*(split['FW']['Est']-(U_choose['FW'])*1e3*a_choose[-1])/(U_choose['FW'])/1e3/(a_choose[0]-a_choose[-1])
    #and how much into DWS?
    split['PWN']['DWS']=U_choose['PWN']-split['PWN']['PWS']
    # How is AWS split between the limbs?
    split['AWS']['PWS']=-(U_choose['PWS']*S['PWS']+split['PWN']['PWS']*S['PWN'])/S['AWS']
    split['AWS']['DWS']=-(U_choose['DWS']*S['DWS']+split['PWN']['DWS']*S['PWN'])/S['AWS']
    split['AWS']['AWN']=U_choose['AWS']-split['AWS']['PWS']-split['AWS']['DWS']
    # What about heat flux?
    split['Q']['Est']=-(U_choose['PWS']*T['PWS']+split['PWN']['PWS']*T['PWN']+split['AWS']['PWS']*T['AWS'])*cp*rhow/1e6
    split['Q']['Ov']=-(U_choose['DWS']*T['DWS']+split['PWN']['DWS']*T['PWN']+split['AWS']['DWS']*T['AWS'])*cp*rhow/1e6
    split['Q']['AW']=-(U_choose['AWN']*T['AWN']+split['AWS']['AWN']*T['AWS'])*cp*rhow/1e6
    split['Q']['check']=split['Q']['Est']+split['Q']['Ov']+split['Q']['AW']-U_choose['Q']*cp*rhow/1e6
    return split

split_base=get_split(Ubase,a_pwmat[:,10,10],b_pwmat[:,10,10])

split_base

split_pw=get_split(U_pw,a_pwmat[:,5,5],b_pwmat[:,5,5])
split_pw['Q']
U_std

split_fw=get_split(U_fwvar,a_fw,b_fw)

split_fw
Ubase['FW']*1e3-split_base['FW']['AW']
split_base
Ubase


split_yr={}
for yrS in ['2015','2016','2017']:
    for yrN in ['2006','2007','2008','2009','2010']:
        split_yr[yrS+'/'+yrN]=get_split(Usol_yr[yrS+'/'+yrN],a[yrS+'/'+yrN],b[yrS+'/'+yrN])

#The reason for the branching into 3 is the difference in AWS salinity in the three years
#Big impact on how much need to dilute AWS to get AWN, which sets overall FW totals
# ii=0
# for yrS in ['2015','2016','2017']:
#     for yrN in ['2006','2007','2008','2009','2010']:
#         plot(ii,Syr[yrS+'/'+yrN]['AWS'],'o')
#         ii+=1


########### Plot the dependence!
basecol='k'
pwcol='#FFBC0A'
fwcol='#04724d'
ash='s'

def plot_adep_base():
    msize=8
    f,axx=subplots(1,2,figsize=(10,3.5),sharex=True)
    f.subplots_adjust(wspace=0.4)
    alf=0.3
    for ii,var in enumerate([a_pwmat,b_pwmat]):
        if ii==0:
            xvar=(1-epsilon)
            xvar2=1
            xvar3=0
        else:
            xvar=epsilon
            xvar2=0
            xvar3=1
            # for yrS in ['2015','2016','2017']:
            #     for yrN in ['2006','2007','2008','2009','2010']:
            #         axx[ii].plot(xvar*Usol_yr[yrS+'/'+yrN]['PWN'],b[yrS+'/'+yrN]*(Usol_yr[yrS+'/'+yrN]['FW']+Usol_yr[yrS+'/'+yrN]['SI'])*1e3,linewidth=1,color=basecol,label='',zorder=2,alpha=alf)
        axx[ii].plot(xvar*Ubase['PWN'],var[:,10,10]*(Ubase['FW'])*1e3,linewidth=4,color=basecol,label='Full mean solution',zorder=5)
        axx[ii].plot(xvar2*Ubase['PWN'],var[0,10,10]*(Ubase['FW'])*1e3,'o',color='w',label='',zorder=5,mec='k',markersize=msize)
        axx[ii].plot(xvar3*Ubase['PWN'],var[-1,10,10]*(Ubase['FW'])*1e3,ash,color='w',label='',zorder=5,mec='k',markersize=msize)
        axx[ii].set_ylim(0,100)
        axx[ii].set_xlim(-0.05,1.3)
        # if ii==0:
            # for yrS in ['2015','2016','2017']:
            #     for yrN in ['2006','2007','2008','2009','2010']:
            #         if (yrS=='2015') & (yrN=='2006'):
            #             axx[ii].plot(xvar*Usol_yr[yrS+'/'+yrN]['PWN'],a[yrS+'/'+yrN]*(Usol_yr[yrS+'/'+yrN]['FW']+Usol_yr[yrS+'/'+yrN]['SI'])*1e3,linewidth=1,color=basecol,label='Combinations of single year means',zorder=2,alpha=alf)
            #         else:
            #             axx[ii].plot(xvar*Usol_yr[yrS+'/'+yrN]['PWN'],a[yrS+'/'+yrN]*(Usol_yr[yrS+'/'+yrN]['FW']+Usol_yr[yrS+'/'+yrN]['SI'])*1e3,linewidth=1,color=basecol,label='',zorder=2,alpha=alf)
    axx[0].plot(split_base['PWN']['PWS'],split_base['FW']['Est'],'*',color='w',label='',markersize=15,zorder=100,mec='k')
    axx[1].plot(split_base['PWN']['DWS'],split_base['FW']['Ov'],'*',color='w',label='',markersize=15,zorder=100,mec='k')
    # axx[0].legend(loc=(0.2,-0.5),ncol=3,fontsize=12)
    axx[0].set_title('a) Estuarine limb',fontsize=14)
    axx[1].set_title('b) Overturning limb',fontsize=14)
    axx[0].set_ylabel('$\mathbf{\delta}\ U_{FW}$\nSFW transport in $\mathbf{PWS}$ [mSv]')
    axx[1].set_ylabel('$\mathbf{\gamma}\ U_{FW}$\nSFW transport in $\mathbf{DWS}$ [mSv]')
    axx[0].set_xlabel('$\mathbf{(1-\epsilon)} \ U_{PWN}$\nPWN transport in $\mathbf{PWS}$ [Sv]')
    axx[1].set_xlabel('$\mathbf{\epsilon} \ U_{PWN}$\nPWN transport in $\mathbf{DWS}$ [Sv]')
    # axx[0].axhline(split_base['FW']['Est'],color='k',linestyle='--',zorder=1)
    # axx[1].axhline(split_base['FW']['Ov'],color='k',linestyle='--',zorder=1)
    savefig(figdir_paper+'/FWfrac_obs.png',bbox_inches='tight')
    savefig(figdir_paper+'/FWfrac_obs.pdf',bbox_inches='tight')

plot_adep_base()

a_pwmat[0,10,10]*(Ubase['FW'])*1e3,
b_pwmat[0,10,10]*(Ubase['FW'])*1e3,

a_pwmat[-1,10,10]*(Ubase['FW'])*1e3,
b_pwmat[-1,10,10]*(Ubase['FW'])*1e3,

basecol='#8b1c62'

def plot_adep_pw():
    msize=8
    f,axx=subplots(1,2,figsize=(11,3.5),sharex=True)
    f.subplots_adjust(wspace=0.3)
    for ii,var in enumerate([a_pwmat,b_pwmat]):
        if ii==0:
            xvar=(1-epsilon)
            xvar2=1
            xvar3=0
            # for yrS in ['2015','2016','2017']:
            #     for yrN in ['2006','2007','2008','2009','2010']:
            #         axx[ii].plot(xvar*Usol_yr[yrS+'/'+yrN]['PWN'],a[yrS+'/'+yrN]*(Usol_yr[yrS+'/'+yrN]['FW']+Usol_yr[yrS+'/'+yrN]['SI'])*1e3,linewidth=1,color=basecol,label='',zorder=2,alpha=0.3)
        else:
            xvar=epsilon
            xvar2=0
            xvar3=1
            # for yrS in ['2015','2016','2017']:
            #     for yrN in ['2006','2007','2008','2009','2010']:
            #         axx[ii].plot(xvar*Usol_yr[yrS+'/'+yrN]['PWN'],b[yrS+'/'+yrN]*(Usol_yr[yrS+'/'+yrN]['FW']+Usol_yr[yrS+'/'+yrN]['SI'])*1e3,linewidth=1,color=basecol,label='',zorder=2,alpha=0.3)
        axx[ii].plot(xvar*Ubase['PWN'],var[:,10,10]*(Ubase['FW'])*1e3,linewidth=4,color=basecol,label='Base case',zorder=5)
        axx[ii].plot(xvar*U_pw['PWN'],var[:,5,5]*(U_pw['FW'])*1e3,color=pwcol,zorder=4,label='Polar Waters fresher by 0.5',linewidth=4)
        axx[ii].plot(xvar2*Ubase['PWN'],var[0,10,10]*(Ubase['FW'])*1e3,'o',color=basecol,label='',zorder=5,markersize=msize,mec='k')
        axx[ii].plot(xvar2*U_pw['PWN'],var[0,5,5]*(U_pw['FW'])*1e3,'o',color=pwcol,zorder=4,label='',mec='k',markersize=msize)
        axx[ii].plot(xvar3*Ubase['PWN'],var[-1,10,10]*(Ubase['FW'])*1e3,ash,color=basecol,label='',zorder=5,markersize=msize,mec='k')
        axx[ii].plot(xvar3*U_pw['PWN'],var[-1,5,5]*(U_pw['FW'])*1e3,ash,color=pwcol,zorder=4,label='',mec='k',markersize=msize)
        axx[ii].set_ylim(-10,140)
    axx[0].plot((1-epsilon)*U_fwvar['PWN'],a_fw*(U_fwvar['FW'])*1e3,linewidth=4,color=fwcol,label='Add 20 mSv of Surface Fresh Water')
    axx[1].plot(epsilon*U_fwvar['PWN'],b_fw*(U_fwvar['FW'])*1e3,linewidth=4,color=fwcol)
    axx[0].plot(U_fwvar['PWN'],a_fw[0]*(U_fwvar['FW'])*1e3,'o',color=fwcol,label='',mec='k',markersize=msize)
    axx[1].plot(0,b_fw[0]*(U_fwvar['FW'])*1e3,'o',color=fwcol,label='',mec='k',markersize=msize)
    axx[0].plot(0,a_fw[-1]*(U_fwvar['FW'])*1e3,ash,color=fwcol,label='',mec='k',markersize=msize)
    axx[1].plot(U_fwvar['PWN'],b_fw[-1]*(U_fwvar['FW'])*1e3,ash,color=fwcol,label='',mec='k',markersize=msize)
    axx[0].plot(split_base['PWN']['PWS'],split_base['FW']['Est'],'*',color=basecol,label='',markersize=15,zorder=100,mec='k')
    axx[1].plot(split_base['PWN']['DWS'],split_base['FW']['Ov'],'*',color=basecol,label='',markersize=15,zorder=100,mec='k')
    axx[0].plot(split_pw['PWN']['PWS'],split_pw['FW']['Est'],'*',color=pwcol,label='',markersize=15,zorder=100,mec='k')
    axx[1].plot(split_pw['PWN']['DWS'],split_pw['FW']['Ov'],'*',color=pwcol,label='',markersize=15,zorder=100,mec='k')
    axx[0].plot(split_fw['PWN']['PWS'],split_fw['FW']['Est'],'*',color=fwcol,label='',markersize=15,zorder=100,mec='k')
    axx[1].plot(split_fw['PWN']['DWS'],split_fw['FW']['Ov'],'*',color=fwcol,label='',markersize=15,zorder=100,mec='k')
    axx[0].legend(loc=(0,-0.5),ncol=3,fontsize=12)
    axx[0].set_title('a) Estuarine limb',fontsize=14)
    axx[1].set_title('b) Overturning limb',fontsize=14)
    axx[0].set_ylabel('$\mathbf{\delta}\ U_{FW}$\nSFW transport in $\mathbf{PWS}$ [mSv]')
    axx[1].set_ylabel('$\mathbf{\gamma}\ U_{FW}$\nSFW transport in $\mathbf{DWS}$ [mSv]')
    axx[0].set_xlabel('$\mathbf{(1-\epsilon)} \ U_{PWN}$\nPWN transport in $\mathbf{PWS}$ [Sv]')
    axx[1].set_xlabel('$\mathbf{\epsilon} \ U_{PWN}$\nPWN transport in $\mathbf{DWS}$ [Sv]')
    for axi in axx[0],axx[1]:
        axi.axhline(0,color='k')
        axi.set_xlim(-0.05,2)
    # axx[0].axhline(split_base['FW']['Est'],color='k',linestyle='--',zorder=50)
    # axx[1].axhline(split_base['FW']['Ov'],color='k',linestyle='--',zorder=50)
    savefig(figdir_paper+'/FWfrac_obs_allsens.png',bbox_inches='tight')
    savefig(figdir_paper+'/FWfrac_obs_allsens.pdf',bbox_inches='tight')

plot_adep_pw()


# #######################################################################################
# ##############  GRAVEYARD  ###########################
# #######################################################################################
#
#
#
#
# ##############  What happens if we add more FW and make PWS fresher?  ###########################
#
# AM=array([[1,1,1,1,1,1,1,0],\
# [S['PWS']-0.5,S['AWS'],S['DWS'],S['PWN']-0.5,S['AWN'],S['FW'],S['SI'],0],\
# [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])
#
# xinit=xbase.copy()
# xinit[5]=xinit[5]+0.02
# dv=-AM.dot(xinit)
#
# Evec=xinit/5
# Evec[5:7]=1e-10
# E=diag(Evec)
# Winv=diag([1,1/Snorm,1/Tnorm])
# Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))
#
# Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
# Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
# xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
# xsol_Ad=E.dot(xsol_prime)
# x_both=xinit+xsol_Ad
# U_both=get_U_from_x(x_both)
#
# S_PW=S.copy()
# S_PW['PWS']=S['PWS']-0.5
# S_PW['PWN']=S['PWN']-0.5
# a_both,b_both=get_a_b_fracs(U_both,S_PW)

# def plot_in_each(axi):
#     axi.plot(S['PWN'],S['PWS'],'ko',markersize=10)
#     axi.plot(S['PWN']+PWN_Svec,S['PWN']+PWN_Svec,'r-',linewidth=3)
#
# def plot_PW_Sdep(Svec,Tvec,simats):
#     f,axx=subplots(2,4,figsize=(15,6),sharex=True,sharey=True)
#     axivec=array([])
#     for axirow in axx:
#         for axi in axirow:
#             axivec=hstack((axivec,axi))
#     for axi,kk in zip(axivec,simats):
#         if (kk=='FW') | (kk=='SI'):
#             climi=20
#             contit=axi.contourf(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk])*1e3,vmin=-climi,vmax=climi,cmap=cm.RdBu)
#             axi.contour(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),levels=[0],colors='k')
#             axi.set_title(kk+' [mSv]')
#             cbar=colorbar(contit,ax=axi,format='%1.0f')
#             plot_in_each(axi)
#         elif kk=='Q':
#             climi=30
#             contit=axi.contourf(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,cp*rhow*(pwmats['Q'].T-Ubase['Q'])/1e6,vmin=-climi,vmax=climi,cmap=cm.PiYG_r)
#             axi.contour(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),levels=[0],colors='k')
#             axi.set_title(kk+' [TW]')
#             cbar=colorbar(contit,ax=axi,format='%2.0f')
#             plot_in_each(axi)
#         else:
#             climi=1.5
#             contit=axi.contourf(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),vmin=-climi,vmax=climi,cmap=cm.PuOr_r)
#             axi.contour(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),levels=[0],colors='k')
#             axi.set_title(kk+' [Sv]')
#             cbar=colorbar(contit,ax=axi,format='%0.2f')
#             plot_in_each(axi)
#         for label in cbar.ax.yaxis.get_ticklabels()[1::2]:
#             label.set_visible(False)
#     axi.set_ylim(S['PWS']+PWS_Svec[0],S['PWS']+PWS_Svec[-1])
#     f.text(0.5, 0, 'PWN salinity', ha='center',fontsize=14)
#     f.text(0.05, 0.5, 'PWS salinity', va='center',rotation='vertical',fontsize=14)
#
#     savefig(figdir_paper+'_extra_2004/PW_Sdep.png',bbox_inches='tight')
#     savefig(figdir_paper+'_extra_2004/PW_Sdep.pdf',bbox_inches='tight')
#
#
# plot_PW_Sdep(PWN_Svec,PWS_Svec,pwmats)


# def plot_PW_Sdep_lines():
#     f,axx=subplots(2,4,figsize=(15,6),sharex=True)
#     axivec=array([])
#     for axirow in axx:
#         for axi in axirow:
#             axivec=hstack((axivec,axi))
#     for axi,kk in zip(axivec,simats):
#             axi.plot(((S['PWN']+PWN_Smat)-(S['PWS']+PWS_Smat))[-2,:],(pwmats[kk].T[-2,:]),label='vary PWN salinity')
#             axi.plot(((S['PWN']+PWN_Smat)-(S['PWS']+PWS_Smat))[:,-3],(pwmats[kk].T[:,-3]),label='vary PWS salinity')
#             axi.plot(((S['PWN'])-(S['PWS'])),(Ubase[kk]),'ko',label='base case')
#             axi.plot(((S['PWN'])-(S['PWS'])),(pwmats[kk].T[5,5]),'ro',label='both 0.5 fresher')
#             axi.plot(((S['PWN'])-(S['PWS'])),(pwmats[kk].T[0,0]),'go',label='both 1 fresher')
#             axi.set_title(kk)
#     axi.legend(loc=(1,0.7))
#     f.text(0.5, 0, 'PWN salinity - PWS salinity', ha='center',fontsize=14)
#     # f.text(0.05, 0.5, 'PWS salinity', va='center',rotation='vertical',fontsize=14)
#
#     # savefig(figdir_paper+'/PW_Sdep.png',bbox_inches='tight')
#     # savefig(figdir_paper+'/PW_Sdep.pdf',bbox_inches='tight')
#
# plot_PW_Sdep_lines()
# Ubase.keys()
