from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Linear/'
figdir_paper='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/paperfigs'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################
WM=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_2004.nc')
WM_mb=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_mb_2004.nc')

cp=3850
rhow=1025
tera=10**12
#Noresm (taking sea ice into account)
Q=-251*tera/rhow/cp/1e6 #for the Sverdrups

def get_U_S_T_from_WM(WM):
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

U,S,T=get_U_S_T_from_WM(WM)
U_mb,S_mb,T_mb=get_U_S_T_from_WM(WM_mb)

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
x0={}

AM['base']=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

x0['base']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],U['Q']]

AM['massbal']=array([[1,1,1,0,0,0.5,0.5,0],\
[0,0,0,1,1,0.5,0.5,0],\
[S_mb['PWS'],S_mb['AWS'],S_mb['DWS'],S_mb['PWN'],S_mb['AWN'],S_mb['FW'],S_mb['SI'],0],\
[T_mb['PWS'],T_mb['AWS'],T_mb['DWS'],T_mb['PWN'],T_mb['AWN'],T_mb['FW'],T_mb['SI'],1]])

x0['massbal']=[U_mb['PWS'],U_mb['AWS'],U_mb['DWS'],U_mb['PWN'],U_mb['AWN'],U_mb['FW'],U_mb['SI'],U_mb['Q']]

#vars that I want to be handy for later calcs
Qvar=10
Snorm=35
Tnorm=5
def run_inverse_model(zz,U,S,T):
    dv=-AM[zz].dot(x0[zz])

    if zz=='base':
        Winv=diag([1,1/Snorm,1/Tnorm])
    elif zz=='massbal':
        Winv=diag([1,1,1/Snorm,1/Tnorm])


    Evec=hstack((abs(array([U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN']])/5),0.02,0.02,Qvar))
    # Evec=hstack((5*[1],0.02,0.02,Qvar))
    E=diag(Evec)

    Umat,D,VmatT=linalg.svd(Winv.dot(AM[zz].dot(E)))

    Lambda_inv = zeros((AM[zz].shape[0], AM[zz].shape[1])).T
    Lambda_inv[:AM[zz].shape[0], :AM[zz].shape[0]] = diag(1/D)
    xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
    xsol_Ad=E.dot(xsol_prime)
    xbase=x0[zz]+xsol_Ad
    P=diag(E-E.dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E.dot(AM[zz].T))+linalg.inv(Winv)).dot(AM[zz].dot(E)))))
    Ubase=get_U_from_x(xbase)
    Ue=get_U_from_x(P)
    return Ubase,Ue,xbase

Ubase,Ue,xbase=run_inverse_model('base',U,S,T)

Umb_sol,Umb_err,xmb=run_inverse_model('massbal',U_mb,S_mb,T_mb)

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'cyan','Q':'limegreen'}

def plot_base_case_simple(Ubase,Ue,plt):
    f,axx=subplots(1,4,figsize=(9,2.5),constrained_layout=True,gridspec_kw=dict(width_ratios=[2,3,1,1]))

    alf=0.75
    capi=7
    #U
    axx[0].bar(range(2),[Ubase[kk] for kk in ['AWS','DWS']],color=[coldic[kk] for kk in ['AWS','DWS']],yerr=[Ue[kk] for kk in ['AWS','DWS']],capsize=capi,alpha=alf)
    axx[0].plot(range(2),[U[kk] for kk in ['AWS','DWS']],'o',color='k')

    ylimi=20
    axx[0].set_ylim(-ylimi,ylimi)
    ylimi=4
    axx[1].set_ylim(-ylimi,ylimi)
    axx[1].bar(range(3),[Ubase[kk] for kk in ['PWS','PWN','AWN']],color=[coldic[kk] for kk in ['PWS','PWN','AWN']],yerr=[Ue[kk] for kk in ['PWS','PWN','AWN']],capsize=capi,alpha=alf)
    axx[1].plot(range(3),[U[kk] for kk in ['PWS','PWN','AWN']],'o',color='k')

    axx[2].bar(range(1),U['SI']+Ubase['FW'],color=coldic['FW'],yerr=Ue['SI']+Ue['FW'],capsize=capi,alpha=alf)
    axx[2].plot(range(1),U['SI']+U['FW'],'o',color='k')
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

    savefig(figdir_paper+'_extra_2004/InvBudSol_'+plt+'.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/InvBudSol_'+plt+'.pdf',bbox_inches='tight')

plot_base_case_simple(Ubase,Ue,'base')

WM.sel(WM='PWS').PSAL.plot()

basediff=[(kk,Ubase[kk]-U[kk]) for kk in Ubase]
basediff

plot_base_case_simple(Umb_sol,Umb_err,'mb')
[(kk,Umb_sol[kk]-U_mb[kk]) for kk in Ubase]
##################################################################################
# Calculate fraction of fresh water vs. other water masses that goes into each limb
#################################################################################
#fraction of PWN in DWS limb
epsilon=arange(0,1.1,0.1)

def get_a_b_fracs(Ubase,S):
    #fraction of FW in PWS, as a function of epsilon
    a=((1-epsilon)*Ubase['PWN']*(S['PWN']/S['AWS']-1)+Ubase['PWS']*(S['PWS']/S['AWS']-1))/(Ubase['FW']+Ubase['SI'])
    #fraction of FW in DWS, as a function of epsilon
    b=(epsilon*Ubase['PWN']*(S['PWN']/S['AWS']-1)+Ubase['DWS']*(S['DWS']/S['AWS']-1))/(Ubase['FW']+Ubase['SI'])
    return a,b



S['PWN']/S['AWS']
S['PWS']/S['AWS']
S['DWS']/S['AWS']

Ubase['PWS']
Ubase['DWS']
Ubase['PWN']*(S['PWN']/S['AWS']-1)
Ubase['PWS']*(S['PWS']/S['AWS']-1)
Ubase['DWS']*(S['DWS']/S['AWS']-1)

(Ubase['FW']+Ubase['SI'])

a={}
b={}
a['base'],b['base']=get_a_b_fracs(Ubase,S)
a['mb'],b['mb']=get_a_b_fracs(Umb_sol,S_mb)
[(kk,S[kk]-S_mb[kk]) for kk in S]
def plot_adep():
    for ii,kk in enumerate(a):
        plot(1-epsilon,a[kk],linewidth=3,label=kk,color='C'+str(ii))

    xlabel('$\mathbf{1-\epsilon}$\nfraction of PWN in PWS')
    ylabel('$\mathbf{a}$\n fraction of (FW + SI) in PWS')
    xlim(0,1)
    axhline(0,color='k')
    legend()
    savefig(figdir_paper+'_extra_2004/FWfrac_mbdep.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/FWfrac_mbdep.pdf',bbox_inches='tight')


plot_adep()

#################################################################################
##### Look into how much Sea ice properties matter
#################################################################################
sivar={}
for S_SI in range(0,10,2):
    sivar[S_SI]={}
    for T_SI in range(-90,5,10):
        AM=array([[1,1,1,1,1,1,1,0],\
        [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S_SI,0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T_SI,1]])

        dv=-AM.dot(xbase)

        Evec=hstack((5*[1],0.02,0.02,Qvar))
        E=diag(Evec)
        Winv=diag([1,1/Snorm,1/Tnorm])
        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))


        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        sivar[S_SI][T_SI]=xbase+xsol_Ad

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

contourf(simats['AWN'].T-Ubase['AWN']+simats['PWN'].T-Ubase['PWN'])
colorbar()

#################################################################################
##### Test dependence on PW salinity (both north and south)
#################################################################################

pwsvar={}
for S_PWNa in arange(-1,0.3,0.1):
    pwsvar[S_PWNa]={}
    for S_PWSa in arange(-1.0,0.2,0.1):
        AM=array([[1,1,1,1,1,1,1,0],\
        [S['PWS']+S_PWSa,S['AWS'],S['DWS'],S['PWN']+S_PWNa,S['AWN'],S['FW'],S['SI'],0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

        dv=-AM.dot(xbase)

        Evec=hstack((5*[1],0.02,0.02,Qvar))
        E=diag(Evec)
        Winv=diag([1,1/Snorm,1/Tnorm])
        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        pwsvar[S_PWNa][S_PWSa]=xbase+xsol_Ad

PWN_Svec,PWS_Svec,pwmats=get_mats_from_dic(pwsvar)

def plot_in_each(axi):
    axi.plot(S['PWN'],S['PWS'],'ko',markersize=10)
    axi.plot(S['PWN']+PWN_Svec,S['PWN']+PWN_Svec,'r-',linewidth=3)

def plot_PW_Sdep(Svec,Tvec,simats):
    f,axx=subplots(2,4,figsize=(15,6),sharex=True,sharey=True)
    axivec=array([])
    for axirow in axx:
        for axi in axirow:
            axivec=hstack((axivec,axi))
    for axi,kk in zip(axivec,simats):
        if (kk=='FW') | (kk=='SI'):
            climi=20
            contit=axi.contourf(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk])*1e3,vmin=-climi,vmax=climi,cmap=cm.RdBu)
            axi.contour(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),levels=[0],colors='k')
            axi.set_title(kk+' [mSv]')
            cbar=colorbar(contit,ax=axi,format='%1.0f')
            plot_in_each(axi)
        elif kk=='Q':
            climi=30
            contit=axi.contourf(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,cp*rhow*(pwmats['Q'].T-Ubase['Q'])/1e6,vmin=-climi,vmax=climi,cmap=cm.PiYG_r)
            axi.contour(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),levels=[0],colors='k')
            axi.set_title(kk+' [TW]')
            cbar=colorbar(contit,ax=axi,format='%2.0f')
            plot_in_each(axi)
        else:
            climi=1.5
            contit=axi.contourf(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),vmin=-climi,vmax=climi,cmap=cm.PuOr_r)
            axi.contour(S['PWN']+PWN_Svec,S['PWS']+PWS_Svec,(pwmats[kk].T-Ubase[kk]),levels=[0],colors='k')
            axi.set_title(kk+' [Sv]')
            cbar=colorbar(contit,ax=axi,format='%0.2f')
            plot_in_each(axi)
        for label in cbar.ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    axi.set_ylim(S['PWS']+PWS_Svec[0],S['PWS']+PWS_Svec[-1])
    f.text(0.5, 0, 'PWN salinity', ha='center',fontsize=14)
    f.text(0.05, 0.5, 'PWS salinity', va='center',rotation='vertical',fontsize=14)

    savefig(figdir_paper+'_extra_2004/PW_Sdep.png',bbox_inches='tight')
    savefig(figdir_paper+'_extra_2004/PW_Sdep.pdf',bbox_inches='tight')


plot_PW_Sdep(PWN_Svec,PWS_Svec,pwmats)


####################################################################################################
########  Response is pretty uniform: try to tease out a pattern (and look at other deps?)  #######
##################################################################################################
PWN_Smat,PWS_Smat=meshgrid(PWN_Svec,PWS_Svec)

def plot_PW_Sdep_lines():
    f,axx=subplots(2,4,figsize=(15,6),sharex=True)
    axivec=array([])
    for axirow in axx:
        for axi in axirow:
            axivec=hstack((axivec,axi))
    for axi,kk in zip(axivec,simats):
            axi.plot(((S['PWN']+PWN_Smat)-(S['PWS']+PWS_Smat))[-2,:],(pwmats[kk].T[-2,:]),label='vary PWN salinity')
            axi.plot(((S['PWN']+PWN_Smat)-(S['PWS']+PWS_Smat))[:,-3],(pwmats[kk].T[:,-3]),label='vary PWS salinity')
            axi.plot(((S['PWN'])-(S['PWS'])),(Ubase[kk]),'ko',label='base case')
            axi.plot(((S['PWN'])-(S['PWS'])),(pwmats[kk].T[5,5]),'ro',label='both 0.5 fresher')
            axi.plot(((S['PWN'])-(S['PWS'])),(pwmats[kk].T[0,0]),'go',label='both 1 fresher')
            axi.set_title(kk)
    axi.legend(loc=(1,0.7))
    f.text(0.5, 0, 'PWN salinity - PWS salinity', ha='center',fontsize=14)
    # f.text(0.05, 0.5, 'PWS salinity', va='center',rotation='vertical',fontsize=14)

    # savefig(figdir_paper+'/PW_Sdep.png',bbox_inches='tight')
    # savefig(figdir_paper+'/PW_Sdep.pdf',bbox_inches='tight')

plot_PW_Sdep_lines()
Ubase.keys()


U_si=get_U_from_x(sivar[6][-30])
U_pw=get_U_from_x(pwsvar[-1][-1])
U_si
U_pw
U_pw['Q']*cp*rhow/1e6

def lineplot_PW_salinity():
    f,axx=subplots(1,3,figsize=(13.5,3.5),sharey=True)
    xind=-2
    yind=-3
    svr=12
    xvar=[(S['PWN']+PWN_Smat)[xind,:],(S['PWS']+PWS_Smat)[:,yind],[(S['PWS']+PWS_Smat)[ii,ii] for ii in range(svr)]]
    ufw_tot=-Ubase['SI']-Ubase['FW']
    yvar_fw=[pwmats['FW'].T[xind,:]+pwmats['SI'].T[xind,:]+ufw_tot,pwmats['FW'].T[:,yind]+pwmats['SI'].T[:,yind]+ufw_tot,array([pwmats['FW'].T[ii,ii]+pwmats['SI'].T[ii,ii]+ufw_tot for ii in range(svr)])]
    yvar_Q=[pwmats['Q'].T[xind,:]-Ubase['Q'],pwmats['Q'].T[:,yind]-Ubase['Q'],array([pwmats['Q'].T[ii,ii]-Ubase['Q'] for ii in range(svr)])]
    xlab=['PWN salinity','PWS salinity','PWS salinity']
    titvec=['PWS = 34.4','PWN = 33.7','']
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
        ax1.set_ylim(-30,30)
        ax2.set_ylim(-30,30)
        axx[ii].set_title(titvec[ii])
        if ii!=2:
            ax1.set_yticklabels('')
            ax2.set_yticklabels('')
        axx[ii].set_xlim(xvar[ii][0],xvar[ii][-1])
    axx[0].set_ylim(-1.1,1.1)
    axx[0].set_yticks(arange(-1,1.1,0.5))
    ax2.spines["right"].set_position(("axes", 1.3))
    axx[0].set_ylabel('Transport anomaly [Sv]')
    ax1.set_ylabel('Fresh water transport anomaly [mSv]',color='c')
    ax2.set_ylabel('Heat flux anomaly [TW]',color='limegreen')
    ax1.tick_params(axis='y', colors='c')
    ax2.tick_params(axis='y', colors='limegreen')
    leg=axx[0].legend(loc=(0.7,-0.4),ncol=5,fontsize=13)
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    axi2=axx[2].twiny()
    axi2.set_xticks(arange(32.75,34,0.25))
    axi2.set_xlim(xvar[0][0],xvar[0][-1])
    axi2.set_xlabel('PWN salinity')


lineplot_PW_salinity()


# Maintain difference between them, but make both fresher:
# INFLOWS
# More FW
# More PWN
# Less AWS
# OUTFLOWS
# Less PWS
# More DWS
# More AWN
# Extract less heat (since bringing less AWS in...?)
# How to understand this:
# More PWS to make than PWN coming in, so need extra freshwater coming in, explains upping FW and PWN, making less PWS, bringing in less salty AWS.
# That water has to go somewhere, so produce more AWN and DWS, which are less fresh.

# Fresher PWN than PWS:
# INFLOWS
# Less FW needed
# Less PWN brought in
# More AWS brought in
# OUTFLOWS
# More PWS produced (to a point)
# Less DWS produced
# Less AWN produced
# More heat extracted (because of more AWS I suspect)

#######################################################################################
##############  What happens if we add more FW? (Like 100mSv)  ###########################
#######################################################################################
Ubase['FW']
Ubase['SI']

xbase
get_U_from_x(xbase)
xbase

fwvar={}
for U_FW in arange(0,0.11,0.01):
        AM=array([[1,1,1,1,1,1,1,0],\
        [S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
        [T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

        xinit=array([-3.31214887e+00,  1.73787555e+01, -1.29083726e+01,  1.81280407e+00,
       -3.06611941e+00,  2.50406419e-02+U_FW,  7.00406419e-02, -6.36153015e+01])

        dv=-AM.dot(xinit)

        Evec=hstack((5*[1],0.02,0.02,Qvar))
        E=diag(Evec)
        Winv=diag([1,1/Snorm,1/Tnorm])
        Umat,D,VmatT=linalg.svd(Winv.dot(AM.dot(E)))

        Lambda_inv = zeros((AM.shape[0], AM.shape[1])).T
        Lambda_inv[:AM.shape[0], :AM.shape[0]] = diag(1/D)
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv.dot(dv))))
        xsol_Ad=E.dot(xsol_prime)
        fwvar[U_FW]=xinit+xsol_Ad

U_fwvar=get_U_from_x(fwvar[0.05])
a_fw,b_fw=get_a_b_fracs(U_fwvar,S)
U['FW']+U['SI']
Ubase['FW']+Ubase['SI']
U_fwvar['FW']+U_fwvar['SI']

pwmats

U_fwvar['Q']*cp*rhow/1e6
#######################################################################################
##############  Now look at consequences  for FW dist  ###########################
#######################################################################################
a_pwmat=zeros((len(epsilon),shape(pwmats['Q'])[1],shape(pwmats['Q'])[0]))
b_pwmat=a_pwmat.copy()
for ii,ee in enumerate(1-epsilon):
    a_pwmat[ii,:,:]=(ee*pwmats['PWN'].T*((S['PWN']+PWN_Smat)/S['AWS']-1)+pwmats['PWS'].T*((S['PWS']+PWS_Smat)/S['AWS']-1))/(pwmats['FW'].T+pwmats['SI'].T)
    b_pwmat[ii,:,:]=((1-ee)*pwmats['PWN'].T*((S['PWN']+PWN_Smat)/S['AWS']-1)+pwmats['DWS'].T*(S['DWS']/S['AWS']-1))/(pwmats['FW'].T+pwmats['SI'].T)
c_pwmat=1-a_pwmat-b_pwmat

epsilon=arange(0,1.1,0.1)

fwcol='#43a2ca'

def plot_adep_pw():
    f,axx=subplots(1,2,figsize=(11,3.2),sharex=True)
    f.subplots_adjust(wspace=0.3)
    for ii,var in enumerate([a_pwmat,b_pwmat]):
        axx[ii].plot(1-epsilon,var[:,-2,-3],linewidth=4,color='k',label='Base case',zorder=5)
        axx[ii].plot(1-epsilon,var[:,5,5],color='purple',zorder=4,label='Freshen both Polar Waters by 0.5',linewidth=3)
        axx[ii].set_ylim(-0.3,1.3)
        axx[ii].set_yticks(arange(0,1.1,0.5))
    axx[0].plot(1-epsilon,a_fw,linewidth=3,color=fwcol,label='Add 50 mSv of Fresh Water')
    axx[1].plot(1-epsilon,b_fw,linewidth=3,color=fwcol)
    axx[0].legend(loc=(0.25,-0.5),ncol=3)
    axx[0].set_title('Estuarine limb',fontsize=14)
    axx[1].set_title('Overturning limb',fontsize=14)
    axx[0].set_ylabel('$\mathbf{\delta}}$,  FW fraction in PWS')
    axx[1].set_ylabel('$\mathbf{\gamma}$,  FW fraction in DWS')
    for axi in axx:
        axi.set_xlabel('$\mathbf{1-\epsilon}$\nFraction of PWN in PWS')
        axi.axhline(0,color='k',linestyle='--')
        # axi.axhline(0.5,color='k',linestyle='--')
        axi.set_xlim(0,1)
    savefig(figdir_paper+'/FWfrac_obs_pwdep.png',bbox_inches='tight')
    savefig(figdir_paper+'/FWfrac_obs_pwdep.pdf',bbox_inches='tight')

plot_adep_pw()
