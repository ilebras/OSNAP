from firstfuncs_1618 import *
from scipy.optimize import minimize
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/Linear/'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################

WM=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')

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
    S['SI']=6
    S['FW']=0
    T['SI']=-90
    T['FW']=0


    return U,S,T

U,S,T=get_U_S_T_from_WM()

Amat=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

xbase=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]
Q=-65

Usum=0
SUsum=0
TUsum=Q
for wm in U:
    Usum=Usum+U[wm]
    SUsum=SUsum+U[wm]*S[wm]
    TUsum=TUsum+U[wm]*T[wm]

dvec=[-Usum,-SUsum,-TUsum]

#############################################################
# Getting SVD
#############################################################

#weights:
Snorm=35
Tnorm=5
W=diag([1,1/Snorm,1/Tnorm])

Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values/10,0.025,0.025,10))
E=diag(Evec)

Umat,D,VmatT=linalg.svd(W.dot(Amat.dot(E)))
Sigma = zeros((Amat.shape[0], Amat.shape[1]))
# populate Sigma with n x n diagonal matrix
Sigma[:Amat.shape[0], :Amat.shape[0]] = diag(D)
Sigma_inv=Sigma.T
Sigma_inv[:Amat.shape[0], :Amat.shape[0]] = diag(1/D)

xsol_prime=VmatT.T.dot(Sigma_inv.dot(Umat.T.dot(W.dot(dvec))))

xsol=E.dot(xsol_prime)

###quick checks:
W.dot(Amat.dot(E.dot(inv(E).dot(xsol))))
W.dot(dvec)

Amat.dot(xsol)
dvec


#############################################################
# Plotting infrastructure
#############################################################
xbase=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]

def BudgetOut(x):
    mall=x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]
    sall=(x[0]*S['PWS']+x[1]*S['AWS']+x[2]*S['DWS']+x[3]*S['PWN']+x[4]*S['AWN']+x[6]*S['SI'])
    tall=(x[0]*T['PWS']+x[1]*T['AWS']+x[2]*T['DWS']+x[3]*T['PWN']+x[4]*T['AWN']+x[5]*T['FW']+x[6]*T['SI']+x[7])
    return mall,sall,tall

mb,sb,tb=BudgetOut(xbase)
mg,sg,tg=BudgetOut(xbase+xsol)

mb,mg
sb,sg
tb,tg

def plot_budget_out(casename):
    f,axx=subplots(3,1,figsize=(8,9),sharex=True)
    ii=0
    for kk in res:
        if casename in kk:
            mall,sall,tall=BudgetOut(res[kk].x)
            axx[0].plot(ii,mall,'o',label=kk)
            axx[1].plot(ii,sall,'o',label=kk)
            axx[2].plot(ii,tall,'o',label=kk)
            ii+=1
    for axi in axx:
        axi.axhline(0,color='grey')
    axx[1].legend(loc=(1.05,0))
    axx[2].set_xticklabels('')
    axx[0].set_ylabel('Volume balance [Sv]')
    axx[1].set_ylabel('Salinity balance [S x Sv]')
    axx[2].set_ylabel('Temperature balance [$^\circ$C x Sv]')
    axx[0].set_title('Balance residuals for '+casename)
    savefig(figdir+'Gamma_dependence_'+casename+'.png',bbox_inches='tight')

colvec='grey'

def get_U_from_x(x):
    U={}
    U['PWS']=x[0]
    U['AWS']=x[1]
    U['DWS']=x[2]
    U['PWN']=x[3]
    U['AWN']=x[4]
    U['FW']=x[5]
    U['SI']=x[6]
    return U

Ue=get_U_from_x(Evec[:-1])

def barplot(casename):
    Ue=get_U_from_x(Evec[:-1])
    Ug=U.copy()
    Sg=S.copy()
    Tg=T.copy()
    f,axx=subplots(3,2,figsize=(16,10),sharex='col')
    varlen=len(Ug)
    alf=0.4
    #U
    axx[0,0].bar(range(varlen)[:-2],[Ug[kk] for kk in Ug][:-2],color=colvec,yerr=[Ue[kk] for kk in Ug][:-2],capsize=10,alpha=alf)
    axx[0,0].set_ylim(-25,25)
    axx[0,0].axvline(4.5,color='k')
    ax1=axx[0,0].twinx()
    ax1.bar([5,6],[Ug[kk] for kk in ['FW','SI']],color=colvec,yerr=[Ue[kk] for kk in ['FW','SI']],capsize=10,alpha=alf)
    ax1.set_ylim(-0.2,0.2)
    #UxS
    axx[1,0].bar(range(varlen),[Ug[kk]*Sg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Sg[kk] for kk in Ug],capsize=10,alpha=alf)
    axx[1,0].set_ylim(-750,750)
    #UxT
    axx[2,0].bar(range(varlen),[Ug[kk]*Tg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Tg[kk] for kk in Ug],capsize=10,alpha=alf)
    axx[2,0].bar(varlen,xbase[-1],color=colvec,yerr=Evec[-1],capsize=10,alpha=alf)
    for axi in [axx[0,0],axx[1,0]]:
        axi.axvline(4.5,color='k')
    # axx[2,0].axvline(5.5,color='k')
    ii=0
    for case in res:
        if casename in case:
            Ux=get_U_from_x(res[case])
            axx[0,0].plot(range(varlen)[:5],[Ux[kk] for kk in U][:5],'o-',label=case,color='C'+str(ii))
            ax1.plot(range(varlen)[5:],[Ux[kk] for kk in U][5:],'o-',label=case,color='C'+str(ii))
            axx[1,0].plot(range(varlen),[Ux[kk]*S[kk] for kk in S],'o-',label=case,color='C'+str(ii))
            axx[2,0].plot(range(varlen),[Ux[kk]*T[kk] for kk in S],'o-',label=case,color='C'+str(ii))
            axx[2,0].plot(varlen,res[case][-1],'o-',label=case,color='C'+str(ii))
            axx[0,1].plot(range(varlen),[Ux[kk]-Ug[kk] for kk in U],'o-',label=case,color='C'+str(ii))
            axx[1,1].plot(range(varlen),[(Ux[kk]*S[kk])-(Ug[kk]*Sg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
            axx[2,1].plot(range(varlen),[(Ux[kk]*T[kk])-(Ug[kk]*Tg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
            axx[2,1].plot(varlen,res[case][-1]-xbase[-1],'o-',label=case,color='C'+str(ii))
            ii+=1
    axx[0,1].set_title('Anomalies from initial conditions')
    # ax2.legend(loc='upper right',)
    axx[0,1].legend(loc=(1.05,0))
    for axesi in axx:
        for axxi in axesi:
            axxi.set_xticks(range(varlen+1))
            axxi.set_xticklabels(Tg.keys())
            axxi.axhline(0,color='k')
    axx[0,0].set_ylabel('Transport, U [Sv]')
    axx[1,0].set_ylabel('Salinity budget term\n UxS [Sv]')
    axx[2,0].set_ylabel('Temperature budget term\n UxT [$^\circ$C Sv]')

    savefig(figdir+'FWBudget_Linear_'+casename+'.png',bbox_inches='tight')


res={}
res['base']=xbase+xsol

barplot('base')
