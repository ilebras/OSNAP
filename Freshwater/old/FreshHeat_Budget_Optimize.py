from firstfuncs_1618 import *
from scipy.optimize import minimize
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/optimize/'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################

WM=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')
xnames=hstack((['U_'+str(ww) for ww in WM.WM.values],'U_FW','U_SI',['S_'+str(ww) for ww in WM.WM.values],'S_SI',['T_'+str(ww) for ww in WM.WM.values],'T_FW','T_SI','q_T'))
# xg = same order, first guess/initial condition for each of the above
# delx = same order, but a vector which describes their allowed +/- range


# first guesses and ranges for each of the above:
xg=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,WM['PSAL'].mean(dim='TIME').values,6,WM['PTMP'].mean(dim='TIME').values,0,-90,-65))


len(xg)
len(delx)
len(xnames)

def get_U_S_T_from_x(x):
# x = [ U_AWS, U_PWS, U_DWS, U_AWN, U_PWN, U_FW, U_SI, S_AWS, S_PWS, S_DWS, S_AWN, S_PWN]
    U={}
    S={}
    T={}
    U['PWS']=x[0]
    U['AWS']=x[1]
    U['DWS']=x[2]
    U['PWN']=x[3]
    U['AWN']=x[4]
    U['FW']=x[5]
    U['SI']=x[6]
    S['PWS']=x[7]
    S['AWS']=x[8]
    S['DWS']=x[9]
    S['PWN']=x[10]
    S['AWN']=x[11]
    S['FW']=0
    S['SI']=x[12]
    T['PWS']=x[13]
    T['AWS']=x[14]
    T['DWS']=x[15]
    T['PWN']=x[16]
    T['AWN']=x[17]
    T['FW']=x[18]
    T['SI']=x[19]
    T['$q_T$']=x[20]
    return U,S,T

Ug,Sg,Tg=get_U_S_T_from_x(xg)
Ue,Se,Te=get_U_S_T_from_x(delx)


def BudgetOut(x):
    mall=x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]
    sall=(x[0]*x[7]+x[1]*x[8]+x[2]*x[9]+x[3]*x[10]+x[4]*x[11]+x[6]*x[12])
    tall=(x[0]*x[13]+x[1]*x[14]+x[2]*x[15]+x[3]*x[16]+x[4]*x[17]+x[5]*x[18]+x[6]*x[19]+x[20])
    return mall,sall,tall

mg,sg,tg=BudgetOut(xg)

#############################################################
# Plotting infrastructure
#############################################################

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

def barplot(casename):
    delx=delxvec[casename]
    Ue,Se,Te=get_U_S_T_from_x(delx)
    f,axx=subplots(5,2,figsize=(16,16),sharex='col')
    varlen=len(Ug)
    alf=0.4
    #U
    axx[0,0].bar(range(varlen)[:-2],[Ug[kk] for kk in Ug][:-2],color=colvec,yerr=[Ue[kk] for kk in Ug][:-2],capsize=10,alpha=alf)
    axx[0,0].set_ylim(-25,25)
    axx[0,0].axvline(4.5,color='k')
    ax1=axx[0,0].twinx()
    ax1.bar([5,6],[Ug[kk] for kk in ['FW','SI']],color=colvec,yerr=[Ue[kk] for kk in ['FW','SI']],capsize=10,alpha=alf)
    ax1.set_ylim(-0.2,0.2)
    #S
    axx[1,0].bar(range(varlen-2),[Sg[kk] for kk in Ug][:-2],color=colvec,yerr=[Se[kk] for kk in Ug][:-2],capsize=10,alpha=alf)
    # axx[1,0].set_ylim(33.5,35.5)
    ax2=axx[1,0].twinx()
    ax2.bar([6],Sg['SI'],color=colvec,yerr=Se['SI'],capsize=10,alpha=alf)
    #T
    axx[2,0].bar(range(varlen-1),[Tg[kk] for kk in Ug][:-1],color=colvec,yerr=[Te[kk] for kk in Ug][:-1],capsize=10,alpha=alf)
    axx[2,0].set_ylim(-10,10)
    ax3=axx[2,0].twinx()
    ax3.bar(6,Tg['SI'],color=colvec,yerr=Te['SI'],capsize=10,alpha=alf)
    ax3.set_ylim(-100,100)
    #UxS
    axx[3,0].bar(range(varlen),[Ug[kk]*Sg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Sg[kk]+Se[kk]*Ug[kk] for kk in Ug],capsize=10,alpha=alf)
    axx[3,0].set_ylim(-750,750)
    #UxT
    axx[4,0].bar(range(varlen),[Ug[kk]*Tg[kk] for kk in Ug],color=colvec,yerr=[Ue[kk]*Tg[kk]+Te[kk]*Ug[kk] for kk in Ug],capsize=10,alpha=alf)
    axx[4,0].bar(varlen,xg[-1],color=colvec,yerr=delx[-1],capsize=10,alpha=alf)
    for axi in [axx[0,0],axx[1,0]]:
        axi.axvline(4.5,color='k')
    axx[2,0].axvline(5.5,color='k')
    ii=0
    for case in res:
        if casename in case:
            U,S,T=get_U_S_T_from_x(res[case].x)
            axx[0,0].plot(range(varlen)[:5],[U[kk] for kk in U][:5],'o-',label=case,color='C'+str(ii))
            ax1.plot(range(varlen)[5:],[U[kk] for kk in U][5:],'o-',label=case,color='C'+str(ii))
            axx[1,0].plot(range(varlen)[:5],[S[kk] for kk in S][:5],'o-',label=case,color='C'+str(ii))
            ax2.plot(range(varlen)[5:],[S[kk] for kk in S][5:],'o-',label=case,color='C'+str(ii))
            axx[2,0].plot(range(varlen)[:6],[T[kk] for kk in S][:6],'o-',label=case,color='C'+str(ii))
            ax3.plot(range(varlen)[6:],[T[kk] for kk in S][6:],'o-',label=case,color='C'+str(ii))
            axx[3,0].plot(range(varlen),[U[kk]*S[kk] for kk in S],'o-',label=case,color='C'+str(ii))
            axx[4,0].plot(range(varlen),[U[kk]*T[kk] for kk in S],'o-',label=case,color='C'+str(ii))
            axx[4,0].plot(varlen,res[case].x[-1],'o-',label=case,color='C'+str(ii))
            axx[0,1].plot(range(varlen),[U[kk]-Ug[kk] for kk in U],'o-',label=case,color='C'+str(ii))
            axx[1,1].plot(range(varlen),[S[kk]-Sg[kk] for kk in S],'o-',label=case,color='C'+str(ii))
            axx[2,1].plot(range(varlen),[T[kk]-Tg[kk] for kk in S],'o-',label=case,color='C'+str(ii))
            axx[3,1].plot(range(varlen),[(U[kk]*S[kk])-(Ug[kk]*Sg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
            axx[4,1].plot(range(varlen),[(U[kk]*T[kk])-(Ug[kk]*Tg[kk]) for kk in S],'o-',label=case,color='C'+str(ii))
            axx[4,1].plot(varlen,res[case].x[-1]-xg[-1],'o-',label=case,color='C'+str(ii))
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
    axx[1,0].set_ylabel('Salinity, S')
    axx[2,0].set_ylabel('Temperature, T [$^\circ$C]')
    axx[3,0].set_ylabel('Salinity budget term\n UxS [Sv]')
    axx[4,0].set_ylabel('Temperature budget term\n UxT [$^\circ$C Sv]')

    savefig(figdir+'FWBudget_Optimize_'+casename+'.png',bbox_inches='tight')


#############################################################
# Get solutions for a couple different gammas:
#############################################################

gammavec=[4,20,100,250,500,1100]

#Scales to normalize S and T equations by so that they are not favored over volume balance:
Snorm=35
Tnorm=10

#Vary the error allowances to see sensitivity to each type of parameter first:
delxvec={}
# delxvec['vary all']=hstack((WM['TRANS'].std(dim='TIME'),0.01,0.01,WM['PSAL'].std(dim='TIME'),0.1,WM['PTMP'].std(dim='TIME'),0.5,1,10))
fac=1e5
delxvec['varysal']=hstack((WM['TRANS'].std(dim='TIME')/fac,0.01/fac,0.01/fac,WM['PSAL'].std(dim='TIME'),0.1,WM['PTMP'].std(dim='TIME')/fac,0.5/fac,1/fac,10/fac))
delxvec['varytmp']=hstack((WM['TRANS'].std(dim='TIME')/fac,0.01/fac,0.01/fac,WM['PSAL'].std(dim='TIME')/fac,0.1/fac,WM['PTMP'].std(dim='TIME'),0.5,1,10))
delxvec['varytrans']=hstack((WM['TRANS'].std(dim='TIME'),0.01,0.01,WM['PSAL'].std(dim='TIME')/fac,0.1/fac,WM['PTMP'].std(dim='TIME')/fac,0.5/fac,1/fac,10/fac))
delxvec['varytrans_plusheat']=hstack((WM['TRANS'].std(dim='TIME'),0.01,0.01,WM['PSAL'].std(dim='TIME')/fac,0.1/fac,WM['PTMP'].std(dim='TIME')/fac,0.5/fac,1/fac,10))




res={}
for dd in delxvec:
    delx=delxvec[dd]
    for gamma in gammavec:
        def CostFunction(x):
            mall=x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]
            sall_norm=(x[0]*x[7]+x[1]*x[8]+x[2]*x[9]+x[3]*x[10]+x[4]*x[11]+x[6]*x[12])/Snorm
            tall_norm=(x[0]*x[13]+x[1]*x[14]+x[2]*x[15]+x[3]*x[16]+x[4]*x[17]+x[5]*x[18]+x[6]*x[19]+x[20])/Tnorm
            J=(mall)**2 + (sall_norm)**2 + (tall_norm)**2 + sum((x-xg)**2/delx**2/gamma**2)
            return J
        res[dd+', gamma = '+str(gamma)]=minimize(CostFunction,xg+delx*(2*random(len(delx))-1),method='powell',options={'xtol':1e-8,'disp': True})
        m,s,t=BudgetOut(res[dd+', gamma = '+str(gamma)].x)
        print(gamma,m,s,t)


plot_budget_out('varytrans_plusheat')

25*13000

barplot('varytrans_plusheat')
