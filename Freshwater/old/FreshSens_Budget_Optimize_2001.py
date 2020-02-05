from firstfuncs_1618 import *
from scipy.optimize import minimize
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/optimize/'

########################################################################################################
########################################################################################################
#### Set up the optimization framework, which allows for varying almost all elements within a prescribed range
########################################################################################################
########################################################################################################

# x = [ U_AWS, U_PWS, U_DWS, U_AWN, U_PWN, U_FW, U_SI, S_AWS, S_PWS, S_DWS, S_AWN, S_PWN]
xnames=[ 'U_AWS', 'U_PWS', 'U_DWS', 'U_AWN', 'U_PWN', 'U_FW', 'U_SI', 'S_AWS', 'S_PWS', 'S_DWS', 'S_AWN', 'S_PWN'] #so i can check
# xg = same order, first guess for each of the above
# delx = same order, but a vector which describes their allowed +/- range

WM=xr.open_dataset(datadir+'OSNAP2016recovery/pickles/gridded/OSNAP2014-16_WM_1912.nc')
WM['TRANS'].mean(dim='TIME')
WM['PSAL'].mean(dim='TIME')
# first guesses and ranges for each of the above:
xg=array([16.9,-3.2,-12.1,-9.6,8.4,0.05,0.06,35.2,34.4,35,35,34.6])

# delx=array([5.4,2.25,6.0,2,2,0.04,0.03,0.06,0.8,0.1,0.2,0.8])
delx=array([2,2,2,2,2,0.01,0.01,0.06,0.8,0.1,0.2,0.8])

gammavec=[4,20,100,250,500,1100]

def get_U_S_from_x(x):
# x = [ U_AWS, U_PWS, U_DWS, U_AWN, U_PWN, U_FW, U_SI, S_AWS, S_PWS, S_DWS, S_AWN, S_PWN]
    U={}
    S={}
    U['AWS']=x[0]
    U['PWS']=x[1]
    U['DWS']=x[2]
    U['AWN']=x[3]
    U['PWN']=x[4]
    U['FW']=x[5]
    U['SI']=x[6]
    S['AWS']=x[7]
    S['PWS']=x[8]
    S['DWS']=x[9]
    S['AWN']=x[10]
    S['PWN']=x[11]
    S['FW']=0
    S['SI']=6
    return U,S

Ug,Sg=get_U_S_from_x(xg)

delx=array([2,2,2,2,2,0.01,0.01,0.06,0.8,0.1,0.2,0.8])

Ue={}
Se={}
Ue['AWS']=delx[0]
Ue['PWS']=delx[1]
Ue['DWS']=delx[2]
Ue['AWN']=delx[3]
Ue['PWN']=delx[4]
Ue['FW']=delx[5]
Ue['SI']=delx[6]
Se['AWS']=delx[7]
Se['PWS']=delx[8]
Se['DWS']=delx[9]
Se['AWN']=delx[10]
Se['PWN']=delx[11]
Se['FW']=0
Se['SI']=0


def BudgetOut(x,throughflow='false'):
    ms=x[0]+x[1]+x[2]+x[5]/2
    mn=x[3]+x[4]+x[6]+x[5]/2
    sall=x[0]*x[7]+x[1]*x[8]+x[2]*x[9]+x[3]*x[10]+x[4]*x[11]+6*x[6]
    return ms,mn,sall

def plot_budget_out(casename,phrase):
    f,axx=subplots(4,1,figsize=(8,12),sharex=True)
    ii=0
    for kk in res:
        if phrase in kk:
            ms,mn,sall=BudgetOut(res[kk].x)
            axx[0].plot(ii,ms,'o',label=kk)
            axx[1].plot(ii,mn,'o',label=kk)
            axx[2].plot(ii,mn+ms,'o',label=kk)
            axx[3].plot(ii,sall,'o',label=kk)
            ii+=1
    for axi in axx:
        axi.axhline(0,color='grey')
    axx[1].legend(loc=(1.05,0))
    axx[2].set_xticklabels('')
    axx[0].set_ylabel('Southern\n volume flux [Sv]')
    axx[1].set_ylabel('Northen\n volume flux [Sv]')
    axx[2].set_ylabel('Total volume divergence [Sv]')
    axx[3].set_ylabel('Salinity balance [sal x Sv]')
    axx[0].set_title('Balance residuals for '+casename)
    savefig(figdir+'Gamma_dependence_'+casename+'.png',bbox_inches='tight')


def plot_budget_out_simple(casename,phrase):
    f,axx=subplots(2,1,figsize=(8,6),sharex=True)
    ii=0
    for kk in res:
        if phrase in kk:
            ms,mn,sall=BudgetOut(res[kk].x)
            axx[0].plot(ii,mn+ms,'o',label=kk)
            axx[1].plot(ii,sall,'o',label=kk)
            ii+=1
    for axi in axx:
        axi.axhline(0,color='grey')
    axx[1].legend(loc=(1.05,0))
    axx[0].set_xticklabels('')
    axx[0].set_ylabel('Total volume divergence [Sv]')
    axx[1].set_ylabel('Salinity balance [sal x Sv]')
    axx[0].set_title('Balance residuals for '+casename)
    savefig(figdir+'Gamma_dependence_simple'+casename+'.png',bbox_inches='tight')


colvec='grey'

def barplot(casename,phrase):
    sref=34.8
    f,axx=subplots(4,2,figsize=(16,11))
    ax1=axx[0,0]
    ax2=axx[1,0]
    ax3=axx[2,0]
    ax4=axx[3,0]
    U,S=get_U_S_from_x(xg)
    Ug=U.copy()
    Sg=S.copy()
    varlen=len(U.keys())
    alf=0.4
    ax1.bar(range(varlen),[U[kk] for kk in U],color=colvec,yerr=[Ue[kk] for kk in U],capsize=10,alpha=alf)
    ax1.set_ylim(-25,25)
    ax1.axvline(4.5,color='k')
    ax1_1=ax1.twinx()
    ax1_1.bar([5,6],[U[kk] for kk in ['FW','SI']],color=colvec,yerr=[Ue[kk] for kk in ['FW','SI']],capsize=10,alpha=alf)
    ax1_1.set_ylim(-0.2,0.2)
    ax2.bar(range(varlen),[S[kk] for kk in S],color=colvec,yerr=[Se[kk] for kk in U],capsize=10,alpha=alf)
    ax3.bar(range(varlen),[U[kk]*S[kk] for kk in S],color=colvec,yerr=[Ue[kk]*S[kk]+Se[kk]*U[kk] for kk in U],capsize=10,alpha=alf)
    ax3.set_ylim(-700,700)
    ax3_1=ax3.twinx()
    ax3_1.bar([5,6],[U[kk]*S[kk] for kk in ['FW','SI']],color=colvec,yerr=[Ue[kk]*S[kk]+Se[kk]*U[kk] for kk in ['FW','SI']],capsize=10,alpha=alf)
    ax3_1.axvline(4.5,color='k')
    ax3_1.set_ylim(-0.5,0.5)
    ax4.bar(range(varlen),[U[kk]*(S[kk]-34.8)/34.8 for kk in S],color=colvec,yerr=[Ue[kk]*(S[kk]-sref)/sref+U[kk]*Se[kk]/sref for kk in S],capsize=10,alpha=alf)
    ax2.set_ylim(31,36)
    ax2.axhline(sref,color='grey')
    # ax1.text(4.6,5,str(U['FW']),fontsize=24)
    # ax1.text(5.6,5,str(U['SI']),fontsize=24)
    ax2.text(4.9,32.6,'0',fontsize=24)
    ax2.text(5.9,32.6,'6',fontsize=24)
    # ax3.text(4.6,100,str(U['SI']*S['SI']),fontsize=24)
    # ax3.text(5.9,100,'0',fontsize=24)
    for case in res:
        if phrase in case:
            U,S=get_U_S_from_x(res[case].x)
            ax1.plot(range(varlen)[:5],[U[kk] for kk in U][:5],'o-',label=case)
            ax1_1.plot(range(varlen)[5:],[U[kk] for kk in U][5:],'o-',label=case)
            ax2.plot(range(varlen)[:5],[S[kk] for kk in S][:5],'o-',label=case)
            ax3.plot(range(varlen)[:5],[U[kk]*S[kk] for kk in S][:5],'o-',label=case)
            ax3_1.plot(range(varlen)[5:],[U[kk]*S[kk] for kk in S][5:],'o-',label=case)
            ax4.plot(range(varlen),[U[kk]*(S[kk]-sref)/sref for kk in S],'o-',label=case)
            axx[0,1].plot(range(varlen),[U[kk]-Ug[kk] for kk in U],'o-',label=case)
            axx[1,1].plot(range(varlen)[:5],[S[kk]-Sg[kk] for kk in S][:5],'o-',label=case)
            axx[2,1].plot(range(varlen),[(U[kk]*S[kk])-(Ug[kk]*Sg[kk]) for kk in S],'o-',label=case)
            axx[3,1].plot(range(varlen),[(U[kk]*(S[kk]-sref)-Ug[kk]*(Sg[kk]-sref))/sref for kk in S],'o-',label=case)
    axx[0,1].set_title('Anomalies from initial conditions')
    # ax2.legend(loc='upper right',)
    axx[0,1].legend(loc=(1.05,0))
    for axxi in [ax1,ax2,ax3,ax4,axx[0,1],axx[1,1],axx[2,1],axx[3,1]]:
        axxi.set_xticks(range(varlen))
        axxi.set_xticklabels(U.keys())
        axxi.axhline(0,color='k')
    ax1.set_title('Freshwater Budget Example')
    ax1.set_ylabel('Transport, U [Sv]')
    ax2.set_ylabel('Salinity, S')
    ax3.set_ylabel('Salinity budget term, U*S')
    ax4.set_ylabel('Using '+str(sref)+' reference salinity,\n U(S-'+str(sref)+')/'+str(sref))
    savefig(figdir+'FWBudget_Optimize_'+casename+'.png',bbox_inches='tight')

#############################################################
# Get solutions for a couple different gammas:
#############################################################
gammavec

res={}
for gamma in gammavec:
    def CostFunction(x,delx=delx): # note: gamma is a measure of the importance of satisfying equations vs. matching first guess
        J=(x[0]+x[1]+x[2]+x[5]/2)**2
        +(x[3]+x[4]+x[6]+x[5]/2)**2 \
        +(x[0]*x[7]+x[1]*x[8]+x[2]*x[9]+x[3]*x[10]+x[4]*x[11]+6*x[6])**2 \
        +sum((x-xg)**2/delx**2/gamma**2) #force each solution to be within delx of the first guess
        return J
        ## Minimize the cost function:
    res['volbal, gamma='+str(gamma)]=minimize(CostFunction,xg,method='powell',options={'xtol':1e-8,'disp': True})

plot_budget_out('volbal','volbal')


barplot('volbal','volbal')
#############################################################
# Get solutions for a couple different reference salinities
#############################################################
for gamma in gammavec:
    for srefy in [34.8]:
        def CostFunction_w_sref(x,delx=delx):
            J=(x[0]+x[1]+x[2]+x[5]/2)**2 \
            +(x[3]+x[4]+x[6]+x[5]/2)**2 \
            +(x[0]*x[7]+x[1]*x[8]+x[2]*x[9]+x[3]*x[10]+x[4]*x[11]+x[6]*6)**2/srefy**2 \
            +sum((x-xg)**2/delx**2/gamma**2) #force each solution to be within delx of the first guess
            return J
        res['sref='+str(srefy)+', gamma='+str(gamma)]=minimize(CostFunction_w_sref,xg,method='powell',options={'xtol':1e-8,'disp': True})
plot_budget_out('ref sal','sref=34.8')
barplot('ref sal','sref=34.8')
#throughflow of 1Sv for simplicity
bs=1

for gamma in gammavec:
    def CostFunction_w_Tflow(x,delx=delx):# note: gamma is a measure of the importance of satisfying equations vs. matching first guess
        sref=34.8
        J=(x[0]+x[1]+x[2]+x[5]/2-bs)**2 \
        +(x[3]+x[4]+x[6]+x[5]/2+bs)**2 \
        +(x[0]*x[7]+(x[1])*x[8]+(x[2])*x[9]+x[3]*x[10]+(x[4])*x[11]+x[6]*6)**2 \
        +sum((x-xg)**2/delx**2/gamma**2) #force each solution to be within delx of the first guess
        return J
    res['throughflow, gamma='+str(gamma)]=minimize(CostFunction_w_Tflow,xg,method='powell',options={'xtol':1e-8,'disp': True})

plot_budget_out('throughflow','throughflow')
barplot('throughflow','throughflow')

for gamma in gammavec:
    def CostFunction_DivOnly(x,delx=delx):# note: gamma is a measure of the importance of satisfying equations vs. matching first guess
        sref=34.8
        J=(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6])**2 \
        +(x[0]*x[7]+(x[1])*x[8]+(x[2])*x[9]+x[3]*x[10]+(x[4])*x[11]+x[6]*6)**2 \
        +sum((x-xg)**2/delx**2/gamma**2) #force each solution to be within delx of the first guess
        return J
    res['div_only, gamma='+str(gamma)]=minimize(CostFunction_DivOnly,xg,method='powell',options={'xtol':1e-8,'disp': True})
plot_budget_out('div_only','div_only')
plot_budget_out_simple('div_only','div_only')
barplot('div_only','div_only')

barplot('gamma04','gamma=4')
barplot('gamma20','gamma=20')

plot_budget_out('gamma20','gamma=20')


#############################################################
# TEST SENSITIVITY TO INIITIAL CONDITIONS
#############################################################
gamma=100
#going to do this for div_only case first, then try others.
for ii in range(10):
    res['div_only, gamma='+str(gamma)+', initest '+str(ii)]=minimize(CostFunction_DivOnly,xg+delx*(2*random(len(delx))-1),method='powell',options={'xtol':1e-8,'disp': True})
    res['volbal, gamma='+str(gamma)+', initest '+str(ii)]=minimize(CostFunction,xg+delx*(2*random(len(delx))-1),method='powell',options={'xtol':1e-8,'disp': True})
    res['throughflow, gamma='+str(gamma)+', initest '+str(ii)]=minimize(CostFunction_w_Tflow,xg+delx*(2*random(len(delx))-1),method='powell',options={'xtol':1e-8,'disp': True})


plot_budget_out('initest_divonly_g100','div_only, gamma=100, initest')
barplot('initest_divonly_g100','div_only, gamma=100, initest')


plot_budget_out('initest_throughflow_g100','throughflow, gamma=100, initest')
barplot('initest_throughflow_g100','throughflow, gamma=100, initest')

plot_budget_out('initest_volbal_g100','volbal, gamma=100, initest')
barplot('initest_volbal_g100','volbal, gamma=100, initest')

Sg['AWS']

def plot_ref_to_input():
    varlen=len(Ug.keys())
    colvec=['red','royalblue','grey','purple','orange','black','cyan']
    bar(range(varlen),[Ug[kk]*(Sg[kk]-Sg['AWS'])/Sg['AWS']*1000 for kk in Sg],color=colvec)#,yerr=[Ue[kk]*(S[kk]-['AWS'])/['AWS']+U[kk]*Se[kk]/sref for kk in S],capsize=10,alpha=alf)
    gca().set_xticks(range(varlen))
    gca().set_xticklabels(Ug.keys())
    ylabel('[mSv]')
    title('Freshwater transport relative to AWS salinity, 35.2',fontsize=14)
    savefig(figdir+'Budget_interpretation_example.png')

plot_ref_to_input()
