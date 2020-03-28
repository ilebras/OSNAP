from firstfuncs_1618 import *
from scipy.optimize import minimize
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


def get_U_var_from_x(x):
    U={}
    U['AWS']=x[0]
    U['PWS']=x[1]
    U['DWS']=x[2]
    U['AWN']=x[3]
    U['PWN']=x[4]
    U['FW']=x[5]
    U['SI']=x[6]
    U['Q']=x[7]
    U['alpha']=x[8]
    U['beta']=x[9]
    U['delta']=x[10]
    U['lambda']=x[11]

    return U


xnames=hstack((['U_'+str(ww) for ww in WM.WM.values],'U_FW','U_SI','q_T','alpha\n(PWN)','beta\n(FW)','delta\n(SI)','lambda\n(Q)'))

# xg = same order, first guess/initial condition for each of the above
# delx = same order, but a vector which describes their allowed +/- range



#############################################################
# Plotting
#############################################################
coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'k','Q':'limegreen'}


def barplot(casename):
    f=figure(figsize=(16,3))
    numvars=[5,2,1,4]
    gs = f.add_gridspec(1, 4,width_ratios=numvars,wspace=0.5)
    ax1 = f.add_subplot(gs[0])
    ax2 = f.add_subplot(gs[1])
    ax3 = f.add_subplot(gs[2])
    ax4 = f.add_subplot(gs[3])
    xi=res[casename].x
    Ug=get_U_var_from_x(xi)
    xinit=xg[casename]
    Uinit=get_U_var_from_x(xinit)
    xe=delxvec[casename]
    alf=0.8
    #U
    varvec={}
    varvec[0]=['AWS','PWS','DWS','AWN','PWN']
    varvec[1]=['FW','SI']
    varvec[2]=['Q']
    varvec[3]=xnames[-4:]
    ax1.bar(range(5),[Ug[kk] for kk in varvec[0]],color=[coldic[kk] for kk in varvec[0]],alpha=alf)
    ax1.plot(range(5),[Uinit[kk] for kk in varvec[0]],'ko')
    [ax1.plot([ww,ww],[xinit[ww]-xe[ww],xinit[ww]+xe[ww]],'k-') for ww in range(5)]
    ax1.axhline(0,color='k')
    ax2.bar(range(2),[Ug[kk] for kk in varvec[1]],color=[coldic[kk] for kk in varvec[1]],alpha=alf)
    ax2.plot(range(2),[Uinit[kk] for kk in varvec[1]],'ko')
    [ax2.plot([ww-5,ww-5],[xinit[ww]-xe[ww],xinit[ww]+xe[ww]],'k-') for ww in range(5,7)]
    ax2.set_ylim(0,0.07)
    ax3.bar(range(1),[Ug['Q']],color=[coldic['Q']],alpha=alf)
    ax3.plot(range(1),[Uinit['Q']],'ko')
    [ax3.plot([ww-7,ww-7],[xinit[ww]-xe[ww],xinit[ww]+xe[ww]],'k-') for ww in range(7,8)]
    ax4.bar(range(4),[Ug[kk] for kk in ['alpha','beta','delta','lambda']],color='grey',alpha=alf)
    ax4.plot(range(4),[Uinit[kk] for kk in ['alpha','beta','delta','lambda']],'ko')
    [ax4.plot([ww-8,ww-8],[xinit[ww]-xe[ww],xinit[ww]+xe[ww]],'k-') for ww in range(8,12)]
    ax4.axhline(0,color='k')
    for ii,axi in enumerate([ax1,ax2,ax3,ax4]):
        axi.set_xticks(range(numvars[ii]))
        axi.set_xticklabels(varvec[ii])
    f.suptitle(casename,fontsize=16)

#############################################################
# New equations including partitioning into 2 boxes
#############################################################

# First try setting up nonlinear system
# Unknowns are all Us,alpha,beta,delta,lambda, which are fractions of PWN,SI,FW, and Q going to PWS box
# Also bring back gamma, which weighs the importance of meeting initial conditions vs. equations -- could make this arbitrarily complicated.
# Could even take this part out-- and rely on first guess to settle things out? Test this.

def WriteOutFunction(x):
    print(x[1]+'+'+x[4]+'*'+x[8]+'+'+x[5]+'*'+x[9]+'+'+x[6]+'*'+x[10])
    print(x[0]+'+'+x[2]+'+'+x[3]+'+'+x[4]+'*(1-'+x[8]+')+'+x[5]+'*(1-'+x[9]+')+'+x[6]+'*(1-'+x[10]+')')

WriteOutFunction(xnames)
Snorm=35
Tnorm=10


# first guesses and ranges for each of the above:
xg={}
xg['base']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.4,0.1,0.8,0.1))
xg['fractest']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.5,0.5,0.5))
xg['FW_SI_rat']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.3,0.7,0.5))
xg['FW_SI_rat_Uless']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.3,0.7,0.5))
xg['U_static']=xg['FW_SI_rat'].copy()
xg['All_U_static']=xg['FW_SI_rat'].copy()
xg['All_U_static_FW_SI']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.2,0.8,0.5))
xg['All_U_static_FW_SI_plus']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.1,0.9,0.5))
xg['All_U_static_FW_SI_wide']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.1,0.9,0.5))
xg['U_static_FW_SI_wide']=hstack((WM['TRANS'].mean(dim='TIME').values,0.05,0.05,Q,0.5,0.1,0.9,0.5))

delxvec={}
fracvar=0.1
delxvec['base']=hstack((WM['TRANS'].std(dim='TIME'),0.01,0.01,10,4*[fracvar]))
fracvar=0.5
delxvec['fractest']=hstack((WM['TRANS'].std(dim='TIME'),0.01,0.01,10,4*[fracvar]))
fracvar=0.3
delxvec['FW_SI_rat']=hstack((WM['TRANS'].std(dim='TIME'),0.01,0.01,10,4*[fracvar]))
delxvec['FW_SI_rat_Uless']=hstack((5*[0.5],0.01,0.01,10,4*[fracvar]))
delxvec['U_static']=hstack((5*[0.01],0.01,0.01,10,4*[fracvar]))
delxvec['All_U_static']=hstack((5*[0.01],0.0001,0.0001,10,4*[fracvar]))
delxvec['All_U_static_FW_SI']=delxvec['All_U_static'].copy()
delxvec['All_U_static_FW_SI_plus']=delxvec['All_U_static'].copy()
fracvar=2
delxvec['All_U_static_FW_SI_wide']=hstack((5*[0.01],0.0001,0.0001,10,4*[fracvar]))
delxvec['U_static_FW_SI_wide']=hstack((5*[0.01],0.01,0.01,10,4*[fracvar]))


res={}
for dd in delxvec:
    delx=delxvec[dd]
    xinit=xg[dd]
    def CostFunction(x):
        m1=x[1]+x[4]*x[8]+x[5]*x[9]+x[6]*x[10]
        m2=x[0]+x[2]+x[3]+x[4]*(1-x[8])+x[5]*(1-x[9])+x[6]*(1-x[10])
        s1=(x[1]*S['PWS']+x[4]*x[8]*S['PWN']+x[5]*x[9]*S['FW']+x[6]*x[10]*S['SI'])/Snorm
        t1=(x[1]*T['PWS']+x[4]*x[8]*T['PWN']+x[5]*x[9]*T['FW']+x[6]*x[10]*T['SI']+x[7]*x[11])/Tnorm
        s2=(x[0]*S['AWS']+x[2]*S['DWS']+x[3]*S['AWN']+x[4]*(1-x[8])*S['PWN']+x[5]*(1-x[9])*S['FW']+x[6]*(1-x[10])*S['SI'])/Snorm
        t2=(x[0]*T['AWS']+x[2]*T['DWS']+x[3]*T['AWN']+x[4]*(1-x[8])*T['PWN']+x[5]*(1-x[9])*T['FW']+x[6]*(1-x[10])*T['SI'])/Tnorm
        J=(m1)**2 + (s1)**2 + (t1)**2 + (m2)**2 + (s2)**2 + (t2)**2 + sum((x-xinit)**2/delx**2)
        return J
    res[dd]=minimize(CostFunction,xinit,method='powell',options={'xtol':1e-8,'disp': True})
    barplot(dd)
