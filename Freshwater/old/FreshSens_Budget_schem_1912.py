from pylab import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/SimpleBudgets/'
##############################################################################
## Helper functions
def usum(Dict,leaveout='na'):
    Utot=0
    for ii in Dict:
        if ii==leaveout:
            Utot=Utot
        else:
            Utot=Utot+Dict[ii]
    return Utot

# What is the freshwater volume you would have to subtract from the transports to get this to work?
def SUsum(sdic,udic,leaveout='na'):
    SUtot=0
    for ii in S:
        if ii==leaveout:
            SUtot=SUtot
        else:
            SUtot=SUtot+sdic[ii]*udic[ii]
    return SUtot

##############################################################################
##############################################################################
#### First, get a mean budget based on literature, simplicity, and approx mass balance at both sections
##############################################################################
##############################################################################

U={}
S={}
##############################################################################
## OSNAP East values

U['AWS']=16
S['AWS']=35.2

U['DWS']=-13
S['DWS']=35

U['PWS']=-3
S['PWS']=34.2

#### Assume FramStrait + Barents Sea Opening are made up of PW and AW with different salinities with suffix N.
# Sum up northward AW through FS and BSO to get a single northward flowing AWN
##############################################################################
## sidebar: originally calculated these from fw and vol trans from literature
U['AWN']=-8
U['PWN']=+8

S['AWN']=35
S['PWN']=34.2
# S['PWN']=S['PWS'] #Override this to be the same salinity as PWS,

##############################################################################
## Freshwater and sea ice estimates
U['SI']=0.06
S['SI']=6

U['FW']=0.05 #mostly precip
usum(U)
S['FW']=0
U
S
U['AWS']=U['AWS']-(U['SI']+U['FW'])/4
U['PWN']=U['PWN']-(U['SI']+U['FW'])/4
U['AWN']=U['AWN']-(U['SI']+U['FW'])/4
U['PWS']=U['PWS']-(U['SI']+U['FW'])/4
usum(U)
U
SUsum(S,U)


ut = SUsum(S,U)/(S['AWS']-S['PWN']-S['PWS']+S['AWN'])

uts = SUsum(S,U)/(S['AWS']-S['PWS'])

utn = SUsum(S,U)/(S['AWN']-S['PWN'])

ut
uts
utn

U['AWS']=U['AWS']-ut
U['PWN']=U['PWN']+ut
U['PWS']=U['PWS']+ut
U['AWN']=U['AWN']-ut

U

SUsum(S,U)
usum(U)
U

U['PWS']+U['AWS']+U['DWS']
U['PWN']+U['AWN']

U.keys()
colvec=['red','darkblue','aquamarine','orange','c','purple','limegreen']


def barplot():
    f,[ax1,ax2,ax3,ax4]=subplots(4,1,figsize=(8,11))
    varlen=len(U.keys())
    ax1.bar(range(varlen),[U[kk] for kk in U],color=colvec)#,yerr=[U[kk]/2 for kk in U],capsize=10)
    ax2.bar(range(varlen),[S[kk] for kk in S],color=colvec)
    ax3.bar(range(varlen),[U[kk]*S[kk] for kk in S],color=colvec)
    ax4.bar(range(varlen),[U[kk]*(S[kk]-34.8)/34.8 for kk in S],color=colvec)
    ax2.set_ylim(33.5,35.3)
    ax2.set_yticks(arange(33.5,35.5,0.5))
    ax2.bar(range(varlen),[S[kk] for kk in S],color=colvec)
    ax2.axhline(34.8,color='grey')

    ax1.text(4.6,5,str(U['SI']),fontsize=24)
    ax1.text(5.6,5,str(U['FW']),fontsize=24)
    ax2.text(4.9,33.75,'6',fontsize=24)
    ax2.text(5.9,33.75,'0',fontsize=24)
    ax3.text(4.6,100,str(U['SI']*S['SI']),fontsize=24)
    ax3.text(5.9,100,'0',fontsize=24)
    for axx in [ax1,ax2,ax3,ax4]:
        axx.set_xticks(range(varlen))
        axx.set_xticklabels(U.keys())
        axx.axhline(0,color='grey')
    ax1.set_title('Simple Freshwater Budget Example')
    ax1.set_ylabel('Transport, U [Sv]')
    ax2.set_ylabel('Salinity, S')
    ax3.set_ylabel('Salinity budget term, U*S')
    ax4.set_ylabel('Using 34.8 reference salinity,\n U(S-34.8)/34.8')
    savefig(figdir+'FWBudget_Example_refcomp.png',bbox_inches='tight')

barplot()


########################################################################################################
########################################################################################################
############## Make a TS diagram that has all OSNAP and FS, BSO elements [just idealized ones for now]
T={}

#mean values from OSNAP breakdown
T['AWS']=8.8
T['PWS']=3.6
T['DWS']=3.3

## First, just do a quick Heat Transport conversion for the Tsubouchi values
rho0=1027
cp=3.987e3

HF_PWN=27e12
HF_AWN=92e12
U_PWN=-6.6 #in Tsubouchi sign convention... still confused about the signs in this conversion, but these results agree with theif Fig 12.
U_AWN=7.7
Tref=1.1

T['AWN']=(HF_AWN/(rho0*cp*U_AWN*1e6)+Tref)
T['AWN']
T['PWN']=(HF_PWN/(rho0*cp*U_PWN*1e6)+Tref)
T['PWN']

TSvec=list(U.keys())[:-2]

TSvec


salvec=linspace(31,36,103)
tmpvec=linspace(-3,16,103)
salmat,tmpmat=meshgrid(salvec,tmpvec)

import gsw
SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),0,65)

SA_vec_1000=gsw.SA_from_SP(salvec,1e3*ones(len(salvec)),0,65)

CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat=zeros((shape(salmat)))

for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(SA_vec[ii],CT_vec[jj])


def TSplot():
    [scatter(S[kk],T[kk],s=abs(U[kk])**2*4,c=colvec[ii],zorder=50,linewidth=3, label=TSvec[ii]) for ii,kk in enumerate(TSvec)]
    contour(salmat,tmpmat,pdenmat,levels=arange(25,29,0.2),colors='k')
    xlim(33.8,35.5)
    lgnd=legend(loc=(1.05,0.2),markerscale=1)
    for ii in range(5):
        lgnd.legendHandles[ii]._sizes = [40]

    ylim(-1,12)
    xlabel('salinity')
    ylabel('pot. temperature [$^\circ$C]')
    savefig(figdir+'TS_INOUT_Example.png',bbox_inches='tight')


U
PWS grew a lot... Re-examine how we got here... maybe slap some actual ranges on above fig?

TSplot()
# Hard to read, good to get a sense though. Strange thing right now is that PWN is saltier than PWS -- doesn't need to be that way for this example (prob not true.)
