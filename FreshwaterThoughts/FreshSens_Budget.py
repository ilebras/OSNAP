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
S['PWN']=34.5
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

##############################################################################
##############################################################################
#### Next, ask the question, how does this change if FW doubles?
##############################################################################
##############################################################################
U_FW=U.copy()

fwmult=4
U_FW['FW']=fwmult*U['FW']

U['FW']

##############################################################################
#### Mass balance has to be satisfied:
for wm in U:
    if wm!='FW':
        if wm !='SI':
            U_FW[wm]=U[wm]-(fwmult-1)*U['FW']/5
U_FW

usum(U_FW)
usum(U)

SUsum(S,U_FW)


##############################################################################
#### Allow both north and south to compensate:

alpha=arange(0,1.1,0.1)
len(alpha)
SUsum(S,U_FW)

uta=[SUsum(S,U_FW)/(aa*(S['PWN']-S['AWN'])+(1-aa)*(S['PWS']-S['DWS'])) for aa in alpha]

U_FWA={}
for wm in U_FW:
    U_FWA[wm]=U_FW[wm]*ones(len(alpha))

U_FWA

alpha*uta

U_FWA['PWN']=array(U_FW['PWN']-alpha*uta)
U_FWA['AWN']=array(U_FW['AWN']+alpha*uta)

U_FWA['PWS']=array(U_FW['PWS']-(1-alpha)*uta)
U_FWA['DWS']=array(U_FW['DWS']+(1-alpha)*uta)

col={}
for ii,wm in enumerate(U_FWA):
    col[wm]='C'+str(ii)



def plot_FW_volresp(alpha,U_FWA,tit='na'):
    figure(figsize=(12,6))
    subplot(211)
    [axhline(U_FW[wm],color=col[wm],linestyle='--',label='') for wm in ['PWN','AWN']]
    [plot(alpha,U_FWA[wm],color=col[wm],marker='o',label=wm) for wm in ['PWN','AWN']]
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    legend(loc=(1.05,0.5))
    # xlim(-0.01,1.01)
    subplot(212)
    [axhline(U_FW[wm],color=col[wm],linestyle='--',label='') for wm in ['PWS','AWS','DWS']]
    [plot(alpha,U_FWA[wm],color=col[wm],marker='o',label=wm) for wm in ['PWS','AWS','DWS']]
    axhline(0,color='k')
    ylabel('Transport [Sv]')
    legend(loc=(1.05,0.5))
    # xlim(-0.01,1.01)
    if tit=='na':
        savefig(figdir+'VolTransResp_'+str(fwmult)+'xFW.png',bbox_inches='tight')
    else:
        savefig(figdir+'VolTransResp_'+tit+'.png',bbox_inches='tight')

plot_FW_volresp(alpha,U_FWA)

fwmult
##############################################################################
#### I actually think a more likely scenario is that the salinities would change... Transports being set by winds....
# However, there are no constraints on salinity of each water mass.
# Start off by assuming that only outflowing salinities would change (DWS,PWS,AWN)
fwmult
U_FW

# First, assume that only one water mass changes its salinity:
stweak={}
for wm in ['AWN','PWS','DWS']:
    stweak[wm]=SUsum(S,U_FW)/U_FW[wm]
svec_pws=linspace(0,stweak['PWS'],10)
svec_dws=linspace(0,stweak['DWS'],10)

smat_pws,smat_dws=meshgrid(svec_pws,svec_dws)

S['AWN']

smat_awn=(SUsum(S,U_FW)-smat_pws*U['PWS']-smat_dws*U['DWS'])/U['AWN']

smat_pws
smat_dws

smat_awn

def plot_salout():
    contourf(S['PWS']-svec_pws,S['DWS']-svec_dws,S['AWN']-smat_awn)
    colorbar(label='Atlantic Water salinity')
    contour(S['PWS']-svec_pws,S['DWS']-svec_dws,S['AWN']-smat_awn,levels=[S['AWN']],colors='k',linewidths=3)
    xlabel('Polar Water salinity')
    ylabel('Deep Water salinity')


plot_salout()

def plot_salEM():
    #plot the three end-member solutions
    plot([S['PWS'],S['PWS'],S['PWS']-stweak['PWS']],color=col['PWS'],label='PWS')
    plot([S['DWS'],S['DWS']-stweak['DWS'],S['DWS']],color=col['DWS'],label='DWS')
    plot([S['AWN']-stweak['AWN'],S['AWN'],S['AWN']],color=col['AWN'],label='AWN')
    legend()
    ylabel('Salinity')
    xticks(range(3))
    gca().set_xticklabels(['AWN response','DWS response','PWS response'])

plot_salEM()


##############################################################################
##############################################################################
#### What about if PW inflow gets fresher??
##############################################################################
##############################################################################
S
# Let's assume that, as before only AWN, PWS, and DWS react to a change in PWN salinity:
SPW_vec=linspace(33,S['PWN'],10)
S_PW={}
for wm in S:
    S_PW[wm]=S[wm]*ones(len(SPW_vec))

S_PW['PWN']=SPW_vec

S_PW

# First, assume that only one water mass changes its salinity:
sresp={}
for wm in ['AWN','PWS','DWS']:
    sresp[wm]=SUsum(S_PW,U)/U[wm]

def plot_salout_PW(ii):
    smat_pws,smat_dws=meshgrid(linspace(0,sresp['PWS'][ii],15),linspace(0,sresp['DWS'][ii],15))
    smat_awn=(SUsum(S,U_FW)-smat_pws*U['PWS']-smat_dws*U['DWS'])/U['AWN']
    contourf(S['PWS']-linspace(0,sresp['PWS'][ii],15),S['DWS']-linspace(0,sresp['DWS'][ii],15),S['AWN']-smat_awn)
    colorbar(label='Atlantic Water salinity')
    contour(S['PWS']-linspace(0,sresp['PWS'][ii],15),S['DWS']-linspace(0,sresp['DWS'][ii],15),S['AWN']-smat_awn,levels=[S['AWN']],colors='k',linewidths=3)
    xlabel('Polar Water salinity')
    ylabel('Deep Water salinity')
    title('Change inflowing Polar Water salinity to '+str(S_PW['PWN'][ii]))


plot_salout_PW(0)

plot_salout_PW(6)
sresp

##############################################################################
##############################################################################
#### What if instead PW inflow increases (and AW outflow increases by same amount)
#### How do the southern outflow transports have to adjust?
##############################################################################
#############################################################################

# Vary PWN from 3Sv to 12 Sv
volvar=arange(-3,6,1)
# And have AWN adjust accordingly
U_Nvar={}
for wm in U:
    U_Nvar[wm]=U[wm]*ones(len(volvar))
U_Nvar['PWN']=U['PWN']+volvar
U_Nvar['AWN']=U['AWN']-volvar

U_PWSDWS=U_Nvar.copy()
PWSDWS=SUsum(U_Nvar,S)/(S['DWS']-S['PWS'])
U_PWSDWS['DWS']=U['DWS']-PWSDWS
U_PWSDWS['PWS']=U['PWS']+PWSDWS
plot_FW_volresp(range(len(volvar)),U_PWSDWS,'NorthVar_PWSDWSresp')

U_PWSAWS=U_Nvar.copy()
PWSAWS=SUsum(U_Nvar,S)/(S['AWS']-S['PWS'])
U_PWSAWS['AWS']=U['AWS']-PWSAWS
U_PWSAWS['PWS']=U['PWS']+PWSAWS
SUsum(S,U_PWSAWS)

plot_FW_volresp(range(len(volvar)),U_PWSAWS,'NorthVar_PWSAWSresp')

U_DWSAWS=U_Nvar.copy()
DWSAWS=SUsum(U_Nvar,S)/(S['AWS']-S['DWS'])
U_DWSAWS['AWS']=U['AWS']-DWSAWS
U_DWSAWS['DWS']=U['DWS']+DWSAWS
SUsum(S,U_DWSAWS)

plot_FW_volresp(range(len(volvar)),U_DWSAWS,'NorthVar_DWSAWSresp')
