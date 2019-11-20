from pylab import *

# Here, I explore parameter space/ assumptions involved in constructing mass and salinity budgets (ignore temperature for now) between OSNAP and the Arctic
# Start off with constant best estimates from observations, then invert/fiddle?

#Sign convention: positive is INTO the box.

#### Try an even simpler version....

U={}
S={}

U['AWS']=15
S['AWS']=35.2

U['DWS']=-11
S['DWS']=35

#### Assume FramStrait + Barents Sea Opening are made up of PW and AW with different salinities with suffix N.
# Sum up northward AW through FS and BSO to get a single northward flowing AWN
##############################################################################
## sidebar: calculating these from fw and vol trans from literature
U['AWN']=-9
U['PWN']=+5

S['AWN']=35
S['PWN']=34.5 #Override this to be the same salinity as PWS,


# Let's see what the totals are right now

def usum(Dict,leaveout='na'):
    Utot=0
    for ii in Dict:
        if ii==leaveout:
            Utot=Utot
        else:
            Utot=Utot+Dict[ii]
    return Utot

Utot=usum(U)

Utot

# What is the freshwater volume you would have to subtract from the transports to get this to work?
def SUsum(sdic,udic,leaveout='na'):
    SUtot=0
    for ii in S:
        if ii==leaveout:
            SUtot=SUtot
        else:
            SUtot=SUtot+sdic[ii]*udic[ii]
    return SUtot

SUtot=SUsum(S,U)

U_tweak={}
for wm in S:
    U_tweak[wm]=-SUsum(S,U)/S[wm]

U

U_tweak

##############################################################################
# Try instead distributing it amongst water masses
uspread=SUsum(S,U)/usum(S)

uspread

U_dist={}
for wm in S:
    U_dist[wm]=U[wm]-uspread

U_dist
usum(U)

14/1000

usum(U_dist)

SUsum(S,U_dist)
##############################################################################
# Try instead adding and subtracting same transport from AWS and PWN:

ut = SUsum(S,U)/(S['AWS']-S['PWN'])

##############################################################################
# How much would you have to change each salinity by?
S_tweak = {}
for wm in S:
    S_tweak[wm] = S[wm]  - SUsum(S,U)/U[wm]

S_tweak
SUsum(S,U)
6*0.06

for wm in S:
    print(wm)
    print(U[wm]*S[wm])

##############################################################################
##############################################################################
# Add sea ice and see what it does:
##############################################################################
##############################################################################
S_SI=S.copy()
U_SI=U.copy()

U_SI['SI']=0.06
S_SI['SI']=6


for wm in U:
    U_SI[wm]=U[wm]-0.06/4

# Make no change to sea ice transport:
uspread_nosi=SUsum(S_SI,U_SI)/usum(S_SI,leaveout='SI')
uspread_nosi

U_dist_noSI={}
for wm in S:
    U_dist_noSI[wm]=U_SI[wm]-uspread_nosi

U_dist_noSI
usum(U_dist_noSI)+U_SI['SI']

# Include sea ice in the change:
uspread_si=SUsum(S_SI,U_SI)/usum(S_SI)

uspread_si
U_dist_SI={}
for wm in S_SI:
    U_dist_SI[wm]=U_SI[wm]-uspread_si

usum(U_dist_SI)


# What about taking out 
ut_SI = SUsum(S,U)/(S['AWS']-S['PWN'])
