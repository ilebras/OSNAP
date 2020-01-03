from pylab import *

# Here, I explore parameter space/ assumptions involved in constructing mass and salinity budgets (ignore temperature for now) between OSNAP and the Arctic
# Start off with constant best estimates from observations, then invert/fiddle?

#Sign convention: positive is INTO the box.

#### Assume OSNAP is made up of 3 water masses (AW, PW, and DW) using the suffix S to convey that they are at the southern end of the domain.
# This first pass is from 2014-16 (21 month mean) OSNAP East estimate
# U_ is transport in Sv. S_ is salinity.
U={}
S={}

U['AWS']=16
S['AWS']=35.2

U['PWS']=-3
S['PWS']=34.5

U['DWS']=-11
S['DWS']=35

#### Assume FramStrait + Barents Sea Opening are made up of PW and AW with different salinities with suffix N.
# Sum up northward AW through FS and BSO to get a single northward flowing AWN
##############################################################################
## sidebar: calculating these from fw and vol trans from literature
#
# U['AWN']=-9
# U['PWN']=+8

U['AWN']=-2.3-5.4
U['PWN']=6.6 #(from Tsubouchi )

sref=34.67

FWT_AWN=(58+15)/1000
FWT_PWN=(21)/1000

S['AWN']=(-FWT_AWN/U['AWN']+1)*sref
S['PWN']=(-FWT_PWN/U['PWN']+1)*sref

S['AWN']

S['AWN']
S['PWN']

S['PWN']=34.5 #Override this to be the same salinity as PWS,

### Another worry here is that the "DW" sigma1 isopycnal in Fram Strait is well below anything we see at OSNAP, maybe our diagrams need re-drawing...
## OR its just recirculating and it doesn't matter.

#### Include a sea ice component with a salinity of 6
U['SI']= 0.06 #From Tsubouchi et al. 2018 (and assuming it all comes in Fram Strait and doesn't leave.)
S['SI']= 6

#### Finally, sources of 0 salinity FW such as precip and Greenland
U['FW'] = 17/1000 + 0.1
# Bamber et al. 2018 estimate: 17mSv from Greenland to Eastern Greenland (Def in the noise of below precip estimate.=)
# Schanze et al. 2010 estimated from their Figure 6
#I'll be testing sensitivity to this...


# Let's see what the totals are right now

def usum(Dict,leaveout='na'):
    Utot=0
    for ii in Dict:
        if ii==leaveout:
            Utot=Utot
        Utot=Utot+Dict[ii]
    return Utot

Utot=usum(U)

Utot

def SUsum(sdic,udic,leaveout='na'):
    SUtot=0
    for ii in S:
        if ii==leaveout:
            SUtot=SUtot
        else:
            SUtot=SUtot+sdic[ii]*udic[ii]
    return SUtot

SUtot=SUsum(S,U)

# Try out a few different solutions, and plot and record them in TS space.

##############################################################################
# Try distributing the volume imbalance evenly (except for FW and seaice)

nus=len(U.keys())
ucorr=Utot/(nus-2)

U_even={}
for ii in U:
    if (ii=='FW') | (ii=='SI'):
        U_even[ii]=U[ii]
    else:
        U_even[ii]=U[ii]-ucorr

U_even

Ute=usum(U_even)
SUte=SUsum(S,U_even)
##############################################################################
# Ask the question, what change in salinity would be required in each WM to make the system balance?
S_tweak={}
for wm in S:
    S_tweak[wm]=-SUsum(S,U_even,leaveout=wm)/U[wm]

U_even

S_tweak

##############################################################################
# None of those numbers are realistic. What if you require each WM's salinity to change by a set amount?
# I don't think you can do this with a mass balanced system. (because SUM(u)=0, can't divide by zero...)
scorr=SUsum(S,U)/usum(U)
scorr
# ok this number is crazy, might have done something wrong.


##############################################################################
# Ask the question, what transport change would be required in each WM to make the salinity budget balance?
# Note that the remaining imbalance would have to be compensated by FW...
U_tweak={}
for wm in S:
    U_tweak[wm]=-SUsum(S,U_even,leaveout=wm)/S[wm]

U_tweak

for wm in S:
    print(wm,U[wm]-U_tweak[wm])
Actually I'm quite confused about what I am doing here...

##############################################################################
# OK, that's about 1 Sv each...
# Try instead distributing it amongst water masses (except seaice!)
uspread=SUsum(S,U_even)/usum(S,leaveout='SI')

uspread

U_dist={}
for wm in S:
    U_dist[wm]=U_even[wm]-uspread

U_dist
usum(U_dist)

Or maybe this is OK now...?


# ###### From Lars H. Smedsrund (NOTE: Isopycnal of maximum overturning changes a lot between GSR and OSNAP due to entrainment)
# # Mean Values on GSR
# U_AWR = 8.5; # Atlantic Water Volume Flux [Sv]
# U_PWR = -2.5; # Polar Water
# U_OWR = -6.0; # Overflow Water
#
# T_AWR = 8.0; # AW temperature - should vary [deg C]
# T_PWR = -1.5; # PW temp - may vary
# T_OWR = 0.5; # OW temp - may vary
#
# S_AWR = 35.2;  # AW salt - fixed [psu]
# S_PWR = 34.0;  # PW salt - fixed [psu]
# S_OWR = 34.85; # OW salt - fixed [psu]
