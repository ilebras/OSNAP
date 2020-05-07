from firstfuncs_1618 import *

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


#############################################################
# Set up the base csae -- close the budget first (just tweaks I.C.s, really)
# will probably want to make this flexible enough to scan a range
#############################################################
AM={}
# base case only
AM['base']=array([[1,1,1,1,1,1,1,0],\
[S['PWS'],S['AWS'],S['DWS'],S['PWN'],S['AWN'],S['FW'],S['SI'],0],\
[T['PWS'],T['AWS'],T['DWS'],T['PWN'],T['AWN'],T['FW'],T['SI'],1]])

x0={}
x0['base']=[U['PWS'],U['AWS'],U['DWS'],U['PWN'],U['AWN'],U['FW'],U['SI'],Q]


dv={}
for kk in AM:
    dv[kk]=-AM[kk].dot(x0[kk])


Winv={}
Snorm=35
Tnorm=5
Winv['base']=diag([1,1/Snorm,1/Tnorm])

E={}
Qvar=5
Evec=hstack((WM['TRANS'].groupby('TIME.month').mean('TIME').std(dim='month').values,0.025,0.025,Qvar))
E['base']=diag(Evec)

### this is set up so you only have to define an E or Winv if it is different than the base case:
xsol_Ad={}
P={}
for zz in AM:
    if zz in Winv:
        if zz in E:
            Umat,D,VmatT=linalg.svd(Winv[zz].dot(AM[zz].dot(E[zz])))
            P[zz]=diag(E[zz]-E[zz].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E[zz].dot(AM[zz].T))+linalg.inv(Winv[zz])).dot(AM[zz].dot(E[zz])))))
        else:
            Umat,D,VmatT=linalg.svd(Winv[zz].dot(AM[zz].dot(E['base'])))
            P[zz]=diag(E['base']-E['base'].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E['base'].dot(AM[zz].T))+linalg.inv(Winv[zz])).dot(AM[zz].dot(E['base'])))))
    else:
        if zz in E:
            Umat,D,VmatT=linalg.svd(Winv['base'].dot(AM[zz].dot(E[zz])))
            P[zz]=diag(E[zz]-E[zz].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E[zz].dot(AM[zz].T))+linalg.inv(Winv['base'])).dot(AM[zz].dot(E[zz])))))
        else:
            Umat,D,VmatT=linalg.svd(Winv['base'].dot(AM[zz].dot(E['base'])))
            P[zz]=diag(E['base']-E['base'].dot(AM[zz].T.dot(linalg.inv(AM[zz].dot(E['base'].dot(AM[zz].T))+linalg.inv(Winv['base'])).dot(AM[zz].dot(E['base'])))))
    Lambda_inv = zeros((AM[zz].shape[0], AM[zz].shape[1])).T
    Lambda_inv[:AM[zz].shape[0], :AM[zz].shape[0]] = diag(1/D)
    if zz in Winv:
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv[zz].dot(dv[zz]))))
    else:
        xsol_prime=VmatT.T.dot(Lambda_inv.dot(Umat.T.dot(Winv['base'].dot(dv[zz]))))
    if zz in E:
        xsol_Ad[zz]=E[zz].dot(xsol_prime)
    else:
        xsol_Ad[zz]=E['base'].dot(xsol_prime)


#ditto for x0
res={}
for xx in AM:
    if xx in x0:
        res[xx]=x0[xx]+xsol_Ad[xx]
    else:
        res[xx]=x0['base']+xsol_Ad[xx]


def BudgetOut(x):
    mall=x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]
    return mall

shape(AM['base'])
shape(VmatT)
shape(D)
shape(Umat)
#############################################################
# Now add an inverse model for finding alpha,beta,gamma,delta proportions
#############################################################

AM_2={}
# base case only
xx='base'
AM_2[xx]=array([[Usol[xx]['PWN'],Usol[xx]['FW'],Usol[xx]['SI'],0],\
[Usol[xx]['PWN']*S['PWN'],0,Usol[xx]['SI']*S['SI'],0],\
[Usol[xx]['PWN']*T['PWN'],Usol[xx]['FW']*T['FW'],Usol[xx]['SI']*T['SI'],Usol[xx]['Q']]])


Usol={}
for xx in res:
    Usol[xx]=get_U_from_x(res[xx])

alpha=abs(Ubase['PWS']/U['PWN'])
beta=0.3
gamma=0.7
delta=0

x0_2={}
x0_2['base']=[alpha,beta,gamma,delta]

Winv_2={}
Snorm=35
Tnorm=5
Winv_2['base']=diag([1,1/Snorm,1/Tnorm])

E_2={}
E_2['base']=diag(4*[0.3])

dv_2={}
dv_2['base']=array([-Usol[xx]['PWS'],-Usol[xx]['PWS']*S['PWS'],-Usol[xx]['PWS']*T['PWS']])



### going to just make this simplest solver.
xsol_Ad_2={}
P_2={}
zz='base'
Umat_2,D_2,VmatT_2=linalg.svd(Winv_2[zz].dot(AM_2[zz].dot(E_2[zz])))
P_2[zz]=diag(E_2[zz]-E_2[zz].dot(AM_2[zz].T.dot(linalg.inv(AM_2[zz].dot(E_2[zz].dot(AM_2[zz].T))+linalg.inv(Winv_2[zz])).dot(AM_2[zz].dot(E_2[zz])))))
Lambda_inv_2 = zeros((AM_2[zz].shape[0], AM_2[zz].shape[1])).T
Lambda_inv_2[:AM_2[zz].shape[0], :AM_2[zz].shape[0]] = diag(1/D_2)
xsol_prime_2=VmatT_2.T.dot(Lambda_inv_2.dot(Umat_2.T.dot(Winv_2[zz].dot(dv_2[zz]))))
xsol_Ad_2[zz]=E_2[zz].dot(xsol_prime_2)

#ditto for x0
res_2={}
for xx in AM_2:
    if xx in x0_2:
        res_2[xx]=x0_2[xx]+xsol_Ad_2[xx]
    else:
        res_2[xx]=x0_2['base']+xsol_Ad_2[xx]

res_2

alpha

#############################################################
# Now add an inverse model for finding alpha,beta,gamma,delta proportions
#############################################################

coldic={'AWS':'red','DWS':'grey','PWS':'royalblue','PWN':'purple','AWN':'orange','SI':'cyan','FW':'k','Q':'limegreen'}


cp=3850
rhow=1000
cp*rhow*(Ux['Q']-Ubase['Q'])/1e6
