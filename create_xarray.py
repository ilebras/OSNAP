from aux_funcs import *

#################################################################################
### Create padded arrays so that all moorings are on the same grid
#################################################################################

# Manually input the min and max pressures/dates to pad the arrays to
pmin=0
pmax=2112

dmin=dt.datetime(2014,8,17)
dmax=dt.datetime(2016,7,28)

plen=int((pmax-pmin)/2+1)

hrlen=int(divmod((dmax-dmin).total_seconds(),60*60*24)[0])+1

def padmax(field_less,field_more):
    origmax=max(field_less.index)
    missvc=arange(origmax+2,max(field_more.index)+1,2)
    missdf=pd.DataFrame(index=missvc,columns=field_less.columns,
                            data=tile(field_less[field_less.index==origmax].values[0],[len(missvc),1]))
    field_final=field_less.append(missdf)
    return field_final

def assmat(field):
    matone=field[field.columns[(field.columns>=dmin)&(field.columns<=dmax)]]
    return matone.values


salmat=999*ma.ones((8,plen,hrlen))
tmpmat=999*ma.ones((8,plen,hrlen))
pdenmat=999*ma.ones((8,plen,hrlen))
uacrossmat=999*ma.ones((8,plen,hrlen))
ualongmat=999*ma.ones((8,plen,hrlen))

vminlist=[]
sminlist=[]
vmaxlist=[]
smaxlist=[]
dminlist=[]
dmaxlist=[]

for moor in range(1,9):
    sal,tmp,pden=pickle.load(open('../pickles/TSinterp/CF'+str(moor)+
                                  '_saltmpinterp_ptcorr_notid.pickle','rb'))

    dminlist=hstack((dminlist,min(sal.columns)))
    dmaxlist=hstack((dmaxlist,max(sal.columns)))

    sminlist=hstack((sminlist,min(sal.index)))
    smaxlist=hstack((smaxlist,max(sal.index)))

    u,v=pd.read_pickle(open('../pickles/VELinterp/CF'+str(moor)+'_uvinterp.pickle','rb'))
    u_across=u*cos(theta)+v*sin(theta)
    u_along=-u*sin(theta)+v*cos(theta)


    dminlist=hstack((dminlist,min(v.columns)))
    dmaxlist=hstack((dmaxlist,max(v.columns)))

    vminlist=hstack((vminlist,min(v.index)))
    vmaxlist=hstack((vmaxlist,max(v.index)))

#     extend bottommost value to match the deepest one (be it sal/tmp/den or vel)
    if max(sal.index)>max(u_across.index):
        plenmoor=len(sal.index)
        u_across=padmax(u_across,sal)
        u_along=padmax(u_along,sal)

    else:
        plenmoor=len(u_across.index)
        sal=padmax(sal,u_across)
        tmp=padmax(tmp,u_across)
        pden=padmax(pden,u_across)

    salmat[moor-1,:plenmoor,:]=assmat(sal)
    tmpmat[moor-1,:plenmoor,:]=assmat(tmp)
    pdenmat[moor-1,:plenmoor,:]=assmat(pden)
    uacrossmat[moor-1,:plenmoor,:]=assmat(u_across)
    ualongmat[moor-1,:plenmoor,:]=assmat(u_along)


def fillmask(field):
    field[field==999]=ma.masked
    return field


salmat=fillmask(salmat)
tmpmat=fillmask(tmpmat)
pdenmat=fillmask(pdenmat)
uacrossmat=fillmask(uacrossmat)
ualongmat=fillmask(ualongmat)

#################################################################################
# Create an xarray
#################################################################################

distvecint=[int(dd) for dd in distvec]

daily=xr.Dataset({'temperature': (['distance', 'depth', 'date'],  tmpmat),
                'salinity': (['distance', 'depth', 'date'],  salmat),
                'potential density': (['distance', 'depth', 'date'],  pdenmat),
               'across track velocity': (['distance', 'depth', 'date'],  uacrossmat),
               'along track velocity': (['distance', 'depth', 'date'],  ualongmat),},
                coords={'distance': distvecint,
                        'depth': sw.dpth(sal.index.values,60),
                        'date': sal.columns[(sal.columns>=dmin)&(sal.columns<=dmax)].values})

pickle.dump(daily,open('../pickles/CF_xarray_notid_newtheta.pickle','wb'))
