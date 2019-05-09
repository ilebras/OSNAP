from aux_funcs import *


##########

# Load the data from all 4 OOI moorings

ooi_pos=io.loadmat(datadir+'IrmingerSea/OOI_CTDs_allM_ilebras1905_pos.mat')

ooi_lat={}
ooi_lat['fla']=ooi_pos['lat'].flatten()[0][0][0][0]
ooi_lat['flb']=ooi_pos['lat'].flatten()[0][1][0][0]
ooi_lat['sfc']=ooi_pos['lat'].flatten()[0][2][0][0]
ooi_lat['prf']=ooi_pos['lat'].flatten()[0][3][0][0]


ooi_lon={}
ooi_lon['fla']=ooi_pos['lon'].flatten()[0][0][0][0]
ooi_lon['flb']=ooi_pos['lon'].flatten()[0][1][0][0]
ooi_lon['sfc']=ooi_pos['lon'].flatten()[0][2][0][0]
ooi_lon['prf']=ooi_pos['lon'].flatten()[0][3][0][0]

dat=io.loadmat(datadir+'IrmingerSea/OOI_CTDs_allM_ilebras1905_griddedonly.mat')

grd_date=date_from_matlab(dat['TIME'].flatten())

moorvec=['fla','flb','sfc','prf']

for moor in moorvec:
    dat[moor+'_PSAL'][isinf(dat[moor+'_PSAL'])]=NaN

ooi_all={}
for moor in moorvec:
    ooi_all[moor]=xr.Dataset({'salinity': (['prs','date'],dat[moor+'_PSAL']),
                         'temperature': (['prs','date'], dat[moor+'_TEMP']),
                         'density': (['prs','date'], dat[moor+'_DENS'])},
                         coords={'date': grd_date, 'prs': dat['Pi'].flatten(),'moor':moor})

ooi_xray=xr.concat([ooi_all[moor] for moor in moorvec],dim='moor').interpolate_na(dim='prs')

ooi_daily=ooi_xray.resample(date='1D').mean(dim='date')

def pltall(moor):
    # figure()
    (ooi_daily['temperature']).sel(moor=moor).T.plot(figsize=(12,3))
    ylim(2500,0)
    # figure()
    (ooi_daily['salinity']).sel(moor=moor).T.plot(figsize=(12,3))
    ylim(2500,0)
    # figure()
    (ooi_daily['density']-1000).sel(moor=moor).T.plot(figsize=(12,3))
    ylim(2500,0)


pltall('fla')

pltall('flb')

pltall('sfc')

SA=gsw.SA_from_SP(ooi_daily.salinity,ooi_daily.prs,ooi_lon['sfc'],ooi_lat['sfc'])
PT=gsw.pt0_from_t(SA,ooi_daily.temperature,ooi_daily.prs)
CT=gsw.CT_from_pt(SA,PT)
PDEN=gsw.sigma0(SA,CT)

ooi_daily

ooi_dnew=xr.Dataset({'sal': (['date','moor','prs'],ooi_daily.salinity.values[:,:3,:]),
                     'ptmp': (['date','moor','prs'],PT[:,:3,:]),
                     'pden': (['date','moor','prs'],PDEN[:,:3,:])},
                     coords={'date': ooi_daily.date.values, 'prs': dat['Pi'].flatten(),'moor':moorvec[:3]})

ooi_dnew

ooi_dnew=ooi_dnew.interpolate_na(dim='prs')

def pltden(moor):
    ooi_dnew['pden'].sel(moor=moor).T.plot(vmin=27,vmax=28,figsize=(12,3))
    axhline(90)
    ylim(2500,0)


pltden('fla')

def pltden_prf(moor):
    plot(ooi_dnew['pden'].sel(moor=moor).T,ooi_dnew.prs,'o')
    xlim()
    ylim(2500,0)
    axhline(90)

pltden_prf('sfc')


pltden_prf('fla')

pltden_prf('flb')

pltden('sfc')
pltden('fla')
pltden('flb')

def compden(pmin,pmax):
    plot(ooi_dnew['pden'].sel(moor='fla').where(ooi_dnew.prs<pmax).where(ooi_dnew.prs>pmin),ooi_dnew['pden'].sel(moor='sfc').where(ooi_dnew.prs<pmax).where(ooi_dnew.prs>pmin),'.');
    plot([27.2,27.8],[27.2,27.8],'k-')

def plt_dendiff(m1,m2):
    figure(figsize=(12,3))
    pcolor(ooi_daily.date,ooi_daily.prs,(ooi_dnew['pden'].sel(moor=m1)-ooi_dnew['pden'].sel(moor=m2)).T,cmap=cm.RdBu_r,vmin=-0.1,vmax=0.1)
    ylim(1000,0)
    colorbar()

plt_dendiff('sfc','fla')
plt_dendiff('sfc','flb')

plt_dendiff('fla','flb')

####################
## Load Femke's ML depth ANALYSIS
## NOTE:

dat=io.loadmat(datadir+'IrmingerSea/M4data_8_8_17_deJong0419.mat')

dat.keys()

# (PT =  potential temp, PD is potential dens ref 0, the rest should be more obvious.

# OOI keys:
# [('P', 'O'), ('T', 'O'), ('SP', 'O'), ('PT', 'O'), ('mtime', 'O'), ('O2', 'O'), ('lon', 'O'), ('lat', 'O'), ('SA', 'O'), ('CT', 'O'), ('PD', 'O'), ('N2', 'O'), ('P_mid', 'O'), ('in_funnel', 'O'), ('PV', 'O'), ('PV2', 'O')])

ooi_date = array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in dat['OOI'][0][0][4].flatten()])

ooi_depth=-gsw.z_from_p(dat['OOI'][0][0][0].flatten(),dat['OOI'][0][0][7][0])


contourf(log10(dat['OOI'][0][0][-2]),levels=51,cmap=cm.rainbow_r,vmin=-11.5,vmax=-10)
colorbar()

ooi=xr.Dataset({'salinity': (['prs','date'],dat['OOI'][0][0][2]),
                     'temperature': (['prs','date'], dat['OOI'][0][0][3]),
                     'potential density': (['prs','date'], dat['OOI'][0][0][10]),
                     'oxygen': (['prs','date'], dat['OOI'][0][0][5]),
                     'PV': (['prs_mid','date'], dat['OOI'][0][0][-2]),
                     'depth':(['prs'],ooi_depth)},
                     coords={'date': ooi_date, 'prs': dat['OOI'][0][0][0].flatten(),
                     'prs_mid': dat['OOI'][0][0][-4][:,0], 'lat':dat['OOI'][0][0][7].flatten(),'lon':dat['OOI'][0][0][6].flatten()})


pickle.dump(ooi,open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','wb'),protocol=2)

####################################################################
#################### save pickle of all 4 OOI moorings ############
####################################################################

ooi_dnew

ooi

prfonly=xr.Dataset({'sal': (['date','prs'],ooi.salinity.values.T),
                     'ptmp': (['date','prs'],ooi.temperature.values.T),
                     'pden': (['date','prs'],ooi['potential density'].values.T)},
                     coords={'date': ooi.date.values, 'prs': ooi.prs.values,'moor':'prf'})


ooi_all=xr.concat([ooi_dnew,prfonly],dim='moor')

ooi_all



####################################################################
#################### Mixed Layer calcs #############################
####################################################################

#LOCO is too sparse in the time period I'm interested in..

#LOCO 2 keys:
# [('P', 'O'), ('T', 'O'), ('SP', 'O'), ('mtime', 'O'), ('lon', 'O'), ('lat', 'O'), ('SP_or', 'O'), ('SP5', 'O'), ('T5', 'O'), ('SA', 'O'), ('CT', 'O'), ('PT', 'O'), ('PD', 'O'), ('CTsm', 'O'), ('SAsm', 'O'), ('PDsm', 'O'), ('N2', 'O'), ('P_mid', 'O'), ('in_funnel', 'O'), ('PV', 'O'), ('PV2', 'O')])

loco_date= array([datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366) for matlab_datenum in dat['L2'][0][0][3].flatten()])

loco=xr.Dataset({'salinity': (['prs','date'],dat['L2'][0][0][2]),
                     'potential temperature': (['prs','date'], dat['L2'][0][0][11]),
                     'potential density': (['prs','date'], dat['L2'][0][0][10]),
                     'PV': (['prs_mid','date'], dat['L2'][0][0][-2])},
                     coords={'date': loco_date, 'prs': dat['L2'][0][0][0].flatten(),
                     'prs_mid': dat['L2'][0][0][-4][:,0], 'lat':dat['L2'][0][0][5].flatten(),'lon':dat['L2'][0][0][4].flatten()})

loco.salinity.plot()
xlim(datetime.datetime(2014,1,1),datetime.datetime(2017,1,1))
