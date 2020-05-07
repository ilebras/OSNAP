from firstfuncs_1618 import *
figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Freshwater/NorESM_testing/'
################################################################################################################################
#sigma2 levels copied over from rho.txt:
sig2lev=array([27.2200000000000, 27.7200000000000,28.2020000000000,28.6810000000000,29.1579999999999,29.6320000000001,
30.1020000000001,30.5670000000000,31.0260000000001,31.4770000000001,31.9200000000001,
32.3520000000001,32.7719999999999,33.1759999999999,33.5640000000001,33.9320000000000,34.2790000000000,
34.6020000000001,34.9000000000001,35.1720000000000,35.4169999999999,35.6369999999999,35.8320000000001,
36.0029999999999,36.1530000000000,36.2840000000001,36.3979999999999,36.4970000000001,
36.5840000000001,36.6600000000001,36.7280000000001,36.7890000000000,36.8430000000001,
36.8930000000000,36.9390000000001,36.9820000000000,37.0219999999999,37.0599999999999,
37.0960000000000,37.1310000000001,37.1659999999999,37.1990000000001,37.2310000000000,
37.2639999999999,37.2950000000001,37.3270000000000,37.3580000000000,37.3879999999999,
37.4190000000001,37.4500000000000,37.4800000000000,37.5799999999999,37.8000000000000])

plot(sig2lev,'.')

################################################################################################################################
################################################################################################################################
###########################################      LOAD     #######################################################
################################################################################################################################
################################################################################################################################

#### load noresm and force formatting to be the same
def load_noresm(whichone):
    noresm2=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*_sigma_2010*.nc')[0])
    noresm1=xr.open_dataset(glob.glob(datadir+'NorESM/*'+whichone+'*_sigma_2000*.nc')[0])
    noresm=xr.concat([noresm1,noresm2],dim='time',data_vars='minimal')
    noresm=noresm.rename({'time': 'TIME','dz':'DZ','cell':'LONGITUDE'})
    noresm=noresm.rename({'salt':'PSAL','temp':'PTMP'})
    noresm=noresm.assign_coords(LONGITUDE=(noresm1.lon.values))
    noresm=noresm.assign_coords(LATITUDE=(noresm1.lat.values))
    noresm=noresm.drop('lat').drop('lon')
    noresm['DEPTH']=noresm['deltaz'].cumsum(dim='DZ')
    startyear=2000
    startmonth=1
    endyear=2018
    endmonth=12
    noresm['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return noresm

osnap=load_noresm('OSNAP')
fs=load_noresm('FS')
bso=load_noresm('BSO')
ns=load_noresm('NS')

def add_PDEN(xray):
    PRES_out=gsw.p_from_z(-xray['DEPTH'],60)
    SA_out=NaN*xray['PSAL'].copy()
    for ii in range(shape(PRES_out)[1]):
        for jj,ll in enumerate(xray.LONGITUDE):
            SA_out[:,ii,jj]=gsw.SA_from_SP(xray['PSAL'][:,ii,jj],PRES_out[:,ii,jj],ll,xray.LATITUDE[jj])
    if 'PTMP' in list(xray.data_vars):
        PT_out=xray['PTMP']
    else:
        PT_out=gsw.pt0_from_t(SA_out,xray['TEMP'],PRES_out)
    CT_out=gsw.CT_from_pt(SA_out,PT_out)
    PD_out=gsw.sigma0(SA_out,CT_out)
    xray['PDEN']=(('TIME','DZ','LONGITUDE'),PD_out)

add_PDEN(osnap)
add_PDEN(fs)
add_PDEN(bso)
add_PDEN(ns)
################################################################################################################################
################################################################################################################################
###################################     Use velocity as given by model to get TRANS!   #######################################################
################################################################################################################################
################################################################################################################################

################################################################################################################################
################################################################################################################################
#####################     Calc perpendicular velocity and Transport   #######################################################
################################################################################################################################
################################################################################################################################
bso_uv=array([(1,0),(0,	1),(1,	0),(0,	1),(1,	0),(1,	0),(0,	1),(1,	0),(0,	1),(1,	0),(1,	0),(0,	1),(1,	0),(0,	1),(1,	0),(1,	0),(0,	1),(1,	0),(0,	1),(1,	0),(1,	0),(0,	1),(1,	0),(1,	0)])
fs_uv=array(shape(fs.LONGITUDE)[0]*[(0,1)])
osnap_uv=array([(0,1),(0,1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(1,	1),(0,	1),(0,	1),(1,	1),(0,	1),(0,	1),(1,	1),(0,	1),(0,	1),(1,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1),(0,	1)])

osnap.deltaz.mean(dim='TIME').sum(dim='DZ').plot()

osnap.PSAL.isel(LONGITUDE=0).plot()

def add_uv_ind(xrw,uv_var):
    xrw['u_ind']=(['LONGITUDE'],uv_var[:,0])
    xrw['v_ind']=(['LONGITUDE'],uv_var[:,1])
    return xrw

osnap=add_uv_ind(osnap,osnap_uv)
bso=add_uv_ind(bso,bso_uv)
fs=add_uv_ind(fs,fs_uv)

def get_TRANS(xrw):
    xrw['TRANS']=(xrw['uvel']*xrw['u_ind']*xrw['dy']+xrw['vvel']*xrw['v_ind']*xrw['dx'])*xrw['deltaz']/1e6
    return xrw

for ff in [osnap,fs,bso]:
    ff=get_TRANS(ff)

################################################################################################################################
###########################  Compare to model diagnostics   #######################################################
################################################################################################################################
startyear=2000
startmonth=1
endyear=2018
endmonth=12

def load_allyears(varlab):
    tmp1=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_200001-200912.nc')[0])
    print(tmp1)
    tmp2=xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_201001-201812.nc')[0])
    print(tmp2)
    varout=xr.concat([tmp1,tmp2],dim='time')
    varout=varout.rename({'time':'TIME'})
    varout['TIME']=array([datetime.datetime(m//12, m%12+1, 15) for m in range(startyear*12+startmonth-1, endyear*12+endmonth)])
    return varout

varlab='volumetransports'
xr.open_dataset(glob.glob(datadir+'NorESM/NorESM2*'+varlab+'_200001-200912.nc')[0])

vts=load_allyears('volumetransports')
hts=load_allyears('heattransports')
sts=load_allyears('salinitytransports')

hts

def comp_voltrans(var,vtsname):
    figure(figsize=(12,3))
    vts[vtsname].plot(label='diagnostic = '+str(int(vts[vtsname].mean()*100)/100)+' Sv')
    var['TRANS'].sum(dim='DZ').sum(dim='LONGITUDE').plot(label='calculated from fields = '+str(int(var['TRANS'].sum(dim='DZ').sum(dim='LONGITUDE').mean()*100)/100)+' Sv')
    ylabel('Volume transport [Sv]')
    legend(loc=(1.05,0.2))
    title(vtsname[7:])
    savefig(figdir+'NetVol_'+vtsname[7:]+'.png',bbox_inches='tight')


comp_voltrans(osnap,'net_vt_OSNAP')
comp_voltrans(bso,'net_vt_BSO')
comp_voltrans(fs,'net_vt_FS')
cp=3850
rhow=1025

def comp_heat_trans(var,vtsname):
    figure(figsize=(12,3))
    (hts[vtsname]/1e12).plot(label='diagnostic = '+str(int(hts[vtsname].mean()/1e12))+' TW')
    (var['PTMP']*var['TRANS']*cp*rhow/1e6).sum(dim='DZ').sum(dim='LONGITUDE').plot(label='calculated from fields = '+str(int((var['TRANS']*var['PTMP']).sum(dim='DZ').sum(dim='LONGITUDE').mean()*cp*rhow/1e6))+' TW')
    ylabel('Temp transport [TW]')
    legend(loc=(1.05,0.2))
    title(vtsname[7:])
    savefig(figdir+'NetHeat_'+vtsname[7:]+'.png',bbox_inches='tight')

comp_heat_trans(osnap,'net_ht_OSNAP')
comp_heat_trans(fs,'net_ht_FS')
comp_heat_trans(bso,'net_ht_BSO')

###################################################################################################################
############################### check out the total convergence ##################################################
###################################################################################################################

latvol={}
latvol['diag']=vts.net_vt_OSNAP-vts.net_vt_BSO-vts.net_vt_FS
latvol['fields']=osnap.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')-fs.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')-bso.TRANS.sum(dim='DZ').sum(dim='LONGITUDE')

latheat={}
latheat['diag']=(hts.net_ht_OSNAP-hts.net_ht_BSO-hts.net_ht_FS)/1e12
latheat['fields']=((osnap.TRANS*osnap.PTMP).sum(dim='DZ').sum(dim='LONGITUDE')-(fs.TRANS*fs.PTMP).sum(dim='DZ').sum(dim='LONGITUDE')-(bso.TRANS*bso.PTMP).sum(dim='DZ').sum(dim='LONGITUDE'))*cp*rhow/1e6

latsal={}
latsal['diag']=sts.net_st_OSNAP-sts.net_st_BSO-sts.net_st_FS
latsal['fields']=((osnap.TRANS*osnap.PSAL).sum(dim='DZ').sum(dim='LONGITUDE')-(fs.TRANS*fs.PSAL).sum(dim='DZ').sum(dim='LONGITUDE')-(bso.TRANS*bso.PSAL).sum(dim='DZ').sum(dim='LONGITUDE'))

def comp_tot(var,tit,tit2):
    figure(figsize=(12,3))
    var['diag'].plot(label='diagnostics = '+str(var['diag'].mean().values))
    var['fields'].plot(label='fields = '+str(var['fields'].mean().values))
    legend(loc=(1.01,0.3))
    ylabel(tit)
    title('Total Convergence')
    savefig(figdir+'NetConv_'+tit2+'.png',bbox_inches='tight')

comp_tot(latvol,'Transport [Sv]','voltrans')

comp_tot(latheat,'Heat transport [TW]','heat')

comp_tot(latsal,'Salinity convervence [S x Sv]','salinity')

osnap.to_netcdf(datadir+'NorESM/NorESM_osnap_xray_18yrs_2004_sigma_UVcorr.nc','w')
fs.to_netcdf(datadir+'NorESM/NorESM_fs_xray_18yrs_2004_sigma_UVcorr.nc','w')
bso.to_netcdf(datadir+'NorESM/NorESM_bso_xray_18yrs_2004_sigma_UVcorr.nc','w')
