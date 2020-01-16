from firstfuncs_1618 import *

figdir='/home/isabela/Documents/projects/OSNAP/figures_1418_merged/'

# In this script I want to assess the current status and lay out the necessary steps to acheive uniformity.

# What does M1 2016-2018 data look like? Might be what I'm trying to match.
m1dir18=glob.glob(datadir+'OSNAP2018recovery/M1_netcdf/*MCAT*')
m1dir18
m1_18=xr.open_dataset(m1dir18[0])

def quickpresplot(xray):
    figure(figsize=(14,6))
    plot(xray.TIME,xray.PRES);

quickpresplot(m1_18)

m1_18_daily=m1_18.resample(TIME='1D').mean(dim='TIME')

add_SA_CT_PT(m1_18_daily)

m1_18_daily.to_netcdf(datadir+'OSNAP2018recovery/Daily_netcdf/M1_mcat_2018recovery_daily.nc','w',format='netCDF4')


#Note: M1, 2014-2016 is two deployments -- merge their netcdfs
m1dir16=glob.glob(datadir+'OSNAP2016recovery/M1_netcdf/*MCAT*')

m1_16=xr.open_dataset(m1dir16[0])
m1_15=xr.open_dataset(m1dir16[1])

m1_16_merged=xr.concat([m1_16,m1_15],dim='DEPTH')

quickpresplot(m1_16_merged)

m1_16_daily=m1_16_merged.resample(TIME='1D').mean(dim='TIME')


m1_16_daily['LATITUDE']=59.9
m1_16_daily['LONGITUDE']=-41.125

add_SA_CT_PT(m1_16_daily)


m1_16_daily.to_netcdf(datadir+'OSNAP2016recovery/Daily_netcdf/M1_mcat_2016recovery_daily.nc','w',format='netCDF4')

figure(figsize=(14,6))
plot(m1_16_daily.TIME,m1_16_daily.TEMP);

# How about CF data?

# 2016-2018 data is dip-calibrated and in netcdf format.
#Here's my data merging strategy so far:
# From FirstLook_TSdata_Grid_2018reco.py
CF_18_daily={}
for ii in range(1,8):
    CF_18_daily[ii]=xr.open_dataset(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(ii)+'_mcat_2018recovery_daily.nc')


#Get 2014-2016 data there:
[date_all,month_all,prs_all,sal_all,tmp_all]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/TSdailydic/TS_daily_dic_wJHIL.pickle','rb'))

#confirm that length/start date lines up
for ii in range(1,8):
    for kk in prs_all[ii]:
        print(ii,shape(prs_all[ii][kk]),date_all[ii][kk][0])
#confirmed this is true (but only within moorings, which is fine)

#make sal, ptmp, and pressure matrices
prs_mat={}
sal_mat={}
ptmp_mat={}
date_vec={}
for ii in range(2,8):
    dpths=list(date_all[ii].keys())
    print(dpths)
    date_vec[ii]=date_all[ii][dpths[-1]]
    prs_mat[ii]=zeros((len(date_vec[ii]),len(dpths)))
    sal_mat[ii]=prs_mat[ii].copy()
    ptmp_mat[ii]=prs_mat[ii].copy()
    for kk in range(len(dpths)):
        tmp_len=len(prs_all[ii][dpths[kk]])
        prs_mat[ii][:tmp_len,kk]=prs_all[ii][dpths[kk]]
        sal_mat[ii][:tmp_len,kk]=sal_all[ii][dpths[kk]]
        ptmp_mat[ii][:tmp_len,kk]=tmp_all[ii][dpths[kk]]


#make an xarray for each mooring
CF_16_daily={}
for ii in range(2,8):
    CF_16_daily[ii]=xr.Dataset({'PRES': (['TIME','DEPTH'], prs_mat[ii] ),
                'PSAL': (['TIME','DEPTH'], sal_mat[ii]),
                'PTMP': (['TIME','DEPTH'],  ptmp_mat[ii]),},
                coords={'TIME': date_vec[ii],
                        'DEPTH': list(prs_all[ii].keys()),
                        'LATITUDE': CFlat[ii-1],
                        'LONGITUDE': CFlon[ii-1]})
    add_SA_CT_PT(CF_16_daily[ii])
    CF_16_daily[ii].to_netcdf(datadir+'OSNAP2016recovery/Daily_netcdf/CF'+str(ii)+'_mcat_2016recovery_daily.nc','w',format='netCDF4')

## note: I had some issues with CF1, so I took it out for now. The already produced xarray should be fine (and in fact I'm not using it atm)
## but user should know that this will not produce it again.

#make the CF1 reconstructed version too.
test=CF_16_daily[5].PSAL

[cf1date,cf1time,cf1prs,cf1mnprs,cf1sal,cf1tmp]=pd.read_pickle(open(datadir+'OSNAP2016recovery/pickles/CF1recon/CF1_recon_JH1810.pickle', 'rb'))
#make sal, ptmp, and pressure matrices

for cc in cf1date:
    print(len(cf1date[cc]))

cf1keys=sort(list(cf1date.keys()))
cf1_datevec=cf1date[cf1keys[0]]

cf1_prs_mat=zeros((len(cf1_datevec),len(cf1keys)))
cf1_sal_mat=cf1_prs_mat.copy()
cf1_ptmp_mat=cf1_prs_mat.copy()
for kk in range(len(cf1keys)):
    cf1_prs_mat[:,kk]=cf1prs[cf1keys[kk]]
    cf1_sal_mat[:,kk]=cf1sal[cf1keys[kk]]
    cf1_ptmp_mat[:,kk]=cf1tmp[cf1keys[kk]]

CF1_16_recon_daily=xr.Dataset({'PRES': (['TIME','DEPTH'], cf1_prs_mat),
            'PSAL': (['TIME','DEPTH'], cf1_sal_mat),
            'PTMP': (['TIME','DEPTH'],  cf1_ptmp_mat),},
            coords={'TIME': cf1_datevec,
                    'DEPTH': cf1keys,
                    'LATITUDE': CFlat[0],
                    'LONGITUDE': CFlon[0]})

add_SA_CT_PT(CF1_16_recon_daily)

CF1_16_recon_daily.to_netcdf(datadir+'OSNAP2016recovery/mcat_nc/CF1_recon_2016recovery_dailymerged.nc','w',format='netCDF4')
