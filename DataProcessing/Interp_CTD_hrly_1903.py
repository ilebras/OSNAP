#############################################################################
##############################################################################
#Interpolate CTD data vertically onto a regular grid
#############################################################################
#############################################################################

from aux_funcs import *

datadir=datadir+'OSNAP2016recovery/'
dmin=datetime.datetime(2014,8,10,0)
dmax=datetime.datetime(2016,8,1,0)
dlen=int(divmod((dmax-dmin).total_seconds(),60*60)[0])+1
datevec=array([dmin + datetime.timedelta(hours=x) for x in range(0, dlen)])

## Start with CF1-CF7
[date_all,month_all,prs_all,sal_all,ptmp_all,pden_all]=pickle.load(open(datadir+'pickles/TSdailydic/TS_15min_dic_wJHdipcal.pickle','rb'))

#get rid of bad instrument at CF7

del date_all[7][500]
del month_all[7][500]
del prs_all[7][500]
del sal_all[7][500]
del ptmp_all[7][500]
del pden_all[7][500]



import calendar

xrint={}
for ii in sal_all:
    dpvec=array([])
    for jj,key in enumerate(sal_all[ii]):
        dpvec=hstack((key,dpvec))
        saltmp=np.interp(toTime(datevec),toTime(date_all[ii][key]),run_ave(sal_all[ii][key],4),left=NaN,right=NaN)
        tmptmp=np.interp(toTime(datevec),toTime(date_all[ii][key]),run_ave(ptmp_all[ii][key],4),left=NaN,right=NaN)
        dentmp=np.interp(toTime(datevec),toTime(date_all[ii][key]),run_ave(pden_all[ii][key],4),left=NaN,right=NaN)
        dphtmp=np.interp(toTime(datevec),toTime(date_all[ii][key]),gsw.z_from_p(run_ave(prs_all[ii][key],4),CFlat[ii-1]),left=NaN,right=NaN)
        if jj==0:
            salaccum=saltmp
            tmpaccum=tmptmp
            dphaccum=dphtmp
            denaccum=dentmp
        else:
            salaccum=vstack((saltmp,salaccum))
            tmpaccum=vstack((tmptmp,tmpaccum))
            dphaccum=vstack((dphtmp,dphaccum))
            denaccum=vstack((dentmp,denaccum))

    xrint[ii]=xr.Dataset({'psal': (['dpvec','date'],salaccum),
                          'ptmp': (['dpvec','date'],tmpaccum),
                          'sigma0': (['dpvec','date'],denaccum),
                          'depth': (['dpvec','date'],dphaccum)},
                          coords={'dpvec':dpvec,'date':datevec})
## Add CF8 and CF9
datint={}
for ii in [8,9]:
    print(ii)
    if ii==8:
        mcatlist=hstack((glob.glob(datadir+'NOC_M1/nocm1_01_2014/microcat/*.nc'),glob.glob(datadir+'NOC_M1/nocm1_02_2015/microcat/*.nc')))
    elif ii==9:
        mcatlist=glob.glob(datadir+'NOC_M2/*MCAT.nc')


    for jj,dd in enumerate(mcatlist):
        print(jj,dd)
        dat=xr.open_dataset(dd)
        datint[jj]=dat.interp(TIME=datevec)
        datint[jj]['DEPTH']=[int(pp)*100 for pp in datint[0].DEPTH.values/100]
        datint[jj]['LONGITUDE']=round(float(datint[jj]['LONGITUDE']),1)
        datint[jj]['LATITUDE']=round(float(datint[jj]['LATITUDE']),1)

    dat_tot=xr.merge([datint[0],datint[1]])

    SAvec=gsw.SA_from_SP(dat_tot['PSAL'],dat_tot['PRES'],dat_tot['LONGITUDE'],dat_tot['LONGITUDE'])

    PTMP=gsw.pt0_from_t(SAvec,dat_tot['TEMP'],dat_tot['PRES'])
    PDEN=gsw.sigma0(SAvec,gsw.CT_from_pt(SAvec,PTMP))
    DEPTHVAR=gsw.z_from_p(dat_tot['PRES'],CFlat[4])

    xrint[ii]=xr.Dataset({'psal': (['dpvec','date'],dat_tot['PSAL'].values.T),
                         'ptmp': (['dpvec','date'],PTMP.T),
                         'sigma0': (['dpvec','date'],PDEN.T),
                         'depth': (['dpvec','date'],DEPTHVAR.T)},
                         coords={'dpvec':dat_tot.DEPTH.values,'date':datevec})


xrint_sorted={}
for ii in range(1,10):
    print(len(xrint[ii].dpvec))
    xrint_sorted[ii]=xrint[ii].sortby(-1*xrint[ii].dpvec)


pickle.dump(xrint_sorted,open(datadir+'/pickles/xarray/CF_M_2014-2016_hourlyTSD_1904_noCF7500.pickle','wb'))


finaldatadic={}
finaldatadic['psal']=NaN*ones((10,9,len(datevec)))
finaldatadic['ptmp']=NaN*ones((10,9,len(datevec)))
finaldatadic['sigma0']=NaN*ones((10,9,len(datevec)))
finaldatadic['depth']=NaN*ones((10,9,len(datevec)))

for ii in range(1,10):
    for key in finaldatadic:
        dplen=len(xrint[ii].dpvec)
        finaldatadic[key][:dplen,ii-1,:]=xrint_sorted[ii][key].values

for key in finaldatadic:
    for ii in range(9):
        figure()
        plot(datevec,finaldatadic[key][:,ii,:].T)
len(datevec)


io.savemat(open('/home/isabela/Documents/projects/OSNAP/data/OSNAP2016recovery/pickles/TSinterp/hourly/CF_M_2014-2016_hourlyTSD_1904_noCF7500.mat','wb'),finaldatadic)



XXXXXXXXXXXX


#     #############################################################################
#     #### Create panda of all data, pivot and interpolate
#     #############################################################################
#
#     pdall=pd.DataFrame({'nominal pressure':meanprs,'temperature':tmp,'salinity':sal,
#                         'pressure':prs,'date bin':datebin})
#
#
#     dmin=min(pdall['date bin'])
#     dmax=max(pdall['date bin'])
#     dlen=int(divmod((dmax-dmin).total_seconds(),60*60*24)[0])+1
#     datevec = array([dmin + datetime.timedelta(days=float(x)) for x in range(0, dlen)])
#
#     prsvec=arange(0,int(max(prs))+1,2)
#
#
#     salinterp=pivspline_extrap('salinity',datevec,prsvec,pdall,moornum)
#     tmpinterp=pivspline_extrap('temperature',datevec,prsvec,pdall,moornum)
#
#     figure()
#     plot(salinterp,salinterp.index);
#     plot(pdall['salinity'],pdall['pressure'],'k.')
#     gca().invert_yaxis()
#     ylabel('pressure')
#     title(moorname+': Salinity: measured and interpolated')
#     savefig('../figures/interpolation/TS/'+moorname+'_sal_measinterp_'+savename+'.png',bbox_inches='tight')
#
#
#     figure()
#     plot(tmpinterp,tmpinterp.index);
#     plot(pdall['temperature'],pdall['pressure'],'k.')
#     gca().invert_yaxis()
#     ylabel('pressure')
#     title(moorname+': Temperature (mcat only): measured and interpolated')
#     savefig('../figures/interpolation/TS/'+moorname+'_tmp_mcat_measinterp_'+savename+'.png',bbox_inches='tight')
#
#     plotcontour(salinterp,cm.YlGnBu_r,30,35,moorname)
#     savefig('../figures/interpolation/TS/'+moorname+'_sal_'+savename+'.png',bbox_inches='tight')
#
#
#     plotcontour(tmpinterp,cm.RdBu_r,-2,8,moorname)
#     savefig('../figures/interpolation/TS/'+moorname+'tmp_mcat_'+savename+'.png',bbox_inches='tight')
#
#
#     salinterp = salinterp.apply(pd.to_numeric, errors='coerce')
#     tmpinterp = tmpinterp.apply(pd.to_numeric, errors='coerce')
#
#     pden=pd.DataFrame(index=salinterp.index,
#                             columns=salinterp.columns)
#
#     for jj in range(shape(pden)[1]):
#             nanind=(~isnan(salinterp.values[:,jj]))
#             sal_nonan=salinterp.values[:,jj][nanind]
#             tmp_nonan=tmpinterp.values[:,jj][nanind]
#             # tmp_wtid_nonan=tmpinterp_wtid.values[:,jj][nanind]
#             SA=gsw.SA_from_SP(sal_nonan,salinterp.index[nanind],CFlon[moornum-1],CFlat[moornum-1])
#             pden.values[nanind,jj]=pd.to_numeric(gsw.sigma0(SA,gsw.CT_from_pt(SA,tmp_nonan)))
#
#
#     figure()
#     plot(pden,tmpinterp.index);
#     gca().invert_yaxis()
#     ylabel('pressure')
#     title(moorname+': Potential density')
#     savefig('../figures/interpolation/TS/'+moorname+'_pden_mcat_measinterp_'+savename+'.png',bbox_inches='tight')
#      ########################################################################################
#      #################################### Save fields   ####################################
#      ########################################################################################
#
#     pickle.dump([salinterp,tmpinterp,pden],
#                 open('../pickles/TSinterp/'+moorname+'_saltmpinterp_'+savename+'_notid.pickle','wb'),protocol=2)
#
# for ii in range(1,8):
#     Interp_CTD(ii)
