#############################################################################
##############################################################################
#Interpolate CTD data vertically onto a regular grid
#############################################################################
#############################################################################

from aux_funcs import *

def pivspline_extrap(afield,datevec,prsvec,pdall,moornum):
    panpiv=pd.DataFrame(index=datevec,
                        columns=prsvec)
    for ii,dd in enumerate(datevec):
        ind=pdall['date bin']==dd
        x=pdall['pressure'][ind]
        y=pdall[afield][ind]
        Y=array([Y for X, Y in sorted(zip(x, y))])
        if sum(isnan(Y))!=len(Y):
            #first, interpolate between measuements we do have so that I can then select
            #consistent pressures to interpolate from
            prsnonan=sort(x)[~isnan(Y)]
            f1=interpolate.interp1d(prsnonan,Y[~isnan(Y)],kind='linear')
            interprs=prsvec[(prsvec<=max(prsnonan))&(prsvec>=min(prsnonan))]
            pan1=f1(interprs)
            #pressures between which to calculate extrapolation gradient
            #trying top 20db that exist (this can be v. variable!)
            #CF5 is worst(top inst can go down to 250m)
            z1=int((min(prsnonan)+3)/2)*2
            z2=z1+20
            # print(z1,min(prsnonan))
            s1=pan1[interprs==z1]
            s2=pan1[interprs==z2]
            s0=s1-(s2-s1)*z1/(z2-z1)
            f2=interpolate.interp1d(hstack((0,interprs)),hstack((s0,pan1)),kind='linear',fill_value='extrapolate')
            pani=f2(prsvec)
            panpiv.iloc[ii]=pd.to_numeric(pani)
    return panpiv.T

def Interp_CTD(moornum):
    savename='1803extrap'


    #############################################################################
    #### Start by loading all temp and sal data
    #############################################################################

    moorname='CF'+str(moornum)

    def flattendic(adic):
        flatdic=[]
        for key in adic:
            flatdic=hstack((flatdic,adic[key]))
        return flatdic

    # CF8 needs to be loaded from original files, did not do any corrections to it.
    if moorname=='CF8':
        ctdlist=hstack((glob.glob(datadir+'NOC_M1/nocm1_02_2015/microcat/*.microcat'),
                        glob.glob(datadir+'NOC_M1/nocm1_01_2014/microcat/*.microcat')))
        time=array([])
        tmp=array([])
        sal=array([])
        prs=array([])
        meanprs=array([])
        datebin=array([])
        aveconst=4*24
        for dd in ctdlist:
            dat=pd.read_csv(dd,header=11,sep='  ')
            prs_hrly=hrly_ave(dat.iloc[:,6],aveconst/4)
            d=unique([datetime.datetime(int(dat.iloc[ii,0]),
                    int(dat.iloc[ii,1]),
                    int(dat.iloc[ii,2])) for ii in range(len(dat))])[:len(prs_hrly)]


            time_hrly=array([datetime.datetime.toordinal(adate) for adate in d])

            tmpins_hrly=hrly_ave(dat.iloc[:,4],aveconst/4)
            c=hrly_ave(dat.iloc[:,5],aveconst/4)
            sal_hrly=con2sal(c,tmpins_hrly,prs_hrly)
            meanprs_hrly=nanmean(prs_hrly)*ones(len(prs_hrly))
            SA_hrly=gsw.SA_from_SP(sal_hrly,prs_hrly,CFlon[moornum-1],CFlat[moornum-1])
            tmp_hrly=gsw.pt0_from_t(SA_hrly,tmpins_hrly,prs_hrly)


            datebin=hstack((datebin,d))
            time=hstack((time,time_hrly))
            prs=hstack((prs,prs_hrly))
            meanprs=hstack((meanprs,meanprs_hrly))
            tmp=hstack((tmp,tmp_hrly))
            sal=hstack((sal,sal_hrly))

    elif moorname=='CF1':#CF1- load a reconstructed version (will eventually be choosing from a couple different options)
        [cf1date,cf1time,cf1prs,cf1mnprs,cf1sal,cf1tmp]=pd.read_pickle(open('../pickles/CF1recon/CF1_recon.pickle', 'rb'))
        tmp=flattendic(cf1tmp)
        sal=flattendic(cf1sal)
        prs=flattendic(cf1prs)
        datebin=flattendic(cf1date)
        meanprs=flattendic(cf1mnprs)
    else:
        [date_all,month_all,prs_all,sal_all,tmp_all]=pickle.load(open('../pickles/TSdailydic/TS_daily_dic_wcorr.pickle','rb'))
        tmp=flattendic(tmp_all[int(moornum)])
        sal=flattendic(sal_all[int(moornum)])
        prs=flattendic(prs_all[int(moornum)])
        datebin=flattendic(date_all[int(moornum)])
        meanprs_dic={}
        for key in prs_all[int(moornum)]:
            meanprs_dic[key]=mean(prs_all[int(moornum)][key])*ones(len(prs_all[int(moornum)][key]))
        meanprs=flattendic(meanprs_dic)

    #############################################################################
    #### Create panda of all data, pivot and interpolate
    #############################################################################

    pdall=pd.DataFrame({'nominal pressure':meanprs,'temperature':tmp,'salinity':sal,
                        'pressure':prs,'date bin':datebin})


    dmin=min(pdall['date bin'])
    dmax=max(pdall['date bin'])
    dlen=int(divmod((dmax-dmin).total_seconds(),60*60*24)[0])+1
    datevec = array([dmin + datetime.timedelta(days=float(x)) for x in range(0, dlen)])

    prsvec=arange(0,int(max(prs))+1,2)


    salinterp=pivspline_extrap('salinity',datevec,prsvec,pdall,moornum)
    tmpinterp=pivspline_extrap('temperature',datevec,prsvec,pdall,moornum)

    figure()
    plot(salinterp,salinterp.index);
    plot(pdall['salinity'],pdall['pressure'],'k.')
    gca().invert_yaxis()
    ylabel('pressure')
    title(moorname+': Salinity: measured and interpolated')
    savefig('../figures/interpolation/TS/'+moorname+'_sal_measinterp_'+savename+'.png',bbox_inches='tight')


    figure()
    plot(tmpinterp,tmpinterp.index);
    plot(pdall['temperature'],pdall['pressure'],'k.')
    gca().invert_yaxis()
    ylabel('pressure')
    title(moorname+': Temperature (mcat only): measured and interpolated')
    savefig('../figures/interpolation/TS/'+moorname+'_tmp_mcat_measinterp_'+savename+'.png',bbox_inches='tight')

    plotcontour(salinterp,cm.YlGnBu_r,30,35,moorname)
    savefig('../figures/interpolation/TS/'+moorname+'_sal_'+savename+'.png',bbox_inches='tight')


    plotcontour(tmpinterp,cm.RdBu_r,-2,8,moorname)
    savefig('../figures/interpolation/TS/'+moorname+'tmp_mcat_'+savename+'.png',bbox_inches='tight')


    salinterp = salinterp.apply(pd.to_numeric, errors='coerce')
    tmpinterp = tmpinterp.apply(pd.to_numeric, errors='coerce')

    pden=pd.DataFrame(index=salinterp.index,
                            columns=salinterp.columns)

    for jj in range(shape(pden)[1]):
            nanind=(~isnan(salinterp.values[:,jj]))
            sal_nonan=salinterp.values[:,jj][nanind]
            tmp_nonan=tmpinterp.values[:,jj][nanind]
            # tmp_wtid_nonan=tmpinterp_wtid.values[:,jj][nanind]
            SA=gsw.SA_from_SP(sal_nonan,salinterp.index[nanind],CFlon[moornum-1],CFlat[moornum-1])
            pden.values[nanind,jj]=pd.to_numeric(gsw.sigma0(SA,gsw.CT_from_pt(SA,tmp_nonan)))


    figure()
    plot(pden,tmpinterp.index);
    gca().invert_yaxis()
    ylabel('pressure')
    title(moorname+': Potential density')
    savefig('../figures/interpolation/TS/'+moorname+'_pden_mcat_measinterp_'+savename+'.png',bbox_inches='tight')
     ########################################################################################
     #################################### Save fields   ####################################
     ########################################################################################

    pickle.dump([salinterp,tmpinterp,pden],
                open('../pickles/TSinterp/'+moorname+'_saltmpinterp_'+savename+'_notid.pickle','wb'),protocol=2)


ctdlist=hstack((glob.glob(datadir+'NOC_M1/nocm1_02_2015/microcat/*.microcat'),
                glob.glob(datadir+'NOC_M1/nocm1_01_2014/microcat/*.microcat')))

ctdlist

Interp_CTD(8)
