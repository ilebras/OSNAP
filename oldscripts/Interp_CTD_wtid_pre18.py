#############################################################################
##############################################################################
#Interpolate CTD data vertically onto a regular grid
#############################################################################
#############################################################################

from aux_funcs import *

def Interp_CTD(moornum):


    savename='atom'


    #############################################################################
    #### Start by loading all temp and sal data
    #############################################################################

    moorname='CF'+str(moornum)

    tidlist=glob.glob(datadir+'TidbiT_Data_CF/OS_OSNAP-'+moorname+'*.nc')

    rbrlist=glob.glob(datadir+'RBR/'+moorname+'*solo*_ilebras.mat')

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
            tmp_hrly=gsw.pt0_from_t(sal_hrly,tmpins_hrly,prs_hrly)


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

    def pivspline_ML(afield):
        panpiv=pd.DataFrame(index=datevec,
                            columns=prsvec)
        for ii,dd in enumerate(datevec):
            ind=pdall['date bin']==dd
            x=pdall['pressure'][ind]
            y=pdall[afield][ind]
            Y=array([Y for X, Y in sorted(zip(x, y))])
            if sum(isnan(Y))!=len(Y):
                f1=interpolate.interp1d(hstack((0,sort(x)[~isnan(Y)])),
                                        hstack((Y[~isnan(Y)][0],Y[~isnan(Y)])),
                                        kind='linear',
                                        fill_value="extrapolate")
                pani=f1(prsvec)
                panpiv.iloc[ii]=pd.to_numeric(pani)
        return panpiv.T

    salinterp=pivspline_ML('salinity')
    tmpinterp=pivspline_ML('temperature')

    figure()
    plot(salinterp,salinterp.index);
    plot(pdall['salinity'],pdall['pressure'],'k.')
    gca().invert_yaxis()
    ylabel('pressure')
    title(moorname+': Salinity: measured and interpolated')
    savefig('../figures/interpolation/TS/'+moorname+'_sal_measinterp'+savename+'.png',bbox_inches='tight')


    figure()
    plot(tmpinterp,tmpinterp.index);
    plot(pdall['temperature'],pdall['pressure'],'k.')
    gca().invert_yaxis()
    ylabel('pressure')
    title(moorname+': Temperature (mcat only): measured and interpolated')
    savefig('../figures/interpolation/TS/'+moorname+'_tmp_mcat_measinterp'+savename+'.png',bbox_inches='tight')

    plotcontour(salinterp,cm.YlGnBu_r,30,35,moorname)
    title(moorname,fontsize=25)
    savefig('../figures/interpolation/TS/'+moorname+'_sal_dailyspline'+savename+'.png',bbox_inches='tight')


    plotcontour(tmpinterp,cm.RdBu_r,-2,8,moorname)
    savefig('../figures/interpolation/TS/'+moorname+'tmp_mcat_dailyspline'+savename+'.png',bbox_inches='tight')

    #############################################################################
    #################### Now add TidBit temperature ##############################
    #############################################################################


    #The depth tidbits (and rbr) are supposed to be at
    nomdpth={}
    nomdpth['CF1']=[60,80,140,120,70]#160,(was in first position)
    nomdpth['CF2']=[160,80,120,140,70,60]
    nomdpth['CF3']=[80,180,120,160,140]#70, (was third to last)
    nomdpth['CF4']=[75,150]
    nomdpth['CF5']=[150]
    nomdpth['CF6']=[200,150,300]
    nomdpth['CF7']=[200,300,150]


    lat=60
    nomprs={}
    for key in nomdpth:
        nomprs[key]=[gsw.p_from_z(-mm,lat) for mm in nomdpth[key]]


    time_tid=array([])
    datebin_tid=array([])
    tmp_tid=array([])
    meanprs_tid=array([])
    ii=0
    tidcon=2*24
    for dd in tidlist:
        print(dd)
        print(nomprs[moorname][ii])
        dat = Dataset(dd, 'r')
        time_tid_hrly=list(dat.variables['TIME'][:])[::tidcon][:-1]
        time_tid=hstack((time_tid,time_tid_hrly))
        datebin_tid_hrly=[0]*len(time_tid_hrly)
        tt=time_tid_hrly[0]
        for jj in arange(len(time_tid_hrly)):
            datebin_tid_hrly[jj]=datetime.datetime(1950,1,1)+datetime.timedelta(days=float(int(tt)+jj))
        datebin_tid_hrly=array(datebin_tid_hrly)
        tmp_tid_hrly=hrly_ave([float(tt) for tt in dat.variables['TEMP'][:].flatten()],tidcon)[:len(time_tid_hrly)]

        #remove data beyond first blowdown event for CF1 (not sure where they are and doesn't add much...)
        if moorname=='CF1':
            if nomprs[moorname][ii]<90:
                cf1blow1=datetime.datetime(2015, 5, 18, 19, 0)
                datebin_tid_tmp=datebin_tid_hrly[datebin_tid_hrly<cf1blow1]
                tmp_tid_tmp=tmp_tid_hrly[datebin_tid_hrly<cf1blow1]
            elif nomprs[moorname][ii]<130:
                cf1blow2=datetime.datetime(2015,8,15,13,0)
                datebin_tid_tmp=datebin_tid_hrly[datebin_tid_hrly<cf1blow2]
                tmp_tid_tmp=tmp_tid_hrly[datebin_tid_hrly<cf1blow2]
            if True==True:
                datebin_tid_hrly=datebin_tid_tmp.copy()
                tmp_tid_hrly=tmp_tid_tmp.copy()

        datebin_tid=hstack((datebin_tid,datebin_tid_hrly))
        tmp_tid=hstack((tmp_tid,tmp_tid_hrly))
        meanprs_tid=hstack((meanprs_tid,[nomprs[moorname][ii]]*ones(len(datebin_tid_hrly))))
        ii+=1

    #################################################################################
    ########    Add RBR if applicable!!        ######################################
    #################################################################################

    time_rbr=array([])
    tmp_rbr=array([])
    meanprs_rbr=array([])
    if (moorname=='CF4') | (moorname=='CF6') | (moorname=='CF7'):
        nomdpth_rbr={}
        nomdpth_rbr['CF4']=[300]
        nomdpth_rbr['CF6']=[400]
        nomdpth_rbr['CF7']=[400]
        nomprs_rbr={}
        nomprs_rbr[moorname]=[gsw.p_from_z(-mm,lat) for mm in nomdpth_rbr[moorname]]
        ii=0
        for dd in rbrlist:
            print(dd,nomprs_rbr[moorname][ii])
            dat = io.loadmat(dd)
            hourdiv=60*24
            time_rbr_hrly=dat['dtnum'].flatten()[::hourdiv]
            time_rbr=hstack((time_rbr,time_rbr_hrly))
            tmp_rbr_hrly=hrly_ave([float(tt) for tt in dat['temp'].flatten()],hourdiv)
            tmp_rbr=hstack((tmp_rbr,tmp_rbr_hrly))
            meanprs_rbr=hstack((meanprs_rbr,[nomprs_rbr[moorname][ii]]*ones(len(time_rbr_hrly))))
            ii+=1

    datebin_rbr=array([datetime.datetime(1,1,1)+datetime.timedelta(days=float(int(tt-366))) for tt in time_rbr])


    tidcorr=0.05

    tidpd=pd.DataFrame({'temperature':hstack((tmp_tid,tmp_rbr)),'temperature with corr':hstack((tmp_tid-tidcorr,tmp_rbr)),
                        'date bin':hstack((datebin_tid,datebin_rbr)),'nominal pressure':hstack((meanprs_tid,meanprs_rbr))})
    tidpd=tidpd[(tidpd['date bin']>=min(pdall['date bin'])) & (tidpd['date bin']<=max(pdall['date bin']))]


    plot(tidpd['date bin'])
    plot(pdall['date bin'])


    tid_wcorr_timebinned=tidpd.pivot_table(columns='nominal pressure',
                                           index='date bin',
                                           values='temperature with corr')
    tid_wcorr_melted=pd.melt(tid_wcorr_timebinned)
    if moorname=='CF8':
        tid_wcorr_melted.columns=['n.a','nominal pressure','temperature']
    else:
        tid_wcorr_melted.columns=['nominal pressure','temperature']

    ############################################################################################
    ##################### Reconstruct tidbit pressure ################################
    ############################################################################################


    mcat_pressures=pdall.pivot_table(index='nominal pressure',columns='date bin',values='pressure')

    axtprs=figure(figsize=(14,3))
    recon_prs=array([])
    for nn in unique(tid_wcorr_melted['nominal pressure']):

        timeind=(mcat_pressures.columns>=min(tid_wcorr_timebinned.index))&(mcat_pressures.columns<=max(tid_wcorr_timebinned.index))

        pseries=zeros(len(mcat_pressures.iloc[0][timeind]))

        ii=0
        for mm in mcat_pressures.index:
            if (nn>mm):
                # print(nn,mm)
                pseries=nn + (mcat_pressures.iloc[ii][timeind]-mm)*(mcat_pressures.index[ii+1]-nn)/(mcat_pressures.index[ii+1]-mm) + (mcat_pressures.iloc[ii+1][timeind]-mcat_pressures.index[ii+1])*(nn-mm)/(mcat_pressures.index[ii+1]-mm)
            ii+=1

        recon_prs=hstack((recon_prs,pseries))
        plot(pseries,label=str(nn))
    legend(loc=(1.05,0))
    plot(mcat_pressures.T,'k');
    savefig('../figures/interpolation/TS/'+moorname+'_pressure_recon_tidbit'+savename+'.png',bbox_inches='tight')

    tid_wcorr_melted['date bin']=tile(tid_wcorr_timebinned.index,len(tidlist)+len(rbrlist))

    tid_wcorr_melted['pressure']=recon_prs

    ################################################################################
    ##### combine mcat and tid data ####################################
    ################################################################################


    def pivspline_tmpadd(apanda):
        panpiv=pd.DataFrame(index=datevec,
                            columns=prsvec)
        alldates=append(pdall['date bin'],apanda['date bin'])
        allpressures=append(pdall['pressure'],apanda['pressure'])
        alltmps=append(pdall['temperature'],apanda['temperature'])
        for ii,dd in enumerate(sort(unique(alldates))):
            ind=alldates==dd
            x=allpressures[ind]
            y=alltmps[ind]
            Y=array([Y for X, Y in sorted(zip(x, y))])
            if sum(isnan(Y))!=len(Y):
                f1=interpolate.interp1d(hstack((0,sort(x)[~isnan(Y)])),
                                        hstack((Y[~isnan(Y)][0],Y[~isnan(Y)])),
                                        kind='linear',fill_value='extrapolate')
                pani=f1(prsvec)
                panpiv.iloc[ii]=pd.to_numeric(pani)
        return panpiv.T


    if moorname!='CF8':
        tmpinterp_wtid=pivspline_tmpadd(tid_wcorr_melted)
        figure()
        plot(tmpinterp_wtid,tmpinterp_wtid.index,'r');
        plot(tmpinterp,tmpinterp.index,'k');
        plot(pdall['temperature'],pdall['pressure'],'k.')
        plot(tid_wcorr_melted['temperature'],tid_wcorr_melted['pressure'],'r.')
        gca().invert_yaxis()
        ylabel('pressure')
        title('Temperature (red=with tidbits): measured and interpolated')
        savefig('../figures/interpolation/TS/'+moorname+'_tmp_tid_measinterp'+savename+'.png',bbox_inches='tight')
        figure()
        plotcontour(tmpinterp_wtid,cm.RdBu_r,-2,8,moorname)
        savefig('../figures/interpolation/TS/'+moorname+'_tmp_wtid_dailyspline'+savename+'.png',bbox_inches='tight')
    elif moorname=='CF8':
        tmpinterp_wtid=tmpinterp.copy()

    salinterp = salinterp.apply(pd.to_numeric, errors='coerce')
    tmpinterp = tmpinterp.apply(pd.to_numeric, errors='coerce')
    tmpinterp_wtid = tmpinterp_wtid.apply(pd.to_numeric, errors='coerce')

    pden=pd.DataFrame(index=salinterp.index,
                            columns=salinterp.columns)
    pden_wtid=pd.DataFrame(index=salinterp.index,
                                 columns=salinterp.columns)
    for jj in range(shape(pden)[1]):
            nanind=(~isnan(salinterp.values[:,jj]))
            sal_nonan=salinterp.values[:,jj][nanind]
            tmp_nonan=tmpinterp.values[:,jj][nanind]
            tmp_wtid_nonan=tmpinterp_wtid.values[:,jj][nanind]
            pden.values[nanind,jj]=pd.to_numeric(gsw.sigma0(gsw.SA_from_SP
            gsw.CT_from_pt(sal_nonan,tmp_nonan)))
            pden_wtid.values[nanind,jj]=pd.to_numeric(gsw.sigma0(sal_nonan,
            gsw.CT_from_pt(sal_nonan,tmp_wtid_nonan)))



    figure()
    plot(pden_wtid,pden_wtid.index,'r');
    plot(pden,pden.index,'k');
    gca().invert_yaxis()
    ylabel('pressure')
    title('Potential density (red=with tidbits)')
    savefig('../figures/interpolation/TS/'+moorname+'_pden_tid_interp'+savename+'.png',bbox_inches='tight')

     ########################################################################################
     #################################### Save fields   ####################################
     ########################################################################################

    pickle.dump([salinterp,tmpinterp,pden],
                open('../pickles/TSinterp/'+moorname+'_saltmpinterp_'+savename+'_notid.pickle','wb'),protocol=2)


    pickle.dump([salinterp,tmpinterp_wtid,pden_wtid],
                open('../pickles/TSinterp/'+moorname+'_saltmpinterp_'+savename+'.pickle','wb'),protocol=2)
