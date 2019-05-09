#############################################################################
#############################################################################
################### Interpolate AQD and ADCP vertically ######################
#############################################################################
#############################################################################

from aux_funcs import *

#############################################################################
################### Start with AQD data, then add ADCP and interp############
#############################################################################

# def extrapsurf(afield,datevec,prsvec,pdall):
#     panpiv=pd.DataFrame(index=datevec,
#                         columns=prsvec)
#     for ii,dd in enumerate(datevec):
#         ind=pdall['date bin']==dd
#         x=pdall['pressure'][ind]
#         y=pdall[afield][ind]
#         Y=array([Y for X, Y in sorted(zip(x, y))])
#         if sum(isnan(Y))!=len(Y):
#                 #first, interpolate between measuements we do have so that I can then select
#                 #consistent pressures to interpolate from
#                 prsnonan=sort(x)[~isnan(Y)]
#                 f1=interpolate.interp1d(prsnonan,Y[~isnan(Y)],kind='linear')
#                 interprs=prsvec[(prsvec<=max(prsnonan))&(prsvec>=min(prsnonan))]
#                 pan1=f1(interprs)
#                 #pressures between which to calculate extrapolation gradient
#                 #trying top 20db that exist (this can be v. variable!)
#                 #CF5 is worst(top inst can go down to 250m)
#                 z1=int((min(prsnonan)+3)/2)*2
#                 z2=z1+20
#                 # print(z1,min(prsnonan))
#                 s1=pan1[interprs==z1]
#                 s2=pan1[interprs==z2]
#                 s0=s1-(s2-s1)*z1/(z2-z1)
#                 f2=interpolate.interp1d(hstack((0,interprs)),hstack((s0,pan1)),kind='linear',fill_value='extrapolate')
#                 pani=f2(prsvec)
#                 panpiv.iloc[ii]=pd.to_numeric(pani)
#     return panpiv.T


def fillsurf(afield,datevec,prsvec,pdall):
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
                f1=interpolate.interp1d(hstack((0,prsnonan)),hstack((Y[~isnan(Y)][0],Y[~isnan(Y)])),kind='linear',fill_value='extrapolate')
                pan1=f1(prsvec)
                panpiv.iloc[ii]=pd.to_numeric(pan1)
    return panpiv.T

def lpfilt(vec,hm):
    #Apply a 40hr lowpass (averaging) filter
    ll=int(hm/2)
    filtvec=zeros(len(vec))
    for ff in range(ll):
        filtvec[ff]=nanmean(vec[:ff])
    for ff in range(ll,len(vec)-ll):
        filtvec[ff]=nanmean(vec[ff-ll:ff+ll])
    for ff in range(len(vec)-ll,len(vec)):
        filtvec[ff]=nanmean(vec[:ff])
    return filtvec



def Interp_AQD_ADCP(moornum):

    moorname='CF'+str(moornum)


    if moorname=='CF8':
        aqdlist=hstack((glob.glob(datadir+'NOC_M1/nocm1_01_2014/nor/*.edt'),
                        glob.glob(datadir+'NOC_M1/nocm1_02_2015/nor/*.edt')))
    else:
        aqdlist=glob.glob(datadir+'AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc')



    figure(figsize=(12,10))

    time=array([])
    u=array([])
    v=array([])
    prs=array([])
    meanprs=array([])
    date=array([])

    for dd in aqdlist:
        print(dd)
        if moorname=='CF8':
            dat=pd.read_csv(dd,header=12,sep='\s*')
            date_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                               int(dat.ix[ii,1]),
                                            int(dat.ix[ii,2])) for ii in range(len(dat))])[:len(prs_tmp)]
            time_tmp=array([datetime.datetime.toordinal(adate) for adate in date_tmp])
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            ptest=array(dat.ix[:,5])
            ufilt=lpfilt(utest)
            vfilt=lpfilt(vtest)
            pfilt=lpfilt(ptest)
            u_tmp=ufilt[:24]
            v_tmp=vfilt[:24]
            prs_tmp=pfilt[::24]
            mprs_tmp=mean(dat.ix[:,5])*ones(len(date_tmp))

        else:
            dat = Dataset(dd, 'r')
            utest=array([float(tt) for tt in dat.variables['UCUR'][:].flatten()])
            vtest=array([float(tt) for tt in dat.variables['VCUR'][:].flatten()])
            ptest=list(dat.variables['PRES'][:].flatten())
            ufilt=lpfilt(utest,80)
            vfilt=lpfilt(vtest,80)
            pfilt=lpfilt(ptest,80)
            u_tmp=ufilt[::48]
            v_tmp=vfilt[::48]
            prs_tmp=pfilt[::48]
            time_tmp=list(dat.variables['TIME'][:])[::hrcon*2][:len(prs_tmp)]
            mprs_tmp=nanmean(prs_tmp)*ones(len(time_tmp))
            date_tmp=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=int(tt))
                            for tt in time_tmp])
            u_tmp=u_tmp[:len(date_tmp)]
            v_tmp=v_tmp[:len(date_tmp)]
            prs_tmp=prs_tmp[:len(date_tmp)]

        prs=hstack((prs,prs_tmp))
        u=hstack((u,u_tmp))
        v=hstack((v,v_tmp))
        meanprs=hstack((meanprs,mprs_tmp))
        date=hstack((date,date_tmp))
        time=hstack((time,time_tmp))

    ##### Choose time and pressure bins
    #units are db
    pstep=2

    ##### Now add ADCP data

    lat=60


    ua=array([])
    va=array([])
    prsa=array([])
    meanprs=array([])
    datea=array([])
    if moorname=='CF8':
        adcplist=hstack((glob.glob(datadir+'NOC_M1/nocm1_01_2014/adcp/*.edt'),
                        glob.glob(datadir+'NOC_M1/nocm1_02_2015/adcp/*.edt')))
        for dd in adcplist:
            dat=pd.read_csv(dd,header=12,sep='\s*')
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            ptest=gsw.p_from_z(-dat.ix[:,4],lat)
            utest[utest<-90]=NaN
            vtest[vtest<-90]=NaN
            ufilt= lpfilt(utest,40)
            vfilt= lpfilt(vtest,40)
            pfilt= lpfilt(ptest,40)
            ua_tmp=ufilt[::24]
            va_tmp=vfilt[::24]
            prs_tmp=pfilt[::24]

            datea_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                                  int(dat.ix[ii,1]),
                                                  int(dat.ix[ii,2])) for ii in range(len(dat))])[:len(prs_tmp)]
            prsa=hstack((prsa,prs_tmp))
            ua=hstack((ua,ua_tmp))
            va=hstack((va,va_tmp))
            datea=hstack((datea,datea_tmp))
            plot(datea_tmp,ua_tmp)
            plot(datea_tmp,va_tmp)

    else:
        dat=io.loadmat(datadir+'ADCP_Data_CF/OSNAP_cf'+str(moornum)+'_Final1_ilebras.mat')
        shapmat=shape(dat['u'][1:,:])
        shapmat1=int(shape(dat['u'])[1]/24)
        datea_tmp=unique([datetime.datetime(int(tt[0]),int(tt[1]),int(tt[2])) for tt in dat['time']])[:shapmat1]
    #     Note that adcp "depths" from Dan Torres are also in db...
        prsa_tmp=zeros((shapmat[0],shapmat1))
        ua_tmp=zeros((shapmat[0],shapmat1))
        va_tmp=zeros((shapmat[0],shapmat1))
        for ii in range(shapmat[0]):
            ptest=array(dat['z'][ii+1,:])
            utest=array(dat['u'][ii+1,:])
            vtest=array(dat['v'][ii+1,:])
            #40 hr filter
            pfilt=lpfilt(ptest,40)
            ufilt=lpfilt(utest,40)
            vfilt=lpfilt(vtest,40)
            #decimate to daily
            prsa_tmp[ii,:]=pfilt[::24][:shapmat1]
            ua_tmp[ii,:]=ufilt[::24][:shapmat1]
            va_tmp[ii,:]=vfilt[::24][:shapmat1]
            print('orig mean: ',nanmean(vtest))
            print('filtered mean: ',nanmean(vfilt))
            print('final mean: ',nanmean(va_tmp[ii,:]))

        ua=ua_tmp.flatten()
        va=va_tmp.flatten()
        prsa=prsa_tmp.flatten()

        datea=tile(datea_tmp,[1,shapmat[0]]).flatten()

    pdall=pd.DataFrame({'pressure':hstack((prsa,prs)),
                                   'u':hstack((ua,u)),
                                   'v':hstack((va,v)),
                                   'date bin':hstack((datea,date))})


    pdall['u'][pdall['u']<-2]=NaN
    pdall['v'][pdall['v']<-2]=NaN

    if moornum==4:
        plim=40
    else:
        plim=0

    pdall=pdall[pdall['pressure']>=plim]

    common_mindate=max(min(datea),min(date))
    common_maxdate=min(max(datea),max(date))

    # datevec=sort(unique(pdall['date bin']))
    prsvec=arange(0,int(max(pdall['pressure']))+1,2)


    dmin=min(pdall['date bin'])
    dmax=max(pdall['date bin'])
    dlen=int(divmod((dmax-dmin).total_seconds(),60*60*24)[0])+1
    datevec = array([dmin + datetime.timedelta(days=x) for x in range(0, dlen)])


    u_interp=fillsurf('u',datevec,prsvec,pdall)
    v_interp=fillsurf('v',datevec,prsvec,pdall)


    u_interp[(u_interp<-2.5) | (u_interp>2)]=nan
    v_interp[(v_interp<-2.5) | (v_interp>2)]=nan

    versname='1808lpfilt'

    figure(figsize=(12,4))
    subplot(121)
    plot(u_interp,prsvec);
    plot(pdall['u'],pdall['pressure'],'k.');
    gca().invert_yaxis()
    title('u, zonal velocity (m/s)')
    ylabel('pressure (db)')
    subplot(122)
    plot(v_interp,prsvec);
    plot(pdall['v'],pdall['pressure'],'k.');
    gca().invert_yaxis()
    title('v, meridional velocity (m/s)')
    savefig('../figures/interpolation/VEL/'+moorname+'_uv_measinterp_'+versname+'.png',bbox_inches='tight')

    plotcontour(u_interp,cm.RdBu_r,-2,2,moorname)
    savefig('../figures/interpolation/VEL/'+moorname+'_ucontour_'+versname+'.png',bbox_inches='tight')


    plotcontour(v_interp,cm.RdBu_r,-2,2,moorname)
    savefig('../figures/interpolation/VEL/'+moorname+'_vcontour_'+versname+'.png',bbox_inches='tight')



    pickle.dump([u_interp,v_interp],open('../pickles/VELinterp/'+moorname+'_uvinterp_'+versname+'.pickle','wb'))

Interp_AQD_ADCP(8)
