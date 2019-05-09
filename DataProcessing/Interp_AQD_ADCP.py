#############################################################################
#############################################################################
################### Interpolate AQD and ADCP vertically ######################
#############################################################################
#############################################################################

from aux_funcs import *

#############################################################################
################### Start with AQD data, then add ADCP and interp############
#############################################################################

import ttide as tt

moorname='CF1'
aqdlist=glob.glob(datadir+'OSNAP2016recovery/AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc')
dat = Dataset(aqdlist[0], 'r')
utest=array([float(tt) for tt in dat.variables['UCUR'][:].flatten()])
vtest=array([float(tt) for tt in dat.variables['VCUR'][:].flatten()])
tfit_u = tt.t_tide(utest+1j*vtest,dt=0.5)
tidepart=tfit_u(arange(0,len(utest)/2,0.5))
ufilt=utest-tidepart.real
vfilt=vtest-tidepart.imag

figure(figsize=(15,5))
plot(utest)
plot(tidepart.real)
figure(figsize=(15,5))
plot(vtest)
plot(tidepart.imag)
mean(vtest)
mean(vfilt)
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
            prs_tmp=hrly_ave(dat.ix[:,5],hrcon)
            date_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                               int(dat.ix[ii,1]),
                                            int(dat.ix[ii,2])) for ii in range(len(dat))])[:len(prs_tmp)]
            time_tmp=array([datetime.datetime.toordinal(adate) for adate in date_tmp])
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            tfit_u = tt.t_tide(utest)
            tfit_v = tt.t_tide(vtest)
            ufilt= utest-tfit_u(arange(len(utest)))
            vfilt= vtest-tfit_v(arange(len(vtest)))
            u_tmp=hrly_ave(ufilt,hrcon)
            v_tmp=hrly_ave(vfilt,hrcon)
            mprs_tmp=mean(dat.ix[:,5])*ones(len(date_tmp))

        else:
            dat = Dataset(dd, 'r')
            utest=array([float(tt) for tt in dat.variables['UCUR'][:].flatten()])
            vtest=array([float(tt) for tt in dat.variables['VCUR'][:].flatten()])
            tfit_u = tt.t_tide(utest+1j*vtest,dt=0.5)
            tidepart=tfit_u(arange(0,len(utest)/2,0.5))
            ufilt=utest-tidepart.real
            vfilt=vtest-tidepart.imag
            u_tmp=hrly_ave(ufilt,hrcon*2)
            v_tmp=hrly_ave(vfilt,hrcon*2)
            prs_tmp=hrly_ave(list(dat.variables['PRES'][:].flatten()),hrcon*2)
            time_tmp=list(dat.variables['TIME'][:])[::hrcon*2][:len(prs_tmp)]
            mprs_tmp=nanmean(prs_tmp)*ones(len(time_tmp))
            date_tmp=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=int(tt))
                            for tt in time_tmp])
            u_tmp=u_tmp[:len(date_tmp)]
            v_tmp=v_tmp[:len(date_tmp)]
            prs_tmp=prs_tmp[:len(date_tmp)]

            figure(figsize=(14,6))
            subplot(211)
            title(str(int(nanmean(prs_tmp))))
            plot(utest)
            plot(tidepart.real)
            subplot(212)
            plot(vtest)
            plot(tidepart.imag)


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
            prs_tmp=hrly_ave(gsw.p_from_z(-dat.ix[:,4],lat),hrcon)
            datea_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                      int(dat.ix[ii,1]),
                                      int(dat.ix[ii,2])) for ii in range(len(dat))])[:len(prs_tmp)]
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            utest[utest<-90]=NaN
            vtest[vtest<-90]=NaN
            tfit_u = tt.t_tide(utest)
            tfit_v = tt.t_tide(vtest)
            ufilt= utest-tfit_u(arange(len(utest)))
            vfilt= vtest-tfit_v(arange(len(vtest)))
            ua_tmp=hrly_ave(ufilt,hrcon)
            va_tmp=hrly_ave(vfilt,hrcon)
            prsa=hstack((prsa,prs_tmp))
            ua=hstack((ua,ua_tmp))
            va=hstack((va,va_tmp))
            datea=hstack((datea,datea_tmp))
            plot(datea_tmp,ua_tmp)
            plot(datea_tmp,va_tmp)

    else:
        dat=io.loadmat(datadir+'ADCP_Data_CF/OSNAP_cf'+str(moornum)+'_Final1_ilebras.mat')
        shapmat=shape(dat['u'][1:,:])
        shapmat1=int(shape(dat['u'])[1]/24)+1
        datea_tmp=unique([datetime.datetime(int(tt[0]),int(tt[1]),int(tt[2])) for tt in dat['time']])[:shapmat1]
    #     Note that adcp "depths" from Dan Torres are also in db...
        prsa_tmp=zeros((shapmat[0],shapmat1))
        ua_tmp=zeros((shapmat[0],shapmat1))
        va_tmp=zeros((shapmat[0],shapmat1))
        for ii in range(shapmat[0]):
            prsa_tmp[ii,:]=hrly_ave(dat['z'][ii+1,:],hrcon)
            utest=array(dat['u'][ii+1,:])
            vtest=array(dat['v'][ii+1,:])
            # if sum(~isnan(utest))>len(utest)/2:
            #     # print(utest)
            #     tfit_u = tt.t_tide(utest+1j*vtest)
            #     tidepart=tfit_u(arange(len(utest)))
            #     ufilt=utest-tidepart.real
            #     vfilt=vtest-tidepart.imag
            #     figure(figsize=(14,6))
            #     subplot(211)
            #     title(str(int(nanmean(prsa_tmp[ii,:]))))
            #     plot(utest)
            #     plot(tidepart.real)
            #     subplot(212)
            #     plot(vtest)
            #     plot(tidepart.imag)
            #     ua_tmp[ii,:]=hrly_ave(utest,hrcon)
            #     va_tmp[ii,:]=hrly_ave(vtest,hrcon)
            #     title(str(int(nanmean(vtest)*100)/100)+'/'+str(int(nanmean(va_tmp[ii,:])*100)/100))
            # else:
            ua_tmp[ii,:]=hrly_ave(utest,hrcon)
            va_tmp[ii,:]=hrly_ave(vtest,hrcon)

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

    versname='1808tidecomp'

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

Interp_AQD_ADCP(1)
