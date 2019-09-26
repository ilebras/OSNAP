#############################################################################
#############################################################################
################### Interpolate AQD and ADCP vertically ######################
#############################################################################
#############################################################################

from aux_funcs import *

datadir=datadir+'OSNAP2016recovery/'


#############################################################################
################### Start with AQD data, then add ADCP ############
#############################################################################
#
# def lpfilt(vec,hm):
#     #Apply a 40hr lowpass (averaging) filter
#     ll=int(hm/2)
#     filtvec=zeros(len(vec))
#     for ff in range(ll):
#         filtvec[ff]=nanmean(vec[:ff])
#     for ff in range(ll,len(vec)-ll):
#         filtvec[ff]=nanmean(vec[ff-ll:ff+ll])
#     for ff in range(len(vec)-ll,len(vec)):
#         filtvec[ff]=nanmean(vec[:ff])
#     return filtvec

import ttide as tt

def HTS(vec): #turn vector from half hour to hour time step by simple quick averaging
    if mod(len(vec),2)==1:
        vec2=(vec[::2][:-1]+vec[1::2])/2
    else:
        vec2=(vec[::2]+vec[1::2])/2
    return vec2


import calendar




dmin=datetime.datetime(2014,8,10,0)
dmax=datetime.datetime(2016,8,1,0)
dlen=int(divmod((dmax-dmin).total_seconds(),60*60)[0])+1
datevec=array([dmin + datetime.timedelta(hours=x) for x in range(0, dlen)])

def ttide_filt(utest,vtest):
    #Subtract TIDES
    if sum(~isnan(utest))<3:
        ufilt=utest
        vfilt=vtest
    else:
        tfit_u = tt.t_tide(utest+1j*vtest,dt=1);
        tidepart=tfit_u(arange(0,len(utest)/2,0.5));
        ufilt=utest-tidepart.real;
        vfilt=vtest-tidepart.imag;
    return ufilt,vfilt


def FilterExtract_AQD_ADCP(moornum):
    moorname='CF'+str(moornum)
    if moorname=='CF8':
        aqdlist=hstack((glob.glob(datadir+'NOC_M1/nocm1_01_2014/nor/*.edt'),
                        glob.glob(datadir+'NOC_M1/nocm1_02_2015/nor/*.edt')))
    elif moorname=='CF9':
        aqdlist=glob.glob(datadir+'NOC_M2/*Nortek.nc')
    else:
        aqdlist=glob.glob(datadir+'AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc')

    time=array([])
    u=array([])
    v=array([])
    dpth=array([])
    date=array([])

    for ii,dd in enumerate(aqdlist):
        print(dd)
        if moorname=='CF8':
            dat=pd.read_csv(dd,header=12,sep='\s*')
            date_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                               int(dat.ix[ii,1]),
                                               int(dat.ix[ii,2]),
                                               int(dat.ix[ii,3])) for ii in range(len(dat))])
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            ptest=array(dat.ix[:,5])
            ufilt,vfilt=ttide_filt(utest,vtest)

            u_tmp=np.interp(toTime(datevec),toTime(date_tmp),ufilt,left=NaN,right=NaN)
            v_tmp=np.interp(toTime(datevec),toTime(date_tmp),vfilt,left=NaN,right=NaN)
            dpth_tmp=np.interp(toTime(datevec),toTime(date_tmp),gsw.z_from_p(ptest,CFlat[moornum-1]),left=NaN,right=NaN)

        elif moorname=='CF9':
            #note: I should model all subsequent manipulations off this, xarray so much better than anyhting else...
            dat=xr.open_dataset(dd)
            utest=dat['UCUR'].values.flatten()
            vtest=dat['VCUR'].values.flatten()
            dtest=gsw.z_from_p(dat['PRES'].values.flatten(),CFlat[-1])
            ufilt,vfilt=ttide_filt(utest,vtest)
            date_tmp=np64ToDatetime(dat['TIME'])
            u_tmp=np.interp(toTime(datevec),toTime(date_tmp),ufilt,left=NaN,right=NaN)
            v_tmp=np.interp(toTime(datevec),toTime(date_tmp),vfilt,left=NaN,right=NaN)
            dpth_tmp=np.interp(toTime(datevec),toTime(date_tmp),dtest,left=NaN,right=NaN)


        else:
            dat = Dataset(dd, 'r')
            utest=array([float(tt) for tt in dat.variables['UCUR'][:].flatten()])
            vtest=array([float(tt) for tt in dat.variables['VCUR'][:].flatten()])
            ptest=array(dat.variables['PRES'][:].flatten())
            ufilt,vfilt=ttide_filt(HTS(utest),HTS(vtest))
            dtest=gsw.z_from_p(HTS(ptest),CFlat[moornum-1])
            time_tmp=list(dat.variables['TIME'][:])[::2][:len(dtest)]
            date_tmp=array([datetime.datetime(1950,1,1)+datetime.timedelta(int(tt*24)/24)
                            for tt in time_tmp])
            u_tmp=np.interp(toTime(datevec),toTime(date_tmp),ufilt,left=NaN,right=NaN)
            v_tmp=np.interp(toTime(datevec),toTime(date_tmp),vfilt,left=NaN,right=NaN)
            dpth_tmp=np.interp(toTime(datevec),toTime(date_tmp),dtest,left=NaN,right=NaN)


        if ii==0:
            dpth=dpth_tmp
            u=u_tmp
            v=v_tmp
        else:
            dpth=vstack((dpth,dpth_tmp))
            u=vstack((u,u_tmp))
            v=vstack((v,v_tmp))

    ##### Now add ADCP data

    ua=array([])
    va=array([])
    dptha=array([])
    meanda=array([])
    datea=array([])
    if moorname=='CF8':
        adcplist=hstack((glob.glob(datadir+'NOC_M1/nocm1_01_2014/adcp/*.edt'),
                        glob.glob(datadir+'NOC_M1/nocm1_02_2015/adcp/*.edt')))
        for iii,dd in enumerate(adcplist):
            dat=pd.read_csv(dd,header=12,sep='\s*')
            utest=array(dat.ix[:,6]/100)
            vtest=array(dat.ix[:,7]/100)
            dtest=-dat.ix[:,4]
            utest[utest<-90]=NaN
            vtest[vtest<-90]=NaN
            ufilt,vfilt=ttide_filt(utest,vtest)
            datea_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                                  int(dat.ix[ii,1]),
                                                  int(dat.ix[ii,2]),int(dat.ix[ii,3])) for ii in range(len(dat))])

            ua_tmp=np.interp(toTime(datevec),toTime(datea_tmp),ufilt,left=NaN,right=NaN)
            va_tmp=np.interp(toTime(datevec),toTime(datea_tmp),vfilt,left=NaN,right=NaN)
            dptha_tmp=np.interp(toTime(datevec),toTime(datea_tmp),dtest,left=NaN,right=NaN)

            if iii==0:
                dptha=dptha_tmp
                ua=ua_tmp
                va=va_tmp
            else:
                dptha=vstack((dptha,dptha_tmp))
                ua=vstack((ua,ua_tmp))
                va=vstack((va,va_tmp))

    elif moorname=='CF9':
        adcplist=glob.glob(datadir+'NOC_M2/*ADCP*.nc')
        for iii,aa in enumerate(adcplist):
            dat=xr.open_dataset(aa)
            datea_tmp=np64ToDatetime(dat['TIME'])
            bindim=len(dat.BINDEPTH)
            datedim=len(datevec)
            dptha_tmp=zeros((bindim,datedim))
            ua_tmp=zeros((bindim,datedim))
            va_tmp=zeros((bindim,datedim))
            for ii,bb in enumerate(dat.BINDEPTH):
                dtest=gsw.z_from_p(dat['PRES'].sel(BINDEPTH=bb).values.flatten(),CFlat[-1])
                utest=dat['UCUR'].sel(BINDEPTH=bb).values.flatten()
                vtest=dat['VCUR'].sel(BINDEPTH=bb).values.flatten()
                ufilt,vfilt=ttide_filt(utest,vtest)

                ua_tmp[ii,:]=np.interp(toTime(datevec),toTime(datea_tmp),ufilt,left=NaN,right=NaN)
                va_tmp[ii,:]=np.interp(toTime(datevec),toTime(datea_tmp),vfilt,left=NaN,right=NaN)
                dptha_tmp[ii,:]=np.interp(toTime(datevec),toTime(datea_tmp),dtest,left=NaN,right=NaN)

            if iii==0:
                dptha=dptha_tmp
                ua=ua_tmp
                va=va_tmp
            else:
                dptha=vstack((dptha,dptha_tmp))
                ua=vstack((ua,ua_tmp))
                va=vstack((va,va_tmp))


    else:
        dat=io.loadmat(datadir+'ADCP_Data_CF/OSNAP_cf'+str(moornum)+'_Final1_ilebras.mat')
        shapmat=shape(dat['u'][1:,:])
        shapmat1=len(datevec)
        # shapmat1=shape(dat['u'])[1]
        datea_tmp=[datetime.datetime(int(tt[0]),int(tt[1]),int(tt[2]),int(tt[3])) for tt in dat['time']]
    #     Note that adcp "depths" from Dan Torres are also in db...
        dptha_tmp=zeros((shapmat[0],shapmat1))
        ua_tmp=zeros((shapmat[0],shapmat1))
        va_tmp=zeros((shapmat[0],shapmat1))
        for ii in range(shapmat[0]):
            dtest=gsw.z_from_p(array(dat['z'][ii+1,:]),CFlat[moornum-1])
            utest=array(dat['u'][ii+1,:])
            vtest=array(dat['v'][ii+1,:])
            #ttide filter
            ufilt,vfilt=ttide_filt(utest,vtest)

            ua_tmp[ii,:]=np.interp(toTime(datevec),toTime(datea_tmp),ufilt,left=NaN,right=NaN)
            va_tmp[ii,:]=np.interp(toTime(datevec),toTime(datea_tmp),vfilt,left=NaN,right=NaN)
            dptha_tmp[ii,:]=np.interp(toTime(datevec),toTime(datea_tmp),dtest,left=NaN,right=NaN)

        dind=nanmean(dptha_tmp,axis=1)<0
        ua=ua_tmp[dind,:]
        va=va_tmp[dind,:]
        dptha=dptha_tmp[dind,:]

    utot=vstack((ua,u))
    vtot=vstack((va,v))
    dtot=vstack((dptha,dpth))

    return utot,vtot,dtot



# utot={}
# vtot={}
# dtot={}

ii=9
utot[ii],vtot[ii],dtot[ii]=FilterExtract_AQD_ADCP(ii)

plot(dtot[ii].T);


for ii in dtot:
    print(shape(dtot[ii]))


savedat=io.loadmat(open('/home/isabela/Documents/projects/OSNAP/data/OSNAP2016recovery/pickles/VELinterp/hourly/CF_M_hourlyvel_ttide_sorted.mat','rb'))

# savedat={}
# savedat['D']=NaN*ones((41,9,len(datevec)))
# savedat['U']=NaN*ones((41,9,len(datevec)))
# savedat['V']=NaN*ones((41,9,len(datevec)))
dtot

for ii in [8,9]:
    savedat['D'][:shape(dtot[ii])[0],ii-1,:]=dtot[ii]
    savedat['U'][:shape(dtot[ii])[0],ii-1,:]=utot[ii]*cos(theta)+vtot[ii]*sin(theta)
    savedat['V'][:shape(dtot[ii])[0],ii-1,:]=-utot[ii]*sin(theta)+vtot[ii]*cos(theta)
    figure()
    plot(savedat['U'][:,ii-1,:])


# io.savemat(open('/home/isabela/Documents/projects/OSNAP/data/OSNAP2016recovery/pickles/VELinterp/Astrid/CF_M_hourlyvel_40hrlpass.mat','wb'),savedat)
# The below is basically a seperate script, just sorting by depth

# from pylab import *
# from scipy import io
# import datetime
#
#
# dmin=datetime.datetime(2014,8,10,0)
# dmax=datetime.datetime(2016,8,1,0)
# dlen=int(divmod((dmax-dmin).total_seconds(),60*60)[0])+1
# datevec=array([dmin + datetime.timedelta(hours=x) for x in range(0, dlen)])
#
# savedat=io.loadmat(open('/home/isabela/Documents/projects/OSNAP/data/OSNAP2016recovery/pickles/VELinterp/Astrid/CF_M_hourlyvel_40hrlpass.mat','rb'))
#
# import xarray as xr
#

dsets={}
for ii in range(9):
    dsets[ii]=xr.Dataset({'U': (['dnum','date'], savedat['U'][:,ii,:]),
                         'V': (['dnum','date'], savedat['V'][:,ii,:]),
                         'D': (['dnum','date'], savedat['D'][:,ii,:])},
                         coords={'date': datevec, 'dnum': nanmean(savedat['D'][:,ii,:],axis=-1)})
dsets_sorted={}
for ii in range(9):
    dsets_sorted[ii]=dsets[ii].sortby(-1*dsets[ii].dnum)

dsets_sorted

for ii in range(9):
    figure()
    plot(dsets_sorted[ii].date,dsets_sorted[ii].D.T)

savedat_sorted={}
savedat_sorted['D']=NaN*ones((41,9,len(datevec)))
savedat_sorted['U']=NaN*ones((41,9,len(datevec)))
savedat_sorted['V']=NaN*ones((41,9,len(datevec)))

for ii in range(9):
    for key in savedat_sorted:
        savedat_sorted[key][:,ii,:]=dsets_sorted[ii][key]

io.savemat(open('/home/isabela/Documents/projects/OSNAP/data/OSNAP2016recovery/pickles/VELinterp/hourly/CF_M_2014-2016_hourlyvel_ttide_sorted_rev1_1903.mat','wb'),savedat_sorted)
