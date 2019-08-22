from aux_funcs import *

def AQDlist(moornum):
    moorname='CF'+str(moornum)
    if moorname=='CF8':
        aqdlist=hstack((glob.glob(datadir+'OSNAP2016recovery/NOC_M1/nocm1_01_2014/nor/*.edt'),
                        glob.glob(datadir+'OSNAP2016recovery/NOC_M1/nocm1_02_2015/nor/*.edt')))
    else:
        aqdlist=glob.glob(datadir+'OSNAP2016recovery/AQD_Data_CF/OS_OSNAP-'+moorname+'*.nc')

    return sort(aqdlist)

AQD={}
AQD['CF5_500m']=AQDlist(5)[-2]
AQD['CF6_500m']=AQDlist(6)[-2]
AQD['M1_500m']=hstack((AQDlist(8)[2],AQDlist(8)[5]))


# Load the AQD files I want to look at
def loadAQD(moorname):

    dd=AQD[moorname+'_500m']
    if moorname=='M1':
        date=array([])
        u=array([])
        v=array([])
        p=array([])
        for d1 in dd:
            dat=pd.read_csv(d1,header=12,sep='\s*')
            date_tmp=[datetime.datetime(int(dat.ix[ii,0]),
                                            int(dat.ix[ii,1]),
                                            int(dat.ix[ii,2]),
                                            int(dat.ix[ii,3])) for ii in range(len(dat))]
            u_tmp=array(dat.ix[:,6]/100)
            v_tmp=array(dat.ix[:,7]/100)
            p_tmp=array(dat.ix[:,5])

            p=hstack((p,p_tmp))
            u=hstack((u,u_tmp))
            v=hstack((v,v_tmp))
            date=hstack((date,date_tmp))

    else:
        dataset=xr.open_dataset(dd)
        date=dataset['TIME'].values.flatten()
        u=dataset['UCUR'].values.flatten()
        v=dataset['VCUR'].values.flatten()
        p=dataset['PRES'].values.flatten()

    cutind=150
    output=xr.Dataset({'u':('date',u[cutind:]),'v':('date',v[cutind:]),'spd':('date',sqrt(u[cutind:]**2+v[cutind:]**2)),'prs':('date',p[cutind:])},coords={'date':date[cutind:]})
    return output

CF5=loadAQD('CF5')
CF6=loadAQD('CF6')
M1=loadAQD('M1')


CF5.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/AQD_cf5_500m.nc','w',format='netCDF4')
CF6.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/AQD_cf6_500m.nc','w',format='netCDF4')
M1.to_netcdf(datadir+'OSNAP2016recovery/gridded_CF-OOI/AQD_m1_500m.nc','w',format='netCDF4')

def plot_allprs():
    figure(figsize=(20,13))
    subplot(311)
    CF5.prs.plot()
    CF6.prs.plot()
    M1.prs.plot()
    subplot(312)
    CF5.u.plot()
    CF6.u.plot()
    M1.u.plot()
    subplot(313)
    CF5.v.plot()
    CF6.v.plot()
    M1.v.plot()


plot_allprs()


def plotall(field):
    figure(figsize=(12,3))
    for ii,xx in enumerate([CF5,CF6,M1]):
        if ii==2:
            dstep=4
        else:
            dstep=8
        Z,X = sig.butter(2,1./20/dstep/6, output='ba')
        xx_resamp=xx[field].rolling(date=dstep).std()
        xx_max=xx_resamp.resample(date='1D').mean()
        xx_max.plot(color='C'+str(ii),alpha=0.3)
        xx_sm=sig.filtfilt(Z,X,xx_resamp[~isnan(xx_resamp)].values)
        plot(xx_resamp.date[~isnan(xx_resamp)],xx_sm,linewidth=2,color='C'+str(ii),zorder=3)



plotall('spd')
ylim(0.02,0.065)

plotall('prs')
