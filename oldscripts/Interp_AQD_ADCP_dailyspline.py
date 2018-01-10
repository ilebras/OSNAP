#############################################################################
#############################################################################
################### Interpolate AQD and ADCP vertically ######################
#############################################################################
#############################################################################

from aux_funcs import *

#############################################################################
################### Start with AQD data, then add ADCP and interp############
#############################################################################

# def Interp_AQD_ADCP(moornum):

moorname='CF7'
moornum=moorname[2]

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
        u_tmp=hrly_ave(dat.ix[:,6]/100,hrcon)
        v_tmp=hrly_ave(dat.ix[:,7]/100,hrcon)
        mprs_tmp=mean(dat.ix[:,5])*ones(len(date_tmp))

    else:
        dat = Dataset(dd, 'r')
        u_tmp=hrly_ave([float(tt) for tt in dat.variables['UCUR'][:].flatten()],hrcon*2)
        v_tmp=hrly_ave([float(tt) for tt in dat.variables['VCUR'][:].flatten()],hrcon*2)
        prs_tmp=hrly_ave(list(dat.variables['PRES'][:].flatten()),hrcon*2)
        time_tmp=list(dat.variables['TIME'][:])[::hrcon*2][:len(prs_tmp)]
        mprs_tmp=nanmean(prs_tmp)*ones(len(time_tmp))
        date_tmp=array([datetime.datetime(1950,1,1)+datetime.timedelta(days=int(tt))
                        for tt in time_tmp])

    prs=hstack((prs,prs_tmp))
    u=hstack((u,u_tmp))
    v=hstack((v,v_tmp))
    meanprs=hstack((meanprs,mprs_tmp))
    date=hstack((date,date_tmp))
    time=hstack((time,time_tmp))


    subplot(311)
    plot(date_tmp,prs_tmp)
    ylabel('pressure')
    subplot(312)
    plot(date_tmp,u_tmp)
    ylabel('u')
    subplot(313)
    plot(date_tmp,v_tmp)
    ylabel('v')


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
        print dd
        dat=pd.read_csv(dd,header=12,sep='\s*')
        prs_tmp=hrly_ave(gsw.p_from_z(-dat.ix[:,4],lat),hrcon)
        datea_tmp=unique([datetime.datetime(int(dat.ix[ii,0]),
                                  int(dat.ix[ii,1]),
                                  int(dat.ix[ii,2])) for ii in range(len(dat))])[:len(prs_tmp)]
        ua_tmp=hrly_ave(dat.ix[:,6]/100,hrcon)
        va_tmp=hrly_ave(dat.ix[:,7]/100,hrcon)
        prsa=hstack((prsa,prs_tmp))
        ua=hstack((ua,ua_tmp))
        va=hstack((va,va_tmp))
        datea=hstack((datea,datea_tmp))
        plot(datea_tmp,ua_tmp)
        plot(datea_tmp,va_tmp)

else:
    dat=io.loadmat(datadir+'ADCP_Data_CF/OSNAP_cf'+str(moornum)+'_Final1_ilebras.mat')
    shapmat=shape(dat['u'][1:,:])
    datea_tmp=unique([datetime.datetime(tt[0],tt[1],tt[2]) for tt in dat['time']])[:shapmat[1]/24]
#     Note that adcp "depths" from Dan Torres are also in db...
    prsa_tmp=zeros((shapmat[0],shapmat[1]/24))
    ua_tmp=zeros((shapmat[0],shapmat[1]/24))
    va_tmp=zeros((shapmat[0],shapmat[1]/24))
    for ii in range(shapmat[0]):
        prsa_tmp[ii,:]=hrly_ave(dat['z'][ii+1,:],hrcon)
        ua_tmp[ii,:]=hrly_ave(dat['u'][ii+1,:],hrcon)
        va_tmp[ii,:]=hrly_ave(dat['v'][ii+1,:],hrcon)

    ua=ua_tmp.flatten()
    va=va_tmp.flatten()
    prsa=prsa_tmp.flatten()

    datea=tile(datea_tmp,[1,shapmat[0]]).flatten()

print shape(prsa),shape(prs),shape(datea),shape(date)
pdall=pd.DataFrame({'pressure':hstack((prsa,prs)),
                               'u':hstack((ua,u)),
                               'v':hstack((va,v)),
                               'date bin':hstack((datea,date))})


# In[56]:


pdall['u'][pdall['u']<-2]=NaN
pdall['v'][pdall['v']<-2]=NaN


# In[58]:


pdall=pdall[pdall['pressure']>=0]


# In[59]:


common_mindate=max(min(datea),min(date))
common_maxdate=min(max(datea),max(date))


# In[60]:


common_maxdate


# In[61]:


# datevec=sort(unique(pdall['date bin']))
prsvec=arange(0,int(max(pdall['pressure']))+1,2)


# In[62]:


dmin=min(pdall['date bin'])
dmax=max(pdall['date bin'])
dlen=int(divmod((dmax-dmin).total_seconds(),60*60*24)[0])+1
datevec = array([dmin + datetime.timedelta(days=x) for x in range(0, dlen)])


# In[63]:


smoothparam=1e-2


# In[64]:


u_interp=pivspline_ML('u')


# In[65]:


v_interp=pivspline_ML('v')


# In[66]:


u_interp[(u_interp<-2) | (u_interp>2)]=nan
v_interp[(v_interp<-2) | (v_interp>2)]=nan


# In[67]:


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
savefig('../figures/interpolation/VEL/'+moorname+'/'+moorname+'_uv_measinterp.png',bbox_inches='tight')


# In[68]:


plotcontour(u_interp,cm.RdBu_r,-2,2)
savefig('../figures/interpolation/VEL/'+moorname+'/'+moorname+'_ucontour.png',bbox_inches='tight')


# In[69]:


plotcontour(v_interp,cm.RdBu_r,-2,2)
savefig('../figures/interpolation/VEL/'+moorname+'/'+moorname+'_vcontour.png',bbox_inches='tight')


# In[70]:


pickle.dump([u_interp,v_interp],open('../pickles/VELinterp/'+moorname+'_uvinterp.pickle','w'))
