## Post EGU re-doing all analyses

from aux_funcs import *

#mooring resolution gridded CF data
dat=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_notid_1808lpfilt.pickle','rb'))

ooi=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_HYPM_xray.pickle','rb'))

# Potential density for properties at 750 - most relevant for ISW
SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_pt(SA_vec,tmpvec)
pdenmat2=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat2[jj,ii]=gsw.pot_rho_t_exact(SA_vec[ii],tmpvec[jj],750,0)-1e3


denvec=arange(27.5,27.85,0.01)
middenvec=denvec[:-1]+diff(denvec)/2
varvec=['sal','tmp','dpth','thick']
denmats={}
for var in varvec:
    denmats[var]=zeros((len(dat.distance),len(middenvec),len(dat.date)))
    for moornum in range(8):
        for dd in range(len(middenvec)):
            if var=='dpth':
                denmats[var][moornum,dd,:]=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).min(dim='depth');
                minnan=dat['potential density'][moornum,:,:].min(dim='depth')>denvec[dd]
                denmats[var][moornum,dd,minnan]=NaN
            elif var=='thick':
                alldpths=dat.depth.where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1])
                denmats[var][moornum,dd,:]=alldpths.max(dim='depth')-alldpths.min(dim='depth');
            else:
                denmats[var][moornum,dd,:]=dat[univec[var][0]][moornum,:,:].where(dat['potential density'][moornum,:,:]>denvec[dd]).where(dat['potential density'][moornum,:,:]<=denvec[dd+1]).mean(dim='depth');



ooidenmat={}
for var in varvec:
    ooidenmat[var]=zeros((len(middenvec),len(ooi.date))) # note: OOI data is every 12 hours instead of 24
    for dd in range(len(middenvec)):
        if var=='dpth':
            ooidenmat[var][dd,:]=ooi['depth'].where(ooi['potential density']>denvec[dd]).min(dim='prs');
            minnan=ooi['potential density'].min(dim='prs')>denvec[dd]
            ooidenmat[var][dd,minnan]=NaN
        elif var=='thick':
            alldpths=ooi.depth.where(ooi['potential density']>denvec[dd]).where(ooi['potential density']<=denvec[dd+1])
            ooidenmat[var][dd,:]=alldpths.max(dim='prs')-alldpths.min(dim='prs');
        else:
            ooidenmat[var][dd,:]=ooi[univec[var][0]].where(ooi['potential density']>denvec[dd]).where(ooi['potential density']<=denvec[dd+1]).mean(dim='prs');

def makedendat():
    dendat=xr.Dataset({'depth': (['distance','den_bnd','date'],  denmats['dpth']),
                       'thickness': (['distance','den','date'],  denmats['thick']),
                       'sal': (['distance','den','date'],  denmats['sal']),
                       'tmp': (['distance','den','date'],  denmats['tmp'])},
                       coords={'den': middenvec,
                               'den_bnd': denvec[:-1],
                               'distance': dat.distance.values,
                               'date': dat.date.values})
# 'PV': (['distance','den','date'],  denmats['PV']),

    return dendat

dendat=makedendat()

ooi_dist=round(sw.dist([CFlat[0],float(ooi.lat[0])],[CFlon[0],float(ooi.lon)])[0][0])
def make_ooi_dendat():
    dendat=xr.Dataset({'depth': (['den_bnd','date'],  ooidenmat['dpth']),
                       'thickness': (['den','date'],  ooidenmat['thick']),
                       'sal': (['den','date'],  ooidenmat['sal']),
                       'tmp': (['den','date'],  ooidenmat['tmp'])},
                       coords={'den': middenvec,
                               'den_bnd': denvec[:-1],
                               'date': ooi.date.values,
                               'distance': ooi_dist})
# 'PV': (['distance','den','date'],  denmats['PV']),

    return dendat

ooi_dendat=make_ooi_dendat()


###############################################################################################################################
###############################################################################################################################
############################################## COLORS   ####################################################
###############################################################################################################################
###############################################################################################################################


m1col='C1'
cf5col='#ae017e'
ooicol='C0'


###############################################################################################################################
###############################################################################################################################
############################################## TS before and after   ####################################################
###############################################################################################################################
###############################################################################################################################


def TSmrange(salvar,tmpvar,axx,d1,d2,col,labit,shade='no'):

    if shade=='yes':
        smin=salvar.sel(date=slice(d1,d2)).min(dim='date').values.flatten()
        smax=salvar.sel(date=slice(d1,d2)).max(dim='date').values.flatten()
        # smin=(salvar.sel(date=slice(d1,d2)).mean(dim='date')-2*salvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
        # smax=(salvar.sel(date=slice(d1,d2)).mean(dim='date')+2*salvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
        sall=hstack((smin[::-1],smax))
        tmin=tmpvar.sel(date=slice(d1,d2)).min(dim='date').values.flatten()
        tmax=tmpvar.sel(date=slice(d1,d2)).max(dim='date').values.flatten()
        # tmin=(tmpvar.sel(date=slice(d1,d2)).mean(dim='date')-2*tmpvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
        # tmax=(tmpvar.sel(date=slice(d1,d2)).mean(dim='date')+2*tmpvar.sel(date=slice(d1,d2)).std(dim='date')).values.flatten()
        tall=hstack((tmin[::-1],tmax))
        p = plt.Polygon(np.column_stack((sall,tall)), facecolor=col, alpha=.25, edgecolor=col)
        axx.add_artist(p)
        axx.plot(smin,tmin,color=col)
        axx.plot(smax,tmax,color=col)

    tsm=axx.plot(salvar.sel(date=slice(d1,d2)).mean(dim='date').values.flatten(),tmpvar.sel(date=slice(d1,d2)).mean(dim='date').values.flatten(),
                 color=col,linewidth=3,label=labit)
    den=axx.contour(salvec[73:81],tmpvec,pdenmat2[:,73:81],colors='k',levels=[27.68,27.74,27.8])

    axx.set_ylim(3,5)
    axx.set_yticks(arange(3,5,0.5))
    axx.set_xticks(arange(34.86,35,0.04))
    axx.set_xlim(34.87,34.98)

    # clabel(den)

    return tsm,den


#
# def TS_BAshaded():
#     f,axx=subplots(1,2,sharex=True,sharey=True,figsize=(8,3))
#     #Before
#     b1='2014-9-10'
#     b2='2014-9-30'
#     l1,na=TSmrange(ooi_dendat['sal'],ooi_dendat['tmp'],axx[0],b1,b2,ooicol,'Irminger gyre interior, OOI')
#     l2,na=TSmrange(dendat['sal'][-1,:,:],dendat['tmp'][-1,:,:],axx[0],b1,b2,m1col,'Offshore of the boundary current, M1',shade='yes')
#     l3,na=TSmrange(dendat['sal'][4,:,:],dendat['tmp'][4,:,:],axx[0],b1,b2,cf5col,'Boundary current maximum, CF5',shade='yes')
#
#     a1='2015-4-10'
#     a2='2015-4-30'
#     TSmrange(ooi_dendat['sal'],ooi_dendat['tmp'],axx[1],a1,a2,ooicol,'')
#     TSmrange(dendat['sal'][-1,:,:],dendat['tmp'][-1,:,:],axx[1],a1,a2,m1col,'',shade='yes')
#     na,d3=TSmrange(dendat['sal'][4,:,:],dendat['tmp'][4,:,:],axx[1],a1,a2,cf5col,'',shade='yes')
#
#
#     axx[0].text(34.955,4.45,'uLSW',fontsize=14)
#     axx[0].text(34.95,3.85,'dLSW',fontsize=14)
#
#     t=axx[1].text(34.95,4.7,'27.68')
#     t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
#     t=axx[1].text(34.95,4.1,'27.74')
#     t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
#     t=axx[1].text(34.95,3.55,'27.8')
#     t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
#
#     axx[0].set_ylabel('Potential temperature [$^\circ$C]',fontsize=12)
#     f.text(0.5, -0.03, 'Salinity', ha='center',fontsize=12)
#     titf=14
#     axx[0].set_title('September 10-30, 2014\n (before convection)',fontsize=titf)
#     axx[1].set_title('April 10-30, 2015\n (after convection)',fontsize=titf)
#     axx[0].legend(loc=(-0.5,-0.4),ncol=3,fontsize=11)
#     savefig(figdir+'MixedLayer/paperfigs/TS_BeforeAfter.pdf',bbox_inches='tight')
#     savefig(figdir+'MixedLayer/paperfigs/TS_BeforeAfter.png',bbox_inches='tight')
#
# TS_BAshaded()

###############################################################################################################################
###############################################################################################################################
############################################## Layer thickness evolution   ####################################################
###############################################################################################################################
###############################################################################################################################

years=matplotlib.dates.YearLocator()
months=matplotlib.dates.MonthLocator()
threemonth=matplotlib.dates.MonthLocator(bymonthday=1,interval=3)
monthFMT=matplotlib.dates.DateFormatter('%B')
yearFMT=matplotlib.dates.DateFormatter('%Y \n ')

middenvec[17]

middenvec[17]
middenvec[18]
middenvec[23]
middenvec[29]

d1=27.65
d2=27.725
d3=27.76


#make an xarray which charters the depths of the bounding isopycnals at all cf + ooi
cf_bnds=xr.concat([dendat.depth[4:,18,:],dendat.depth[4:,24,:],dendat.depth[4:,30,:]],dim='den_bnd')
ooi_dpth_daily=ooi_dendat.depth.resample(date='1D').mean(dim='date')
ooi_bnds=xr.concat([ooi_dpth_daily[18,:],ooi_dpth_daily[24,:],ooi_dpth_daily[30,:]],dim='den_bnd')
all_bnds=xr.concat([cf_bnds,ooi_bnds],dim='distance')

figure(figsize=(12,3))
plot(ooi_dpth_daily.date,ooi_bnds[:2,:].T);
axvline(datetime.datetime(2015,4,14))
axvline(datetime.datetime(2015,4,20))
axvline(datetime.datetime(2015,4,26))



def ThickComp_Tseries(m1,m2,c1,c2,d1,d2,ax1,m3=0,c3=0):
        N  = 2    # Filter order
        Wn = 0.02 # Cutoff frequency (50 days)
        B, A = sig.butter(N, Wn, output='ba')

        hfa=0.3
        thick1=nansum(dendat.thickness[m1,d1:d2,:],axis=0)
        ax1.plot(dendat.date,thick1,alpha=hfa,color=c1)
        ax1.plot(dendat.date,sig.filtfilt(B,A,thick1),color=c1,linewidth=3)
        if m2=='ooi':
            # ax2=ax1.twinx()
            thick2=nansum(ooi_dendat.thickness[d1:d2,:],axis=0)
            ax1.plot(ooi_dendat.date,thick2,label='',color=c2,alpha=hfa)
            ax1.plot(ooi_dendat.date[::2],sig.filtfilt(B,A,thick2[::2]),color=c2,linewidth=3)
            ax1.set_ylim(400,1700)
            # ax1.set_yticks(range(500,2000,500))
            # ax2.set_ylim(400,1750)
            # def color_y_axis(ax, color):
            #     """Color your axes."""
            #     for t in ax.get_yticklabels():
            #         t.set_color(color)
            #     return None
            # color_y_axis(ax1, c1)
            # color_y_axis(ax2, c2)
        else:
            thick2=nansum(dendat.thickness[m2,d1:d2,:],axis=0)
            ax1.plot(dendat.date,thick2,label='',color=c2,alpha=hfa)
            ax1.plot(dendat.date,sig.filtfilt(B,A,thick2),color=c2,linewidth=3)
        if m3!=0:
            thick3=nansum(dendat.thickness[m3,d1:d2,:],axis=0)
            ax1.plot(dendat.date,thick3,label='',color=c3,alpha=hfa)
            ax1.plot(dendat.date,sig.filtfilt(B,A,thick3),color=c3,linewidth=3)


dlist=['2014-9-20','2015-2-20','2015-4-20','2015-9-20','2016-2-1','2016-4-1']
dtitlist=['Sep 20, 2014','Feb 20, 2015','Apr 20, 2015','Sep 20, 2015','Feb 1, 2016','Apr 1, 2016']
clist=['#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b']


def ThickEvol():

    f=figure(figsize=(20,6))
    gs = gridspec.GridSpec(2, 6)

    gs.update(hspace=0.5,wspace=0.3)

    ax1 = plt.subplot(gs[0, :3])
    ax2 = plt.subplot(gs[0, 3:])

    ThickComp_Tseries(4,7,cf5col,m1col,18,24,ax1)
    ThickComp_Tseries(7,'ooi',m1col,ooicol,24,30,ax2)

    ax1.set_ylabel('Layer thickness [m]')
    ax1.set_title('upper LSW \n',fontsize=18)
    ax2.set_title('deep LSW \n',fontsize=18)

    ax1.set_ylim(100,900)

    for axx in [ax1,ax2]:
        axx.xaxis.tick_top()
        axx.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        axx.xaxis.set_minor_formatter(monthFMT)
        axx.xaxis.set_major_formatter(yearFMT)

    # ax2.xaxis.tick_top()
    # ax2.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
    # ax2.xaxis.set_major_locator(years)
    # ax2.xaxis.set_minor_locator(threemonth)
    # ax2.xaxis.set_minor_formatter(monthFMT)
    # ax2.xaxis.set_major_formatter(yearFMT)

    [ax1.axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]
    [ax2.axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]


    for ii,d1 in enumerate(dlist):
        axx=plt.subplot(gs[1,ii])
        dmin=datetime.datetime.strptime(d1, '%Y-%m-%d')-datetime.timedelta(days=5)
        dmax=datetime.datetime.strptime(d1, '%Y-%m-%d')+datetime.timedelta(days=5)
        axx.plot(all_bnds.distance,all_bnds.sel(date=slice(dmin,dmax)).mean(dim='date').T,color=clist[ii],linewidth=3)
        axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
        axx.axvline(distvec[4],color=cf5col,linewidth=4)
        axx.axvline(distvec[7],color=m1col,linewidth=4)
        axx.axvline(ooi_dist,color=ooicol,linewidth=4)
        if ii!=0:
            gca().set_yticklabels('')
        else:
            ylabel('depth [m]',fontsize=12)
            text(120,700,'uLSW',fontsize=16)
            text(120,1300,'dLSW',fontsize=16)
        axx.set_title(dtitlist[ii])
        axx.set_ylim(2000,0)
        axx.set_xlim(30,190)
    f.text(0.5,-0.0025,'distance [km]', ha='center',fontsize=12)
    savefig(figdir+'MixedLayer/paperfigs/ThicknessEvol.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/ThicknessEvol.png',bbox_inches='tight')


ThickEvol()


def TS_triptich(axx,b1,b2):
    l1,na=TSmrange(ooi_dendat['sal'],ooi_dendat['tmp'],axx,b1,b2,ooicol,'Irminger gyre interior, OOI')
    l2,na=TSmrange(dendat['sal'][-1,:,:],dendat['tmp'][-1,:,:],axx,b1,b2,m1col,'Offshore of the boundary current, M1')
    l3,na=TSmrange(dendat['sal'][4,:,:],dendat['tmp'][4,:,:],axx,b1,b2,cf5col,'Boundary current maximum, CF5',shade='yes')
    return l1,l2,l3



def Thick_Tseries(m1,c1,d1,d2,ax1):
        N  = 2    # Filter order
        Wn = 0.02 # Cutoff frequency (50 days)
        B, A = sig.butter(N, Wn, output='ba')
        hfa=0.3

        if m1=='ooi':
            thick1=nansum(ooi_dendat.thickness[d1:d2,:],axis=0)
            ax1.plot(ooi_dendat.date,thick1,label='',color=c1,alpha=hfa)
            ax1.plot(ooi_dendat.date[::2],sig.filtfilt(B,A,thick1[::2]),color=c1,linewidth=3)
        else:
            thick1=nansum(dendat.thickness[m1,d1:d2,:],axis=0)
            ax1.plot(dendat.date,thick1,alpha=hfa,color=c1)
            ax1.plot(dendat.date,sig.filtfilt(B,A,thick1),color=c1,linewidth=3)

def ThickEvol_wTS():

    f=figure(figsize=(20,12))
    gs = gridspec.GridSpec(4, 6)

    gs.update(hspace=0.5,wspace=0.3)

    ax1 = plt.subplot(gs[0, :3])
    ax2 = plt.subplot(gs[0, 3:])

    Thick_Tseries('ooi',ooicol,18,24,ax1)
    Thick_Tseries('ooi',ooicol,24,30,ax2)

    ax11 = plt.subplot(gs[1, :3])
    ax22 = plt.subplot(gs[1, 3:])


    Thick_Tseries(4,cf5col,18,24,ax11)
    Thick_Tseries(7,m1col,18,24,ax11)

    Thick_Tseries(4,cf5col,24,30,ax22)
    Thick_Tseries(7,m1col,24,30,ax22)

    ax1.set_ylabel('Layer thickness [m]')
    ax11.set_ylabel('Layer thickness [m]')
    ax1.set_title('upper LSW \n',fontsize=18)
    ax2.set_title('deep LSW \n',fontsize=18)

    ax1.set_ylim(0,1200)
    ax2.set_ylim(400,2000)
    ax11.set_ylim(0,1100)
    ax22.set_ylim(0,1100)

    for axx in [ax1,ax2,ax11,ax22]:
        axx.xaxis.tick_top()
        axx.set_xlim([datetime.datetime(2014,8,1),datetime.datetime(2016,8,15)])
        axx.xaxis.set_major_locator(years)
        axx.xaxis.set_minor_locator(threemonth)
        if (axx==ax1) | (axx==ax2):
            axx.xaxis.set_minor_formatter(monthFMT)
            axx.xaxis.set_major_formatter(yearFMT)
        [axx.axvline(dd,color=clist[ii],linewidth=3) for ii,dd in enumerate(dlist)]

    ax11.set_xticklabels('')
    ax22.set_xticklabels('')

    for ii,d1 in enumerate(dlist):
        axx=plt.subplot(gs[2,ii])
        dmin=datetime.datetime.strptime(d1, '%Y-%m-%d')-datetime.timedelta(days=5)
        dmax=datetime.datetime.strptime(d1, '%Y-%m-%d')+datetime.timedelta(days=5)
        axx.plot(all_bnds.distance,all_bnds.sel(date=slice(dmin,dmax)).mean(dim='date').T,color=clist[ii],linewidth=3)
        axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
        axx.axvline(distvec[4],color=cf5col,linewidth=4)
        axx.axvline(distvec[7],color=m1col,linewidth=4)
        axx.axvline(ooi_dist,color=ooicol,linewidth=4)
        if ii!=0:
            gca().set_yticklabels('')
        else:
            ylabel('depth [m]',fontsize=12)
            text(120,700,'uLSW',fontsize=16)
            text(120,1300,'dLSW',fontsize=16)
        axx.set_title(dtitlist[ii])
        axx.set_ylim(2000,0)
        axx.set_xlim(30,190)

    for ii,d1 in enumerate(dlist):
        axx=plt.subplot(gs[3,ii])
        dmin=datetime.datetime.strptime(d1, '%Y-%m-%d')-datetime.timedelta(days=5)
        dmax=datetime.datetime.strptime(d1, '%Y-%m-%d')+datetime.timedelta(days=5)
        l1,l2,l3=TS_triptich(axx,dmin,dmax)
        if ii!=0:
            gca().set_yticklabels('')
        else:
            ylabel('pot. temperature [$^\circ$C]',fontsize=12)
            axx.text(34.95,4.45,'uLSW',fontsize=14)
            axx.text(34.945,3.85,'dLSW',fontsize=14)

    f.text(0.5,0.29,'distance [km]', ha='center',fontsize=12)
    f.text(0.5,0.08,'salinity', ha='center',fontsize=12)
    savefig(figdir+'MixedLayer/paperfigs/ThicknessEvol_wTS.pdf',bbox_inches='tight')
    savefig(figdir+'MixedLayer/paperfigs/ThicknessEvol_wTS.png',bbox_inches='tight')


ThickEvol_wTS()


###############################################################################################################################
###############################################################################################################################
############################################## Geostrophic velocity calc   #################################################
###############################################################################################################################
###############################################################################################################################


ooi_den_daily=ooi['potential density'].resample(date='1D').mean(dim='date')
ooi_dpth.min()

ooi_dpth=-gsw.z_from_p(ooi['prs'],60)
dmin=200
dmax=1975
grid_den={}
grid_dpth=range(dmin,dmax,25)
grid_den=zeros((2,len(grid_dpth),len(ooi_den_daily.date)))
grid_date=ooi_den_daily.date.values
for ii,dd in enumerate(grid_date):
    ## ooi
    ooiden=ooi_den_daily.sel(date=dd)
    if sum(~isnan(ooiden))>1:
        nini=~isnan(ooiden).values
        nanind=argwhere(nini)[0][0]
        # #pressures between which to calculate extrapolation gradient
        # z1=ooi_dpth[nanind]
        # z2=ooi_dpth[nanind+1]
        # # print(z1,min(prsnonan))
        # s1=ooiden[nanind]
        # s2=ooiden[nanind+1]
        # s0=s1-(s2-s1)*z1/(z2-z1)
        # fooi_1=interp1d(ooi_dpth[nini],ooiden[nini],bounds_error=False)
        # grid_den['ooi no extrap'][:,ii]=fooi_1(grid_dpth)
        # fooi_2=interp1d(hstack((0,ooi_dpth[nini])),hstack((s0,ooiden[nini])),bounds_error=False)
        # grid_den['ooi'][:,ii]=fooi_2(grid_dpth)
        fooi=interp1d(ooi_dpth[nini],ooiden[nini],bounds_error=False)
        grid_den[0,:,ii]=fooi(grid_dpth)
    else:
        grid_den[0,:,ii]=NaN*ones(len(grid_dpth))
        # grid_den['ooi no extrap'][:,ii]=NaN*ones(len(grid_dpth))

    ## m1
    m1den=dat['potential density'][-1,:,:].sel(date=dd)
    if sum(~isnan(m1den))>1:
        nini=~isnan(m1den).values
        nanind=argwhere(nini)[0][0]
        #pressures between which to calculate extrapolation gradient
        # z1=ooi_dpth[nanind]
        # z2=ooi_dpth[nanind+1]
        # # print(z1,min(prsnonan))
        # s1=ooiden[nanind]
        # s2=ooiden[nanind+1]
        # s0=s1-(s2-s1)*z1/(z2-z1)
        # fm1=interp1d(hstack((0,dat['depth'][nini])),hstack((s0,m1den[nini])),bounds_error=False)
        # grid_den['m1'][:,ii]=fm1(grid_dpth)
        fm1=interp1d(dat['depth'][nini],m1den[nini],bounds_error=False)
        grid_den[1,:,ii]=fm1(grid_dpth)
    else:
        grid_den[1,:,ii]=NaN*ones(len(grid_dpth))

datelen=len(grid_date)
## now go through every depth level and inerpolate horizontally
for jj,zz in enumerate(grid_dpth):
    ##ooi
    nani=~isnan(grid_den[0,jj,:])
    fooi_time=interp1d(arange(datelen)[nani],grid_den[0,jj,nani],bounds_error=False)
    grid_den[0,jj,:]=fooi_time(arange(datelen))

    ##m1
    nani=~isnan(grid_den[1,jj,:])
    fm1_time=interp1d(arange(datelen)[nani],grid_den[1,jj,nani],bounds_error=False)
    grid_den[1,jj,:]=fm1_time(arange(datelen))



f,[ax1,ax2]=subplots(1,2,sharey=True)#sharex=True,
ax1.plot(grid_den[0,:,:],grid_dpth);
ax2.plot(grid_den[1,:,:],grid_dpth);
ylim(2000,0)

fcor=gsw.f(60)
distdiff=90


shape(geoshear)

geoshear=squeeze(diff(grid_den,axis=0)*9.8/fcor/1028/distdiff/1e3)

geovel=cumsum(geoshear[::-1,:],axis=0)[::-1,:]*diff(grid_dpth)[0]


figure()
plot(geovel,grid_dpth);
plot(nanmean(geovel,axis=1),grid_dpth,'k-',linewidth=3);
ylim(2000,200)



figure(figsize=(12,3))
contourf(grid_date,grid_dpth,geoshear,levels=31,vmin=-0.0002,vmax=0.0002,cmap=cm.RdBu_r)
colorbar()

figure(figsize=(12,3))
contourf(grid_date,grid_dpth,geovel,levels=51,vmin=-0.05,vmax=0.05,cmap=cm.RdBu_r)
colorbar()

figure(figsize=(12,3))
plot(grid_date,sum(geovel,axis=0))
axhline(0)

XXXXXXXXXXX

###############################################################################################################################
###############################################################################################################################
############################################## Mean vel and density contours   #################################################
###############################################################################################################################
###############################################################################################################################
######## This may move to another script with a map!


gridded=pickle.load(open(datadir+'OSNAP2016recovery/pickles/xarray/CF_xarray_gridplot_notid_1810JHIL.pickle','rb'))
gridded['across track velocity']=-1*gridded['across track velocity']


def onecont(field,tit,vrange,coloor,hlevs,axx=0,nomoorlines=1,addcont=0):
    if axx==0:
        fig,axx=subplots(1,1)
    ax1=axx.contourf(gridded.distance,gridded.depth,field,vrange,cmap=coloor,extend='both')
    if addcont==1:
        ax2=contour(gridded.distance,gridded.depth,field,levels=hlevs,colors='k')
        clabel(ax2,fmt='%1.1f')
    axx.fill_between(bathdist,bathbath,2500*ones(len(bathbath)),color='k',zorder=22)
    xlabel('distance [km]',fontsize=16)
    ylabel('depth [m]',fontsize=16)
    xlim([25,190])
    ylim([2000,0])
    title(tit,fontsize=22)
    if nomoorlines==0:
        [axvline(mm,color='w',linewidth=2) for mm in distvec]
        [axvline(mm,color='k',linewidth=0.8) for mm in distvec]
    axx.set_yticks(range(500,2100,500))
    return ax1


#plot mean velocity and density levels for intro slides
def makesec():
        fig,axx=subplots(1,1,figsize=(5,4))
        field='uacross'
        axvel=onecont(gridded[univec[field][0]].mean(dim='date').T,'',univec[field][1][:-3],cm.GnBu,univec[field][3],axx=axx,nomoorlines=0)
        cbaxes = fig.add_axes([0.95, 0.15, 0.04, 0.7])
        cbar=colorbar(axvel,ticks=univec[field][3],label='[m/s]',cax=cbaxes)
        axx.plot(all_bnds.distance,all_bnds.mean(dim='date').T,color='k',linewidth=3)
        # axden=axx.contour(dat.distance,dat.depth,dat['potential density'].mean(dim='date').T,levels=[27.68,27.74,27.8],colors='k',linewidths=2)
        # manual_locations = [(70,500), (70, 1000), (70,1500)]
        # clabel(axden,fmt='%1.2f',manual=manual_locations)
        # axx.text(distvec[4]-30,-100,'Boundary current max',color='black',fontsize=18,zorder=200)
        axx.axvline(distvec[4],color='white',linewidth=9)
        axx.axvline(distvec[4],color=cf5col,linewidth=5)
        axx.axvline(distvec[-1],color='white',linewidth=9)
        # axx.text(distvec[-1]-10,-100,'Just offshore',color='black',fontsize=18,zorder=200)
        axx.axvline(distvec[-1],color=m1col,linewidth=5)
        axx.axvline(ooi_dist,color=ooicol,linewidth=5)
        savefig(figdir+'MixedLayer/paperfigs/Meansec_velden_placeholder.png',bbox_inches='tight')

makesec()


XXXXXXXXXXXXXXXXXX
