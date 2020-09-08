from firstfuncs_1618 import *

#### Redo this - get a new aux_funcs-like script that is updated for current needs (July 2020)
from AR45_funcs import *
uni={}
uni['sal']={}
uni['sal']['cmap']=cm.viridis
uni['sal']['vmin']=34.8
uni['sal']['vmax']=35.05
uni['ptmp']={}
uni['ptmp']['cmap']=cm.RdYlBu_r
uni['ptmp']['vmin']=0
uni['ptmp']['vmax']=6
uni['pden']={}
uni['pden']['cmap']=pden_cmap
uni['pden']['vmin']=27.3
uni['pden']['vmax']=27.9
uni['drhodz']={}
uni['drhodz']['cmap']=cm.Greens
uni['drhodz']['vmin']=0
uni['drhodz']['vmax']=0.001

uni['log(PPV)']={}
uni['log(PPV)']['cmap']=cm.plasma_r
uni['log(PPV)']['vmin']=-12.5
uni['log(PPV)']['vmax']=-9.5

datlist=sort(glob.glob(datadir+'aux_data/OOI/Irminger_fromJToole/*.mat'))

figdir='/home/isabela/Documents/projects/OSNAP/figures_OSNAPwide/Oxygen/OOI_HYPM/'


import h5py

grdat={}
for ii,dd in enumerate(datlist):
    f = h5py.File(dd,'r')
    #### For taking a look at the keys, etc
    dat={}
    for k, v in f.items():
        dat[k] = np.array(v)
    pydate=date_from_matlab(nanmean(dat['TIME'],axis=0))
    grdat[ii]=xr.Dataset({'sal':(['prs','date'],dat['S']),'ptmp':(['prs','date'],dat['THETA']),'pden':(['prs','date'],dat['SIGTH']),'u':(['prs','date'],dat['U']),'v':(['prs','date'],dat['V']),'w':(['prs','date'],dat['W']),},coords={'prs':dat['pgrid'].flatten(),'date':pydate})

rawdatlist=sort(glob.glob(datadir+'aux_data/OOI/HYPM_direct-OOI-download/*.nc'))

rawdat={}
for ii,dd in enumerate(rawdatlist):
    rawdat[ii]=xr.open_dataset(dd).swap_dims({'obs':'time'}).rename({'time':'date'})
    rawdat[ii]['ptmp']=('date',sw.ptmp(rawdat[ii]['practical_salinity'].values,rawdat[ii]['ctdpf_ckl_seawater_temperature'].values,rawdat[ii]['ctdpf_ckl_seawater_pressure'].values))




[ooi_lat,ooi_lon]=pickle.load(open(datadir+'OSNAP2016recovery/pickles/OOI/OOI_locs.pickle','rb'))

 # Calculate vertical density gradient, PPV (fN2), and "mixed layer"
for gg in grdat:
    grdat[gg]['pden_50m']=grdat[gg]['pden'].rolling(prs=100,center=True).mean()
    grdat[gg]['depth']=(['prs'],-gsw.z_from_p(grdat[gg]['prs'],ooi_lat['prf']))
    grdat[gg]['drhodz']=grdat[gg]['pden_50m'].diff(dim='prs')/grdat[gg]['depth'].diff(dim='prs')
    # Remove all points that have inversions:
    # grdat[gg]=grdat[gg].where(grdat[gg]['drhodz']>0)
    #calculate PPV
    grdat[gg]['PPV']=grdat[gg]['drhodz']*sw.f(ooi_lat['prf'])/(grdat[gg]['pden']+1e3)
    grdat[gg]['log(PPV)']=log10(grdat[gg]['drhodz']*sw.f(ooi_lat['prf'])/(grdat[gg]['pden']+1e3))

MLmat={}
# MLpden={}
# MLpden_var={}
# MLsal_var={}
# MLtmp_var={}
MLmat_T={}
MLmat_S={}
denthreshvec=[0.001,0.005,0.01]
for gg in grdat:
    print(gg)
    MLmat[gg]={}
    # MLpden[gg]={}
    # MLpden_var[gg]={}
    # MLsal_var[gg]={}
    # MLtmp_var[gg]={}
    MLmat_T[gg]={}
    MLmat_S[gg]={}
    # MLstd_S=NaN*zeros((len(grdat[gg].date),len(grdat[gg].prs)))
    # MLstd_T=NaN*zeros((len(grdat[gg].date),len(grdat[gg].prs)))
    pdeninterp=grdat[gg].pden.interpolate_na(dim='date')
    salinterp=grdat[gg].sal.interpolate_na(dim='date')
    tmpinterp=grdat[gg].ptmp.interpolate_na(dim='date')
    # for jj in range(len(grdat[gg].date)):
    #     MLstd_T[jj,:]=array([nanstd(grdat[gg].ptmp[int(mind/2):ii,jj]) for ii in range(len(grdat[gg].prs))])
    #     MLstd_S[jj,:]=array([nanstd(grdat[gg].sal[int(mind/2):ii,jj]) for ii in range(len(grdat[gg].prs))])
    for dd in denthreshvec:
        print(dd)
        MLmat[gg][dd]=NaN*zeros((len(grdat[gg].date)))
        # MLpden[gg][dd]=NaN*zeros((len(grdat[gg].date)))
        # MLpden_var[gg][dd]=NaN*zeros((len(grdat[gg].date)))
        # MLsal_var[gg][dd]=NaN*zeros((len(grdat[gg].date)))
        # MLtmp_var[gg][dd]=NaN*zeros((len(grdat[gg].date)))
        MLmat_T[gg][dd*10]=NaN*zeros((len(grdat[gg].date)))
        MLmat_S[gg][dd]=NaN*zeros((len(grdat[gg].date)))
        mind=200 #pressure calculating "mixed layer" from
        for jj in range(len(grdat[gg].date)):
            denvec=pdeninterp[grdat[gg].prs>=mind,jj]
            tmptmp=tmpinterp[grdat[gg].prs>=mind,jj]
            saltmp=salinterp[grdat[gg].prs>=mind,jj]
            if ~isnan(denvec[0]):
                diff_from_top=denvec-denvec[0]
                dft_S=abs(saltmp-saltmp[0])
                dft_T=abs(tmptmp-tmptmp[0])
                MLvec=grdat[gg].depth[grdat[gg].prs>=mind][diff_from_top>dd]
                MLvec_T=grdat[gg].depth[grdat[gg].prs>=mind][dft_T>dd*10]
                MLvec_S=grdat[gg].depth[grdat[gg].prs>=mind][dft_S>dd]
                if sum(MLvec)!=0:
                    MLmat[gg][dd][jj]=MLvec[0]
                if sum(MLvec_T)!=0:
                    MLmat_T[gg][dd*10][jj]=MLvec_T[0]
                if sum(MLvec_S)!=0:
                    MLmat_S[gg][dd][jj]=MLvec_S[0]
                    # MLpden[gg][dd][jj]=grdat[gg].pden[grdat[gg].prs>=mind,jj][diff_from_top>dd][0]
                    # MLpden_var[gg][dd][jj]=nanstd(grdat[gg].pden[grdat[gg].prs>=mind,jj][diff_from_top<=dd])
                    # MLsal_var[gg][dd][jj]=nanstd(grdat[gg].sal[grdat[gg].prs>=mind,jj][diff_from_top<=dd])
                    # MLtmp_var[gg][dd][jj]=nanstd(grdat[gg].ptmp[grdat[gg].prs>=mind,jj][diff_from_top<=dd])
                #make a new T/S based ML based on std thresholds instead!!
                # if sum(grdat[gg].depth[MLstd_T[jj,:]>dd*10])!=0:
                #     MLmat_T[gg][dd*10][jj]=grdat[gg].depth[MLstd_T[jj,:]>dd*10][0]
                # if sum(grdat[gg].depth[MLstd_S[jj,:]>dd])!=0:
                #     MLmat_S[gg][dd][jj]=grdat[gg].depth[MLstd_S[jj,:]>dd][0]


for gg in grdat:
    figure(figsize=(12,5))
    for dd in denthreshvec:
        plot(MLmat[gg][dd])
        title('density')
    figure(figsize=(12,5))
    for dd in denthreshvec:
        plot(MLmat_T[gg][dd*10])
        title('temperature')
    figure(figsize=(12,5))
    for dd in denthreshvec:
        plot(MLmat_S[gg][dd])
        title('salinity')


def plot_allML():
    f,axx=subplots(4,1,figsize=(15,12),sharey=True)
    for gg,axi in enumerate(axx):
        axi.plot(grdat[gg].date,MLmat[gg][0.005],label='density (0.005)')
        axi.plot(grdat[gg].date,MLmat_S[gg][0.005],label='salinity (0.005)')
        axi.plot(grdat[gg].date,MLmat_T[gg][0.05],label='temperature (0.05)')
        axi.set_ylabel('depth [m]')
        axi.set_xlim(datetime.datetime(2014+gg,7,10),datetime.datetime(2015+gg,8,15))
    axx[0].legend()
    axx[0].set_title('Mixed Layer with different definitions')
    savefig(figdir+'ML_alldefs.png',bbox_inches='tight')

plot_allML()

for gg in grdat:
    grdat[gg]['ML_den']=(['date'],MLmat[gg][0.005])
    grdat[gg]['ML_tmp']=(['date'],MLmat_T[gg][0.05])
    grdat[gg]['ML_sal']=(['date'],MLmat_S[gg][0.005])
    grdat[gg]['ML_TSD']=(['date'],xr.concat([grdat[gg]['ML_den'],grdat[gg]['ML_tmp'],grdat[gg]['ML_sal']],dim='ML_type').min(dim='ML_type'))


#Plot these as 4 panel figures showing each winter in the center
#Datasets go from
# 9/12/14 - 8/12/15
# 8/17/15 - 6/28/16
# 7/12/16 - 7/12/17
# 8/08/17 - 8/06/18

def plot_4_panel(var):
    f,axx=subplots(4,1,figsize=(12,12),sharey=True)
    f.subplots_adjust(hspace=0.5)
    for ii,axi in enumerate(axx):
        axi.contour(grdat[ii]['date'].values,grdat[ii]['prs'].values,grdat[ii]['pden'].values,colors='k',levels=[27.73,27.76],zorder=3)#My dISIW bounds were 27.73, 27.77
        grdat[ii][var].plot(ax=axi,vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'],zorder=1)
        # axi.contour(grdat[ii]['date'].values,grdat[ii]['prs'].values,grdat[ii]['pden'].values,colors='k',levels=[27.74,27.8],linestyles='--')#Pickart/Rhein dLSW bounds
        axi.set_xlim(datetime.datetime(2014+ii,7,10),datetime.datetime(2015+ii,8,15))
        if 'log(PPV)'==var:
            grdat[ii]['ML_den'].plot(ax=axi,color='darkblue',label='Mixed layer (density)')
            grdat[ii]['ML_TSD'].plot(ax=axi,color='C2', label='Mixed layer (with T and S)')
        axi.set_xlabel('')
        axi.set_ylabel('pressure [db]')
    axx[0].set_ylim(2625,100)
    axx[0].legend()
    savefig(figdir+'OOIHYPM_4panel_'+var+'.png',bbox_inches='tight')

plot_4_panel('log(PPV)')



for var in ['sal','ptmp','pden','drhodz','log(PPV)']:
    plot_4_panel(var)



##############################################################################################################################
##############################################################################################################################
#################################### Bin data in density space, look at TS
##############################################################################################################################
##############################################################################################################################
denvec=arange(27.5,27.85,0.005)
middenvec=denvec[:-1]+diff(denvec)/2
varvec=['sal','ptmp','dpth','thick']

def makedendat():
    dendat=xr.Dataset({'depth': (['den_bnd','date'],  denmats['dpth']),
                       'thickness': (['den','date'],  denmats['thick']),
                       'psal': (['den','date'],  denmats['sal']),
                       'ptmp': (['den','date'],  denmats['ptmp'])},
                       coords={'den': middenvec,
                               'den_bnd': denvec[:-1],
                               'date': grdat[gg].date.values})

    return dendat

dendat={}
for gg in grdat:
    denmats={}
    for var in varvec:
        denmats[var]=zeros((len(middenvec),len(grdat[gg].date)))
        for dd in range(len(middenvec)):
            if var=='dpth':
                denmats[var][dd,:]=grdat[gg].depth.where(grdat[gg]['pden'][:,:]>denvec[dd]).min(dim='prs');
                minnan=grdat[gg]['pden'].min(dim='prs')>denvec[dd] #
                denmats[var][dd,minnan]=0
            elif var=='thick':
                alldpths=grdat[gg].depth.where(grdat[gg]['pden']>denvec[dd]).where(grdat[gg]['pden']<=denvec[dd+1])
                denmats[var][dd,:]=alldpths.max(dim='prs')-alldpths.min(dim='prs');
            else:
                denmats[var][dd,:]=grdat[gg][var].where(grdat[gg]['pden']>denvec[dd]).where(grdat[gg]['pden']<=denvec[dd+1]).mean(dim='prs');

    dendat[gg]=makedendat()


colvec=['red','orange','blue','purple']
def TSmrange(da1,da2,gg,yr,colcol=0,lablab=0):
    if colcol==0:
        colcol=colvec[gg]
        lablab=str(yr+gg)
    # buf=0.02
    salsel=dendat[gg].psal.sel(date=slice(da1,da2))#.sel(den=slice(d1-buf,d3+buf))
    tmpsel=dendat[gg].ptmp.sel(date=slice(da1,da2))#.sel(den=slice(d1-buf,d3+buf))
    thicksel=dendat[gg].thickness.sel(date=slice(da1,da2))#.sel(den=slice(d1-buf,d3+buf))

    salvar=(salsel*thicksel).sum(dim='date')/thicksel.sum(dim='date')
    tmpvar=(tmpsel*thicksel).sum(dim='date')/thicksel.sum(dim='date')
    thickvar=thicksel.mean(dim='date')

    tsm=scatter(salvar.values.flatten(),
                 tmpvar.values.flatten(),thickvar.values.flatten()/2,
                 linewidth=3,label='',zorder=3,alpha=0.8,color=colcol)

    plot(salvar,tmpvar,color=colcol,label=lablab,linewidth=3)



def TS_layerbubs(yr,m1,m2,tit):
    figure(figsize=(6,5))
    for gg in range(3,-1,-1):
        d1=str(yr+gg)+'-'+str(m1)+'-01'
        d2=str(yr+gg)+'-'+str(m2)+'-01'
        TSmrange(d1,d2,gg,yr)
        contour(salvec,tmpvec,pdenmat2,colors='grey',linewidths=1,levels=denvec[::2])
        contour(salvec,tmpvec,pdenmat2,colors='k',linewidths=3,levels=[27.73,27.76])
        ylim(2.75,4.25)
        xlim(34.855,34.94)
        legend(loc=(1.05,0.3))
        xlabel('practical salinity')
        ylabel('pot. temperature [$^\circ$C]')
        title(tit,fontsize=14)
        savefig(figdir+'TS_layerthick_'+str(m1)+'-'+str(m2)+'.png',bbox_inches='tight')



# TS_layerbubs(2014,10,11,'Before convection (October)')
# TS_layerbubs(2014,10,12,'Before convection (October+November)')
# TS_layerbubs(2014,11,12,'(November)')
# TS_layerbubs(2015,1,2,'(January)')
# TS_layerbubs(2015,2,3,'Just before convection (February)')
# TS_layerbubs(2015,3,4,'During convection (March)')
# TS_layerbubs(2015,4,5,'During convection II (April)')
# TS_layerbubs(2015,6,7,'Just after convection (June)')
# TS_layerbubs(2015,5,6,'End of convection (May)')

import calendar
monthlist=[calendar.month_name[ii] for ii in hstack((range(7,13),range(1,9)))]
len(monthlist)

#Now I want to plot the TS progression by month for each deployment
def TS_layer_monthly():
    figure(figsize=(12,9))
    for gg in range(4):
        for ii in range(14):
            d1=datetime.datetime(2014+gg,7,1)+datetime.timedelta(days=30.4*ii)
            d2=datetime.datetime(2014+gg,8,1)+datetime.timedelta(days=30.4*ii)
            subplot(2,2,gg+1)
            TSmrange(d1,d2,gg,2014,colcol=cm.rainbow(ii/12),lablab=monthlist[ii])
            if gg==0:
                ylabel('pot. temperature [$^\circ$C]')
            if gg==2:
                legend(loc=(0.1,-0.35),ncol=7)
                ylabel('pot. temperature [$^\circ$C]')
            if gg>=2:
                xlabel('practical salinity')
        title(str(2014+gg)+'-'+str(2015+gg))
        contour(salvec,tmpvec,pdenmat2,colors='grey',linewidths=1,levels=denvec[::2])
        contour(salvec,tmpvec,pdenmat2,colors='k',linewidths=3,levels=[27.73,27.76])
        ylim(3,4)
        xlim(34.855,34.94)
    savefig(figdir+'TS_layerthick_allmonthly.png',bbox_inches='tight')

TS_layer_monthly()

### Quick comparison of recovered/deployed records:
def comp_end_beg(gg,enddates,begdates,comment,depshift=1,fm=''):
    figure(figsize=(12,4))
    subplot(1,2,1)
    plot(dendat[gg+depshift].psal.sel(date=begdates),dendat[gg+depshift].ptmp.sel(date=begdates),'.',color='r')
    plot(dendat[gg].psal.sel(date=enddates),dendat[gg].ptmp.sel(date=enddates),'.',color='b')
    plot(dendat[gg+depshift].psal.sel(date=begdates).mean(dim='date'),dendat[gg+depshift].ptmp.sel(date=begdates).mean(dim='date'),'-',color='r',label='first two weeks',linewidth=3)
    plot(dendat[gg].psal.sel(date=enddates).mean(dim='date'),dendat[gg].ptmp.sel(date=enddates).mean(dim='date'),'-',color='b',label='last two weeks',linewidth=3)
    xlabel('salinity')
    ylabel('pot. temperature [$^\circ$C]')
    legend()
    if depshift==1:
        suptitle(str(2015+gg)+comment)
    elif depshift==0:
        suptitle('Deployment '+str(1+gg)+comment)
    contour(salvec,tmpvec,pdenmat2,colors='grey',linewidths=1,levels=denvec[::2])
    contour(salvec,tmpvec,pdenmat2,colors='k',linewidths=3,levels=[27.73,27.76])
    ylim(3,4)
    xlim(34.82,34.94)
    subplot(1,2,2)
    contour(salvec,tmpvec,pdenmat2,colors='grey',linewidths=1,levels=denvec[::2])
    contour(salvec,tmpvec,pdenmat2,colors='k',linewidths=3,levels=[27.73,27.76])
    plot(rawdat[gg+depshift]['practical_salinity'].sel(date=begdates),rawdat[gg+depshift]['ptmp'].sel(date=begdates),'-',color='r',label='first two weeks (raw)',linewidth=1)
    plot(rawdat[gg]['practical_salinity'].sel(date=enddates),rawdat[gg]['ptmp'].sel(date=enddates),'-',color='b',label='last two weeks (raw)',linewidth=1)
    ylim(3,4)
    xlim(34.82,34.94)
    xlabel('salinity')
    legend()
    savefig(figdir+'DeploymentComparison_'+str(2015+gg)+'_dep'+str(depshift)+fm+'_withraw.png',bbox_inches='tight')


comp_end_beg(0,slice('2015-8-1','2015-8-15'),slice('2015-8-15','2015-8-31'),': 5 days between')

comp_end_beg(1,slice('2016-6-15','2016-6-30'),slice('2016-7-10','2016-7-25'),': two weeks between')
comp_end_beg(2,slice('2017-7-1','2017-7-15'),slice('2017-8-1','2017-8-15'),': two weeks between')

#now make similar plots but for the beginning and end of the deployment
#Datasets go from
# 9/12/14 - 8/12/15
comp_end_beg(0,slice('2015-8-1','2015-8-15'),slice('2014-9-15','2014-9-30'),': beg/end of deployment',0)
# 8/17/15 - 6/28/16
comp_end_beg(1,slice('2016-6-15','2016-6-30'),slice('2015-8-15','2015-8-30'),': beg/end of deployment',0)
# 7/12/16 - 7/12/17
comp_end_beg(2,slice('2017-7-1','2017-7-15'),slice('2016-7-15','2016-7-30'),': beg/end of deployment',0)
# 8/08/17 - 8/06/18

comp_end_beg(3,slice('2018-6-1','2018-6-15'),slice('2017-8-1','2017-8-15'),': beg/end of deployment',depshift=0)


#quickly, check if there is a significant adjustment in the first month
comp_end_beg(0,slice('2014-10-1','2014-10-15'),slice('2014-9-15','2014-9-30'),': first month of deployment',0,fm='firstmon')
# 8/17/15 - 6/28/16
comp_end_beg(1,slice('2015-9-1','2015-9-15'),slice('2015-8-15','2015-8-30'),': first month of deployment',0,fm='firstmon')
# 7/12/16 - 7/12/17
comp_end_beg(2,slice('2016-8-1','2016-8-15'),slice('2016-7-15','2016-7-30'),': first month of deployment',0,fm='firstmon')
comp_end_beg(2,slice('2016-9-15','2016-9-30'),slice('2016-7-15','2016-7-30'),': first two months of deployment',0,fm='first2mon')
comp_end_beg(2,slice('2017-6-1','2017-6-15'),slice('2017-5-1','2017-5-15'),': property shift in late May 2017',0,fm='propshift')
# 8/08/17 - 8/06/18

comp_end_beg(3,slice('2017-8-15','2017-9-1'),slice('2017-8-1','2017-8-15'),': first month of deployment',0,fm='firstmon')

comp_end_beg(3,slice('2017-9-15','2017-10-1'),slice('2017-8-1','2017-8-15'),': first two months of deployment',0,fm='firsttwomon')


uni['psal']=uni['sal']
uni['thickness']={}
uni['thickness']['cmap']=cm.inferno_r
uni['thickness']['vmin']=0
uni['thickness']['vmax']=750

def plot_4_panel_denspace(var):
    f,axx=subplots(4,1,figsize=(12,12),sharey=True)
    f.subplots_adjust(hspace=0.5)
    for ii,axi in enumerate(axx):
        dendat[ii][var].plot(ax=axi,vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'],zorder=1)
        axi.set_xlim(datetime.datetime(2014+ii,7,10),datetime.datetime(2015+ii,8,15))
        axi.set_xlabel('')
        axi.set_ylabel('pot. density [kg m$^{-3}$]')
        axi.axhline(27.73,color='k')
        axi.axhline(27.76,color='k')
        if 'anom' in var:
            axi.axvline(str(2014+ii)+'-10-01',color='k')
            axi.axvline(str(2014+ii)+'-10-31',color='k')
    axx[0].set_ylim(27.8,27.6)
    savefig(figdir+'OOIHYPM_4panel_denspace_'+var+'.png',bbox_inches='tight')

plot_4_panel_denspace('thickness')

uni['psal_anom']={}
uni['psal_anom']['cmap']=cm.RdBu_r
salanom=0.05
uni['psal_anom']['vmin']=-salanom
uni['psal_anom']['vmax']=salanom



uni['ptmp_anom']={}
uni['ptmp_anom']['cmap']=cm.RdBu_r
tmpanom=0.1
uni['ptmp_anom']['vmin']=-tmpanom
uni['ptmp_anom']['vmax']=tmpanom

uni['thickness_anom']={}
uni['thickness_anom']['cmap']=cm.RdBu_r
thickanom=500
uni['thickness_anom']['vmin']=-thickanom
uni['thickness_anom']['vmax']=thickanom


########Calculate anomalies from October 2014 average:
for dd in dendat:
    for var in dendat[dd]:
        dendat[dd][var+'_anom']=dendat[dd][var]-dendat[dd][var].sel(date=slice(str(2014+dd)+'-10-01',str(2014+dd)+'-10-31')).mean(dim='date')

plot_4_panel_denspace('ptmp_anom')

plot_4_panel_denspace('thickness_anom')

plot_4_panel_denspace('psal_anom')

###############################################################################################################
#####################################     Graveyard       #####################################################
###############################################################################################################
#
# At one point was looking at
# def checkout_Stds():
#     #standard deviation within ML plotted with Femke's thresholds
#     threshvec=[0.001,0.005,0.010000000000000002]
#     for gg in grdat:
#         figure()
#         title('ML depth')
#         for tt in threshvec:
#             plot(MLmat[gg][tt],'.')
#         axhline(0.05,color='k')
#         figure()
#         title('pot. density std')
#         for tt in threshvec:
#             plot(MLpden_var[gg][tt],'.')
#         axhline(0.05,color='k')
#         figure()
#         title('pot. temperature std')
#         for tt in threshvec:
#             plot(MLtmp_var[gg][tt],'.')
#         axhline(0.05,color='k')
#         figure()
#         title('salinity std')
#         for tt in threshvec:
#             plot(MLsal_var[gg][tt],'.')
#         axhline(0.005,color='k')
#
# checkout_Stds()
#
# MLxr={}
# for gg in grdat:
#     MLxr[gg]=xr.Dataset(coords={'date':grdat[gg].date.values})
#     for dd in denthreshvec:
#         MLxr[gg]['ML_'+str(dd)]=(['date'],MLmat[gg][dd])
#         MLxr[gg]['ML_T'+str(dd)]=(['date'],MLmat_T[gg][dd])
#         MLxr[gg]['ML_S'+str(dd)]=(['date'],MLmat_T[gg][dd])
#         MLxr[gg]['PDEN_'+str(dd)]=(['date'],MLpden[gg][dd])
#
# for gg in grdat:
#     # grdat[gg]['ML']=MLxr[gg]['ML_0.005'].where(MLtmp_var[gg][0.005]<0.05).where(MLsal_var[gg][0.005]<0.005).interpolate_na(dim='date')
#     figure()
#     MLxr[gg]['ML_0.005'].plot(marker='.',linewidth=0,label='pden only')
#     MLxr[gg]['ML_wTS_0.005'].plot(marker='.',linewidth=0,label='with TS')
#     legend()

# def TSplot():
#     for gg in grdat:
#         d1=str(2015+gg)+'-01-01'
#         d2=str(2015+gg)+'-05-01'
#         figure()
#         hexbin(grdat[gg].sal.sel(date=slice(d1,d2)).values.flatten(),grdat[gg].ptmp.sel(date=slice(d1,d2)).values.flatten(),mincnt=1,cmap=cm.hot,vmin=0,vmax=2000)
#         colorbar()
#         contour(salvec,tmpvec,pdenmat2,colors='b',linewidths=2,levels=[27.73,27.76])
#         ylim(2,5)
#         xlim(34.8,35)
#
# TSplot()
#
# def TSplot_means():
#     for gg in range(3,-1,-1):
#         d1=str(2015+gg)+'-03-01'
#         d2=str(2015+gg)+'-05-01'
#         plot(grdat[gg].sal.sel(date=slice(d1,d2)).mean(dim='date'),grdat[gg].ptmp.sel(date=slice(d1,d2)).mean(dim='date'),'.',label=str(2015+gg),color=colvec[gg],linewidth=2)
#         contour(salvec,tmpvec,pdenmat2,colors='k',linewidths=2,levels=[27.73,27.76])
#         ylim(2,5)
#         xlim(34.85,34.95)
#         legend()
#
# TSplot_means()
