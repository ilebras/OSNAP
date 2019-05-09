from AR30_funcs import *

datadir='/home/isabela/Documents/cruises/OSNAP2018_AR30-06/Chipods/AR30/Data/'

ctd=pickle.load(open('OSNAP2018cruise/pickles/CTD_1mbin.pickle','rb'))


ctd

stanum=len(ctd.sta)
zsm=20
fmax=20
respcorr=0
specstr='_256w'
specname=str(respcorr)+'corr_'+str(zsm)+'sm_'+str(fmax)+'Hz'+specstr
dat=io.loadmat(datadir+'AR30_XC_zsm'+str(zsm)+'m_fmax'+str(fmax)+'Hz_respcorr'+str(respcorr)+'_fc_99hz_gamma20'+specstr+'.mat')
#dtype=[('lat', 'O'), ('lon', 'O'), ('dnum', 'O'), ('TPvar', 'O'), ('eps', 'O'), ('chi', 'O'), ('KT', 'O'), ('dTdz', 'O'), ('N2', 'O'), ('t', 'O'), ('s', 'O'), ('P', 'O')]))],
#dtype=[('Name', 'O'), ('ChiDataDir', 'O'), ('MakeInfo', 'O'), ('BinInfo', 'O'), ('allSNs', 'O'), ('castnames', 'O'), ('SN2027_up_T1', 'O'), ('SN2027_down_T1', 'O'), ('SN2030_up_T1', 'O'), ('SN2030_down_T1', 'O')])

#ONLY LOOK AT THE UPCAST!!

#Indices will depend on whether I saved the downcast or not
upsave=0
if upsave==0:
    ind27=-2
    ind30=-1
if upsave==1:
    ind27=-4
    ind30=-2



data={}
data['SN2027']=dat['XC'][0][0][ind27][0][0]
# data['SN2027_dn']=dat['XC'][0][0][-3][0][0]
data['SN2030']=dat['XC'][0][0][ind30][0][0]
# data['SN2030_dn']=dat['XC'][0][0][-1][0][0]

sdv=['SN2027','SN2030']

ind27

# salex=data[sdv[0]][4]
chiprs=data[sdv[0]][-1].flatten()

maxprs=ones(stanum)*NaN

for ii in range(stanum):
    if sum(~isnan(ctd.sal[:,ii]))>1:
        maxprs[ii]=max(ctd.prs[~isnan(ctd.sal[:,ii])])

data[sdv[0]][4]

grdat=xr.Dataset({'epsilon': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][4],data[sdv[1]][4]))),
                  'chi': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][5],data[sdv[1]][5]))),
                  'KT': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][6],data[sdv[1]][6]))),
                  'dTdz': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][7],data[sdv[1]][7]))),
                  'N2': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][8],data[sdv[1]][8]))),
                  'tmp': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][9],data[sdv[1]][9]))),
                  'sal': (['prs','sta','sn_dir'],  dstack((data[sdv[0]][10],data[sdv[1]][10]))),
                  'tmp_1m': (['prs_ctd','sta'],ctd.tmp),
                  'sal_1m': (['prs_ctd','sta'],ctd.sal),
                  'den_1m': (['prs_ctd','sta'],ctd.den),
                  'deni_1m': (['prsi_ctd','sta'],ctd.deni),
                  'turner_1m': (['prsi_ctd','sta'],ctd.turner),
                  'turner_wden_1m': (['prsi_ctd','sta'],ctd.turner_wden),
                  'lat': (['sta'], ctd.lat),#data[sdv[0]][0][0]
                  'lon': (['sta'], ctd.lon),
                  'dist': (['sta'], ctd.dist),
                  'maxprs': (['sta'], maxprs)},#data[sdv[0]][1][0]
                coords={'sta': range(stanum),
                        'prs': chiprs,
                        'prs_ctd': ctd.prs.values,
                        'prsi_ctd': ctd.prsi.values,
                        'sn_dir': sdv})

grdat

grdat.prs[8:13]

def pline(var,ss,axxx):
    if ss=='section 5':
        distcorr1=20
    elif ss=='section 2':
        distcorr1=30

    highp=7

    if (var == 'chi') | (var=='epsilon'):
        # ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,log10(grdat[var][15,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T),color='b',linewidth=3)
        ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,log10(grdat[var][highp,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T),color='k',linewidth=3)
        # ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,log10(grdat[var][5,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T),color='g',linewidth=3)
        ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,log10(grdat[var][:15,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T))
        axxx.set_ylim(-10,-2)
    else:
        # ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,grdat[var][15,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T,color='b',linewidth=3)
        ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,grdat[var][highp,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T,color='k',linewidth=3)
        # ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,grdat[var][5,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir').T,color='g',linewidth=3)
        ret=axxx.plot(sort(grdat.dist[seclab[ss]])-distcorr1,grdat[var][:15,seclab[ss]].sortby(grdat.dist).mean(dim='sn_dir').T)
        axxx.set_ylim(33,35)

    return ret

def qfig():
    f,axi=subplots(4,1,sharex=True,figsize=(8,12))
    pline('sal','section 2',axi[0])
    pline('chi','section 2',axi[1])
    pline('sal','section 5',axi[2])
    pline('chi','section 5',axi[3])
    xlim(0,50)


qfig()

# def plotsections(var):
#     for ss in seclab:
#         f,axx=subplots(2,1,sharex=True,sharey=True,figsize=(10,8))
#         for ii in range(2):
#             px=axx[ii].pcolor(sort(grdat.dist[seclab[ss]]),grdat.prs,log10(grdat[var][:,seclab[ss],ii].sortby(grdat.dist)),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#             colorbar(px,ax=axx[ii])
#             axx[ii].set_title(sdv[ii])
#         # print(grdat.sta[seclab[ss]].sortby(grdat.dist).values)
#         # gca().invert_yaxis()
#         axx[1].set_xlabel('distance [km]')
#         # axx[3].set_ylim(500,0)
#         suptitle('Log('+str(var)+'), '+ss,fontsize=16)
#         # xlim(-5,max(grdat.dist[seclab[ss]])+10)
#         axx[1].set_ylim(max(maxprs[seclab[ss]]),0)
#         if var=='dTdz':
#             savefig(figdir+'AR30/Sections/'+var+'_'+ss+'.png',bbox_inches='tight')
#         else:
#             savefig(figdir+'AR30/Sections/'+var+'_'+ss+'_'+specname+'.png',bbox_inches='tight')
#
# plotsections('KT')
#
# plotsections('chi')
#
# plotsections('dTdz')

### For proposal: I want to make mean (between two sensors) sections of chi for those sections that cover the shelf to deep transition well.
# May also include T/S if there is room
# Sections to include:
# exsecs=['section 1', 'section 2', 'section 3', 'section 5']
# also, will want custom distance and depth ranges to put these on similar footing / draw out my point
bathy={}
bathy['section 5']=io.loadmat(datadir+'ar30_bathy/sect5_clean.mat')
bathy['section 2']=io.loadmat(datadir+'ar30_bathy/sect2_long_clean.mat')
bathy['section 3']=io.loadmat(datadir+'ar30_bathy/sect3_clean.mat')



uni['sal_1m']['vmin']=31
uni['sal_1m']['vmax']=35
uni['tmp_1m']['vmin']=1
uni['tmp_1m']['vmax']=11
uni['chi']['vmin']=-10
uni['chi']['vmax']=-5


uni['epsilon']={}
uni['epsilon']['vmin']=-9
uni['epsilon']['vmax']=-5
uni['epsilon']['cmap']=cm.magma_r
uni['chi']['cmap']=cm.magma_r

uni['sal']['levs']=array([33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35])
uni['tmp']['levs']=linspace(1,9,31)

def pvar(var,ss,axxx):
    if ss=='section 5':
        distcorr1=20
    elif ss=='section 2':
        distcorr1=30
    elif ss=='section 3':
        distcorr1=0

    if (var == 'chi') | (var=='epsilon'):
        ret=axxx.pcolor(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs,log10(grdat[var][:,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir')),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
    else:
        ret=axxx.contourf(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),
        levels=uni[var]['levs'],cmap=uni[var]['cmap'])
        # ret=axxx.pcolor(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),
        # vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
        axxx.contour(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat['den_1m'][:,seclab[ss]].sortby(grdat.dist),
        levels=[27.2,27.5,27.6,27.65,27.8],colors='k')
    axxx.fill_between(bathy[ss]['dist'].flatten()-distcorr1,bathy[ss]['bath'].flatten(),3e3,color='grey')

    return ret

def make_cbar(varg,tit,axxx,fig):
    p2 = axxx.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p2[0], 0.9, p2[2]-p2[0], 0.04])
    cbb=colorbar(varg, cax=ax_cbar, orientation='horizontal',ticklocation='top')
    cbb.set_label(tit)
    return cbb




def propsec():
    f,axx=subplots(2,4,sharey=True,sharex=True,figsize=(9,4))

    t1=pvar('tmp_1m','section 5',axx[1,0])
    c1=make_cbar(t1,'Temperature [$^\circ$C]',axx[0,0],f)
    c1.set_ticks(range(0,9,2))

    s1=pvar('sal_1m','section 5',axx[1,1])
    c2=make_cbar(s1,'Salinity',axx[0,1],f)
    c2.set_ticks(uni['sal']['levs'][::4])

    x1=pvar('chi','section 5',axx[1,2])
    make_cbar(x1,'log10($\chi$) [K$^2$ s$^{-1}$]',axx[0,2],f)

    e1=pvar('epsilon','section 5',axx[1,3])
    make_cbar(e1,'log10($\epsilon$) [W/kg$^{-1}$]',axx[0,3],f)

    pvar('tmp_1m','section 2',axx[0,0])
    pvar('sal_1m','section 2',axx[0,1])
    pvar('chi','section 2',axx[0,2])
    pvar('epsilon','section 2',axx[0,3])

    axx[0,0].set_xlim(0,53)
    axx[0,0].set_ylim(2000,0)

    fss=12
    f.text(0.5, 0, 'distance [km]', ha='center',fontsize=fss)
    f.text(0.05, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=fss)
    f.text(0, 0.75, 'Eastern Section', va='center',rotation='vertical',fontsize=fss+2)
    f.text(0, 0.25, 'Western Section',va='center', rotation='vertical',fontsize=fss+2)


propsec()
savefig('/home/isabela/Documents/proposals/2019Feb_OSNAP_mixing/OSNAPmixing-proposal/figures/PrelimSections.png',bbox_inches='tight',dpi=1e3)
savefig('/home/isabela/Documents/proposals/2019Feb_OSNAP_mixing/OSNAPmixing-proposal/figures/PrelimSections.pdf',bbox_inches='tight')

grdat

(grdat.sal[:,seclab['section 2'],:]).min()

uni['sal']['levs']=array([33, 34,  34.4,  34.8, 34.9, 34.92,34.94,34.96,34.98, 35])

def schemplot():
    f,ax1=subplots(2,1,figsize=(6,4),sharex=True,gridspec_kw={'height_ratios':[1,4]})
    var='sal_1m'
    ss='section 2'
    if ss=='section 5':
        distcorr1=20
    elif ss=='section 2':
        distcorr1=30
    elif ss=='section 3':
        distcorr1=0
    for axx in ax1:
        ret=axx.contourf(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),levels=uni[var]['levs'],cmap=uni[var]['cmap'],extend='both')
        axx.contour(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat['den_1m'][:,seclab[ss]].sortby(grdat.dist),levels=[27.2,27.5,27.6,27.65,27.8],colors='k',linewidth=2)
        axx.fill_between(bathy[ss]['dist'].flatten()-distcorr1,bathy[ss]['bath'].flatten(),3e3,color='grey')
    ax1[0].plot(sort(grdat.dist[seclab[ss]])-distcorr1, [0]*len(seclab[ss]),'ko')
    ax1[0].set_ylim(50,-5)
    ax1[1].set_ylim(500,50)
    ax1[1].set_ylabel('depth [m]')
    ax1[1].set_xlabel('distance [km]')
    xlim(0,80)
    cbax=f.add_axes([0.95,0.15,0.03,0.7])
    colorbar(ret,cax=cbax)

schemplot()


def scatterlines(var,axi,ss,prsmax,prsch):
    prange=(grdat.maxprs<prsmax)
    disss=sort(grdat.dist[seclab[ss]].where(prange))
    if nansum(sort(grdat.dist[seclab[ss]].where(prange)))>1:
        if '1m' in var:
            for ii in range(len(disss)):
                if 'deni' in var:
                    plti=axi.scatter(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange)[:,ii]+disss[ii]/mean(grdat[var][:,seclab[ss]].where(prange)),
                         prsch,c=grdat['turner_wden_1m'][:,seclab[ss]].sortby(grdat.dist).where(prange)[:,ii],vmin=-90,vmax=90,cmap=uni['turner']['cmap'],s=1)
                else:
                    if 'tmp' in var:
                        vmini=3
                        vmaxi=7
                    else:
                        vmini=uni[var[:3]]['vmin']
                        vmaxi=uni[var[:3]]['vmax']
                    plti=axi.scatter(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange)[:,ii]+disss[ii]/mean(grdat[var][:,seclab[ss]].where(prange)),
                         prsch,c=grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange)[:,ii],vmin=vmini,vmax=vmaxi,cmap=uni[var[:3]]['cmap'],s=1)
    return plti

def doublediffsec():
    f,axx=subplots(1,3,sharey=True,figsize=(10,2.5))
    ylim(150,0)
    prslim=800
    axx[0].set_ylabel('Depth [m]')

    psal=scatterlines('sal_1m',axx[0],'section 5',prslim,grdat.prs_ctd)
    make_cbar(psal,'Salinity',axx[0],f)
    axx[0].set_xticks([])

    ptmp=scatterlines('tmp_1m',axx[1],'section 5',prslim,grdat.prs_ctd)
    make_cbar(ptmp,'Temperature [$^\circ$C]',axx[1],f)
    axx[1].set_xticks([])

    pturn=scatterlines('deni_1m',axx[2],'section 5',prslim,grdat.prsi_ctd)
    cturn=make_cbar(pturn,'Turner angle',axx[2],f)
    axx[2].set_xticks(arange(26,27.5,0.5))
    axx[2].set_xticklabels(['26','26.5','27'],fontsize=10)
    cturn.set_ticks([-90,-45,0,45,90])


    axx[2].set_xlabel('Density [$\sigma_0$, kg/m$^3$]                                 ')

doublediffsec()
savefig('/home/isabela/Documents/proposals/OSNAP_mixing/OSNAPmixing-proposal/figures/doublediffsec.pdf',bbox_inches='tight')

#Check out staircases on shelf (no shift)
for ss in seclab:
    f,axx=subplots(2,2,sharey=True,figsize=(10,10))
    prslim=450
    lineplotsec('sal_1m',axx[0,0],ss,prslim,grdat.prs_ctd)
    lineplotsec('tmp_1m',axx[0,1],ss,prslim,grdat.prs_ctd)
    lineplotsec('den_1m',axx[1,0],ss,prslim,grdat.prs_ctd)
    lineplotsec('turner_1m',axx[1,1],ss,prslim,grdat.prsi_ctd)
    axx[0,0].set_ylabel('pressure [db]')
    axx[1,0].set_ylabel('pressure [db]')
    # axx[1].set_xlabel('distance [km]')
    suptitle(ss)
    ylim(100,0)
    savefig(figdir+'AR30/Profiles/StaircaseSearch_'+ss+'.png',bbox_inches='tight')



def getdist(lonex,latex):
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([latex[0],latex[ii]],[lonex[0],lonex[ii]])[0][0]
    return distout

vmadcp={}
vmadcp['section 2']=io.loadmat(datadir+'ADCP/vmadcp/vmpro/sect2os/sect2os_nb_dt.mat')['vm_data'][0][0]
vmadcp['section 5']=io.loadmat(datadir+'ADCP/vmadcp/vmpro/sect5os/sect5os_nb_dt.mat')['vm_data'][0][0]
#dt = de-tided
#[('cid', 'O'), ('sectid', 'O'), ('u', 'O'), ('v', 'O'), ('depth', 'O'), ('lon', 'O'), ('lat', 'O'), ('stID', 'O'), ('sampN', 'O'), ('u_dt', 'O'), ('v_dt', 'O'), ('DTmodel', 'O'), ('DTmodelZ', 'O'), ('DTmodelu', 'O'), ('DTmodelv', 'O'), ('trueBT', 'O')])

sect2dist=getdist(vmadcp['section 2']['lon'][0],vmadcp['section 2']['lat'][0])

pcolor(sect2dist-30,vmadcp['section 2']['depth'],-1*vmadcp['section 2']['v'],vmin=-0.6,vmax=0.6,cmap=cm.RdBu_r)
ylim(5e2,0)
xlim(0,53)
colorbar()

sect5dist=getdist(vmadcp['section 5']['lon'][0],vmadcp['section 5']['lat'][0])
pcolor(sect5dist-20,vmadcp['section 5']['depth'],vmadcp['section 5']['v'],vmin=-0.6,vmax=0.6,cmap=cm.RdBu_r)
ylim(5e2,0)
xlim(0,53)
colorbar()

plot(vmadcp['section 5'][0][0]['u'],'b');
plot(vmadcp['section 5'][0][0]['v'],'r');

plot(vmadcp['section 5'][0][0]['lon'].flatten()-360,vmadcp['section 5'][0][0]['lat'].flatten(),'o-');




# def propsec():
#     f,axx=subplots(3,2,sharey=True,sharex=True,figsize=(7,7))
#     uni['sal_1m']['vmin']=31
#     uni['sal_1m']['vmax']=35
#     uni['tmp_1m']['vmin']=1
#     uni['tmp_1m']['vmax']=9
#     uni['chi']['vmin']=-10
#     uni['chi']['vmax']=-6
#     ss='section 5'
#     var='chi'
#     distcorr1=10
#     px=axx[0,0].pcolor(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs,log10(grdat[var][:,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir')),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#     axx[0,0].fill_between(bathy[ss]['dist'].flatten()-distcorr1,bathy[ss]['bath'].flatten(),2e3,color='grey')
#     axx[0,0].set_title('Section '+ss[-1],fontsize=16)
#     var='tmp_1m'
#     px=axx[1,0].pcolor(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#     axx[1,0].fill_between(bathy[ss]['dist'].flatten()-distcorr1,bathy[ss]['bath'].flatten(),2e3,color='grey')
#     var='sal_1m'
#     px=axx[2,0].pcolor(sort(grdat.dist[seclab[ss]])-distcorr1,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#     axx[2,0].fill_between(bathy[ss]['dist'].flatten()-distcorr1,bathy[ss]['bath'].flatten(),2e3,color='grey')
#     xlim(0,62)
#     ylim(2000,0)
#     ss='section 2'
#     var='chi'
#     distcorr2=30
#     px=axx[0,1].pcolor(sort(grdat.dist[seclab[ss]])-distcorr2,grdat.prs,log10(grdat[var][:,seclab[ss],:].sortby(grdat.dist).mean(dim='sn_dir')),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#     axx[0,1].fill_between(bathy[ss]['dist'].flatten()-distcorr2,bathy[ss]['bath'].flatten(),2e3,color='grey')
#     c1=colorbar(px,ax=axx[0,1])
#     c1.set_label(label='log($\chi$) [$^\circ$C$^2$ s$^{-1}$]',fontsize=14)
#     axx[0,1].set_title('Section '+ss[-1],fontsize=16)
#     var='tmp_1m'
#     px=axx[1,1].pcolor(sort(grdat.dist[seclab[ss]])-distcorr2,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#     axx[1,1].fill_between(bathy[ss]['dist'].flatten()-distcorr2,bathy[ss]['bath'].flatten(),2e3,color='grey')
#     c2=colorbar(px,ax=axx[1,1])
#     c2.set_label(label='temperature [$^\circ$ C]',fontsize=14)
#     var='sal_1m'
#     px=axx[2,1].pcolor(sort(grdat.dist[seclab[ss]])-distcorr2,grdat.prs_ctd,grdat[var][:,seclab[ss]].sortby(grdat.dist),vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#     axx[2,1].fill_between(bathy[ss]['dist'].flatten()-distcorr2,bathy[ss]['bath'].flatten(),2e3,color='grey')
#     c3=colorbar(px,ax=axx[2,1])
#     c3.set_label(label='salinity',fontsize=14)
#     f.text(0.5, -0.01, 'distance [km]', ha='center',fontsize=14)
#     f.text(-0.01, 0.5, 'depth [m]', va='center', rotation='vertical',fontsize=14)
#     f.tight_layout()
#
# propsec()
# savefig('/home/isabela/Documents/proposals/OSNAP_Xpods/figures/PrelimSections_X.png',bbox_inches='tight')
# savefig('/home/isabela/Documents/proposals/OSNAP_Xpods/figures/PrelimSections_X.pdf',bbox_inches='tight')

XXXXXXXXXXXXXXX

#
# def OneProf(axi,var,ss,dirstr,prsch=grdat.prs):
#         if dirstr=='all':
#             sndirvec=[0,1,2,3]
#         elif dirstr=='up':
#             sndirvec=[0,2]
#         elif dirstr=='down':
#             sndirvec=[1,3]
#         elif dirstr=='2027':
#             sndirvec=[0,1]
#         elif dirstr=='2030':
#             sndirvec=[2,3]
#         # for tt in sndirvec:
#         #         axi.plot(log(grdat[var][:,seclab[ss],tt].sortby(grdat.dist)).where(grdat.maxprs>1500),prsch,label='')
#         pmax=300
#         if '1m' in var:
#             # for ii in range(len(seclab[ss])):
#             #     axi.plot(run_ave(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(grdat.maxprs>pmax)[:,ii],10),prsch)#takes too long to do in here. If/when I want to examine more closely, save a smooth version of turner.
#             axi.plot(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(grdat.maxprs>pmax),prsch,color='lightgrey')
#             axi.plot(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(grdat.maxprs>pmax).mean(dim='sta'),prsch,color='grey')
#             axi.plot(run_ave(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(grdat.maxprs>pmax).mean(dim='sta'),10),prsch,linewidth=2,color='k')
#             axi.set_title(var)
#             axi.set_xlim(uni[var]['vmin'],uni[var]['vmax'])
#             if 'turner' in var:
#                 xlabi=[-90,-45,0,45,72,90]
#                 [axi.axvline(xx,color='k') for xx in xlabi]
#                 axi.set_xticks(xlabi)
#                 axi.set_xlim(-120,120)
#         else:
#             for tt in sndirvec:
#                 axi.plot(log(grdat[var][:,seclab[ss],tt].sortby(grdat.dist)).where(grdat.maxprs>[pmax]).mean(dim='sta'),grdat.prs,linewidth=2,label=sdv[tt])
#                 axi.set_title('log('+var+')')
#                 axi.set_xlim(uni[var]['vmin'],uni[var]['vmax'])
#
#
# def plotprofs():
#     for ss in seclab:
#         f, axx = plt.subplots(1,5, sharey=True,figsize=(20,5))
#         dirstr='up'
#         OneProf(axx[0],'chi',ss,dirstr)
#         OneProf(axx[1],'N2',ss,dirstr)
#         OneProf(axx[2],'dTdz',ss,dirstr)
#         OneProf(axx[3],'KT',ss,dirstr)
#         OneProf(axx[4],'turner_1m',ss,dirstr,prsch=grdat.prsi_ctd)
#         gca().invert_yaxis()
#         ylim(2000,0)
#         suptitle(ss)
#         axx[0].legend(loc=[0.7,0.6])
#         savefig(figdir+'AR30/Profiles/ChiTurner_'+ss+'.png',bbox_inches='tight')
#
#
# plotprofs()
#
#
#
# def plotSalTmp():
#     for ss in seclab:
#         f, axx = plt.subplots(1,5, sharey=True,figsize=(20,5))
#         dirstr='all'
#         uni['sal']['vmin']=34.5
#         OneProf(axx[0],'N2',ss,dirstr)
#         OneProf(axx[1],'dTdz',ss,dirstr)
#         OneProf(axx[2],'sal_1m',ss,dirstr,prsch=grdat.prs_ctd)
#         OneProf(axx[3],'tmp_1m',ss,dirstr,prsch=grdat.prs_ctd)
#         OneProf(axx[4],'KT',ss,dirstr)
#         gca().invert_yaxis()
#         suptitle('Section '+str(ss))
#         axx[0].legend(loc=[0.7,0.6])
#
# plotSalTmp()
#
#
# # Note: should probably do mean in density space rather than pressure.
# def JustMean(axi,var,col,ss):
#         axi.plot(log(grdat[var][:,seclab[ss]].sortby(grdat.dist)).where(grdat.maxprs>1500).mean(dim='sta'),grdat.prs,linewidth=2,label=var,color=col)
#
#
# def comp_profs():
#     for ss in seclab:
#         f,ax1=subplots(1,1,figsize=(4,5))
#         JustMean(ax1,'chi','k',ss)
#         ax2=ax1.twiny()
#         JustMean(ax2,'N2','orange',ss)
#         ax3=ax1.twiny()
#         JustMean(ax3,'dTdz','red',ss)
#         gca().invert_yaxis()
#         ax1.set_xticklabels('')
#         ax2.set_xticklabels('')
#         ax3.set_xticklabels('')
#         title('Section '+ss)
#
#
# comp_profs()
#
# def pcolorsec(var,axi,ss,prsmax,prsch=grdat.prs):
#     prange=(grdat.maxprs<prsmax)
#     if nansum(sort(grdat.dist[seclab[ss]].where(prange)))>1:
#         if '1m' in var:
#             hh=axi.pcolor(sort(grdat.dist[seclab[ss]].where(prange)),prsch,
#                           grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange),
#                           vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#             if 'turner' in var:
#                 colorbar(hh,ax=axi,ticks=[-90,-45,0,45,90])
#             else:
#                 colorbar(hh,ax=axi)
#         else:
#             hh=axi.pcolor(sort(grdat.dist[seclab[ss]].where(prange)),prsch,log(grdat[var][:,seclab[ss],:].where(prange).mean(dim='sn_dir').sortby(grdat.dist)),
#                           vmin=uni[var]['vmin'],vmax=uni[var]['vmax'],cmap=uni[var]['cmap'])
#             colorbar(hh,ax=axi)
#         axi.set_title(var)
#
def lineplotsec(var,axi,ss,prsmax,prsch=grdat.prs):
    prange=(grdat.maxprs<prsmax)
    disss=sort(grdat.dist[seclab[ss]].where(prange))
    if nansum(sort(grdat.dist[seclab[ss]].where(prange)))>1:
        if '1m' in var:
            if 'turner' in var:
                axi.plot(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange),prsch,'-o',alpha=0.5)
                # axi.plot(run_ave(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange).mean(dim='sta'),10),prsch,'k',linewidth=3)
                xlim(-90,90)
                xticks(range(-90,91,45))
                [axvline(xx,color='darkgrey') for xx in [-90,-45,0,45,90] ]
            else:
                [axi.plot(grdat[var][:,seclab[ss]].sortby(grdat.dist).where(prange)[:,ii]+disss[ii]/mean(grdat[var][:,seclab[ss]].where(prange)),prsch) for ii in range(len(disss))]#
        else:
            [axi.plot(log(grdat[var][:,seclab[ss],:].where(prange).mean(dim='sn_dir').sortby(grdat.dist))[:,ii]+disss[ii]/2,prsch) for ii in range(len(disss))]
        axi.set_title(var)
#
# seclab['section 3']

#Check out staircases on shelf (no shift)
for ss in seclab:
    f,axx=subplots(2,2,sharey=True,figsize=(10,10))
    prslim=450
    lineplotsec('sal_1m',axx[0,0],ss,prslim,grdat.prs_ctd)
    lineplotsec('tmp_1m',axx[0,1],ss,prslim,grdat.prs_ctd)
    lineplotsec('den_1m',axx[1,0],ss,prslim,grdat.prs_ctd)
    lineplotsec('turner_1m',axx[1,1],ss,prslim,grdat.prsi_ctd)
    axx[0,0].set_ylabel('pressure [db]')
    axx[1,0].set_ylabel('pressure [db]')
    # axx[1].set_xlabel('distance [km]')
    suptitle(ss)
    ylim(100,0)
    savefig(figdir+'AR30/Profiles/StaircaseSearch_'+ss+'.png',bbox_inches='tight')

# # Hone in on slope. Line plots not too helpful here...
# for ss in seclab:
#     f,axx=subplots(1,4,sharey=True,figsize=(25,4))
#     prslim=2500
#     lineplotsec('sal_1m',axx[0],ss,prslim,grdat.prs_ctd)
#     lineplotsec('chi',axx[1],ss,prslim)
#     lineplotsec('dTdz',axx[2],ss,prslim)
#     lineplotsec('KT',axx[3],ss,prslim)
#     axx[0].set_ylabel('pressure [db]')
#     suptitle('Section '+ss)
#     axx[0].set_ylim(100,0)

uni['KT']['cmap']=cm.rainbow

for ss in seclab:
    f,axx=subplots(1,4,sharex=True,sharey=True,figsize=(25,4))
    prslim=3000
    pcolorsec('sal_1m',axx[0],ss,prslim,grdat.prs_ctd)
    pcolorsec('tmp_1m',axx[1],ss,prslim,grdat.prs_ctd)
    pcolorsec('turner_1m',axx[2],ss,prslim,grdat.prsi_ctd)
    pcolorsec('KT',axx[3],ss,prslim)
    axx[0].set_ylabel('pressure [db]')
    axx[1].set_xlabel('distance [km]')
    suptitle(ss)
    axx[3].set_ylim(max(maxprs[seclab[ss]]),0)
    savefig(figdir+'AR30/Sections/CompVars_'+ss+'.png',bbox_inches='tight')
#
# # # Check out the shelf... not much to see.. look for staircases instead.
# # for ss in seclab:
# #     f,axx=subplots(1,4,sharex=True,sharey=True,figsize=(25,4))
# #     prslim=400
# #     pcolorsec('sal_1m',axx[0],ss,prslim,grdat.prs_ctd)
# #     pcolorsec('tmp_1m',axx[1],ss,prslim,grdat.prs_ctd)
# #     pcolorsec('turner_1m',axx[2],ss,prslim,grdat.prsi_ctd)
# #     pcolorsec('KT',axx[3],ss,prslim)
# #     axx[0].set_ylabel('pressure [db]')
# #     axx[1].set_xlabel('distance [km]')
# #     suptitle('Section '+ss)
# #     ylim(prslim,0)
# #     xlim(0,60)
#
# grdat
# grdat
#
# # Plot chi and KT for different distances from the bottom.
#
# def botdistline(var,botz,axi,ss):
#     for tt in range(4):
#         for ii,dd in enumerate(seclab[ss]):
#             mp=maxprs[dd]
#             if ii==0:
#                 axi.plot(grdat.dist[dd],log(grdat[var][:,dd,tt].where(grdat.prs>mp-botz).mean(dim='prs')),'o',color=cm.get_cmap('rainbow',4)(tt),label=sdv[tt])
#             else:
#                 axi.plot(grdat.dist[dd],log(grdat[var][:,dd,tt].where(grdat.prs>mp-botz).mean(dim='prs')),'o',color=cm.get_cmap('rainbow',4)(tt))
#     axi.set_title('log('+var+')')
#     # axi.set_ylim(uni[var]['vmin']-1,uni[var]['vmax']+1)
#
# def maxdpthline(axi,ss):
#     for ii,dd in enumerate(seclab[ss]):
#         mp=maxprs[dd]
#         axi.plot(grdat.dist[dd],mp,'ko')
#         axi.set_title('bottom depth [m]')
#     axi.set_ylim(max(maxprs[seclab[ss]]),0)
#
# def plot_botdist(botz):
#     for ss in seclab:
#         if 'section' in ss:
#             f,axx=subplots(3,1,sharex=True,figsize=(5,11))
#             maxdpthline(axx[0],ss)
#             botdistline('chi',botz,axx[1],ss)
#             botdistline('KT',botz,axx[2],ss)
#             axx[1].legend(loc=(1.05,0.2))
#             axx[2].set_xlabel('distance [km]')
#             suptitle(ss,fontsize=18)
#             savefig(figdir+'AR30/Profiles/BottomProps_'+str(botz)+'m_'+ss+'.png',bbox_inches='tight')
#
# plot_botdist(20)
#
# plot_botdist(200)
