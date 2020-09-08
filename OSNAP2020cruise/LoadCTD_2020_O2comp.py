#### SHOULD PROBABLY MAKE AN AR45 FUNCS FILE TO STORE SOME OF THESE DEETS

from AR45_funcs import *

figdir

#### Load all
bin_list=sort(glob.glob('/home/isabela/Documents/cruises/OSNAP2020_oxygen_AR45/data/ctd/ar45_*'))

def getdist(lon0,lat0,lonex,latex):
    distout=zeros(len(lonex))
    for ii in range(len(lonex)):
        distout[ii]=sw.dist([lat0,latex[ii]],[lon0,lonex[ii]])[0][0]
    return distout


###### Section numbering details + metadata
seclab={} #short for section key
# seclab['section 0']=array([0])
seclab['section 1']=range(1,29) # this is the section at the LS mooring line
seclab['section 2']=range(29,51)
seclab['section 3']=range(51,71)
seclab['section 3b']=range(71,86)
seclab['section 4 offshore']=range(86,94) # Includes CF7 location
seclab['CF tripods']=range(94,97)
seclab['section 4 onshore']=range(97,113) # Rest of the CF line
seclab['section 5']=range(113,137)
seclab['section EGC 1']=range(137,151)
seclab['section EGC 2']=range(151,163)

#### Onshore-most data point (to calc distance from)
startsec={}
startsec['section 1']=14
startsec['section 2']=43
startsec['section 3']=51
startsec['section 3b']=71
startsec['section 4 offshore']=93
startsec['CF tripods']=96
startsec['section 4 onshore']=97
startsec['section 5']=113
startsec['section EGC 1']=137
startsec['section EGC 2']=151

endsec={}
endsec['section 1']=1
endsec['section 2']=29
endsec['section 3']=70
endsec['section 3b']=85
endsec['section 4 offshore']=86
endsec['CF tripods']=94
endsec['section 4 onshore']=112
endsec['section 5']=136
endsec['section EGC 1']=150
endsec['section EGC 2']=162

#####Load it all into an xarray:
prslen=1750 #Nothing is deeper than 3500m
prsvec=range(2,2*prslen+1,2)

dat=fCNV(bin_list[0])

bin_list

def loadCTD_bin(datalist):
    stanum=len(datalist)
    lat=zeros(stanum)
    lon=zeros(stanum)
    date=[1]*stanum
    prs=NaN*ones((prslen,stanum))
    sal=NaN*ones((prslen,stanum))
    SA=NaN*ones((prslen,stanum))
    tmp=NaN*ones((prslen,stanum))
    CT=NaN*ones((prslen,stanum))
    den=NaN*ones((prslen,stanum))
    o2_sbe43=NaN*ones((prslen,stanum))
    v7=NaN*ones((prslen,stanum))
    o2_aa4831=NaN*ones((prslen,stanum))
    prsvec=array(range(2,2*prslen+1,2))
    flourmat=NaN*ones((prslen,stanum))
    maxprsvec=NaN*ones((stanum))
    for ii,dd in enumerate(datalist):
        dat=fCNV(dd)
        lat[ii]=dat.attributes['LATITUDE']
        lon[ii]=dat.attributes['LONGITUDE']
        date[ii]=dat.attributes['datetime']
        pl=len(dat['PRES'].data)
        maxprsvec[ii]=pl*2
        pst=int((dat['PRES'].data[0]-2)/2)
        prs[pst:pl+pst,ii]=dat['PRES'].data
        sal[pst:pl+pst,ii]=dat['PSAL'].data
        tmp[pst:pl+pst,ii]=dat['TEMP'].data
        flourmat[pst:pl+pst,ii]=dat['flECO-AFL'].data
        SA[pst:pl+pst,ii]=gsw.SA_from_SP(sal[:pl,ii],prs[:pl,ii],lon[ii],lat[ii])
        CT[pst:pl+pst,ii]=gsw.CT_from_pt(SA[:pl,ii],tmp[:pl,ii])
        den[pst:pl+pst,ii]=gsw.sigma0(SA[:pl,ii],CT[:pl,ii])
        ## convert from mL/L to uM (umol/l)
        ## 1uM = 0.022391 mL
        o2_sbe43[pst:pl+pst,ii]=dat['oxygen_ml_L']/0.022391
        #remove very small noisy voltages
        v7_clean=dat['v7'].copy()
        v7_clean[dat['v7']<0.1]=NaN
        v7[pst:pl+pst,ii]=v7_clean
        obase=aaoptode_ctdconvert(dat['PSAL'],dat['TEMP'],dat['PRES'],v7_clean)
        #remove noisy low bound aandera o2 values and interpolate over them (not sure why they persist)
        nanind=obase>265
        prel=prsvec[pst:pl+pst]
        ofunc=interpolate.interp1d(prel[nanind],obase[nanind],kind='linear',fill_value='NaN',bounds_error=False)
        o2_aa4831[:,ii]=ofunc(prsvec)

    distvec=array([NaN]*stanum)
    secvec=array([])
    for ss in seclab:
        secvec=hstack((secvec,[ss]*len(seclab[ss])))
        if len(seclab[ss])>1:
            #The -1 is there to deal with the fact that I am not using section 0!!
            news=[sss-1 for sss in seclab[ss]]
            distvec[news]=getdist(lon[startsec[ss]-1],lat[startsec[ss]-1],lon[news],lat[news])
        else:
            distvec[news]=zeros(1)

    grdat=xr.Dataset({'sal': (['prs','sta'],  sal),
                      'tmp': (['prs','sta'],  tmp),
                      'SA': (['prs','sta'],  SA),
                      'CT': (['prs','sta'],  CT),
                      'den': (['prs','sta'],  den),
                      'o2_sbe43': (['prs','sta'],  o2_sbe43),
                      'v7_o2': (['prs','sta'],  v7),
                      'o2_aa4831': (['prs','sta'],  o2_aa4831),
                      'odiff': (['prs','sta'], o2_aa4831-o2_sbe43),
                      'flourescence': (['prs','sta'], flourmat),
                      'max_prs': (['sta'],  maxprsvec),
                      'lat': (['sta'],  lat),
                      'lon': (['sta'],  lon),
                      'sec': (['sta'],  secvec),
                      'dist': (['sta'], distvec)},
                    coords={'sta': range(1,stanum+1),
                            'prs': prsvec})

    return grdat


grdat=loadCTD_bin(bin_list[1:])

datadir='/home/isabela/Documents/projects/OSNAP/data/Shipboard/netcdf/'
grdat.to_netcdf(datadir+'AllCTD20_2m_Leahproc.nc')

def plot_a_sec(var,secnum):
    figure()
    if 'o2' in var:
        varname='o2'
    else:
        varname=var
    pcolor(sort(grdat.dist.where(grdat.sec=='section '+secnum)),grdat.prs,(grdat.where(grdat.sec=='section '+secnum)[var]).sortby(grdat.dist.where(grdat.sec=='section '+secnum)),vmin=uni[varname]['vmin'],vmax=uni[varname]['vmax'],cmap=uni[varname]['cmap'])
    colorbar(label=var)
    ylim(3000,0)
    ylabel('pressure [db]')
    xlabel('distance [km]')
    if secnum=='1':
        title('West Greenland, LS line')
        savefig(figdir+'LSline_'+var+'.png',bbox_inches='tight')
    elif secnum=='4':
        title('East Greenland, CF line')
        savefig(figdir+'CFline_'+var+'.png',bbox_inches='tight')


colors = [(12,44,132),(78,179,211) ,(255,237,160),(217,95,14),(240,59,32)]
sal_cmap = make_cmap(colors,position=[0,0.9,0.96,0.99,1],bit=True)

uni['sal']={}
uni['sal']['cmap']=sal_cmap
uni['sal']['vmin']=32
uni['sal']['vmax']=35.2

uni['den']={}
uni['den']['cmap']=pden_cmap
uni['den']['vmin']=26
uni['den']['vmax']=28

uni['flourescence']={}
uni['flourescence']['cmap']=cm.Greens
uni['flourescence']['vmin']=0
uni['flourescence']['vmax']=6


# For each section, plot a pcolor section of sal, tmp, pden, O2 (aa +sbe), o2diff
# Then make profile comparisons for O2 (1 for each section - then a cumulative direct comp) using only deep sites for both
def plot_all_secs(secnum):
    for var in ['sal','tmp','den','o2_sbe43','o2_aa4831','odiff','flourescence']:
        plot_a_sec(var,secnum)

plot_all_secs('1')
plot_all_secs('4 offshore')
plot_all_secs('4 onshore')
plot_all_secs('5')

def comp_o2_profs(secnum):
    l1=plot(grdat.where(grdat.sec=='section '+secnum)['o2_sbe43'],grdat.prs,color='r',label='SBE43')
    l2=plot(grdat.where(grdat.sec=='section '+secnum)['o2_aa4831'],grdat.prs,color='b',label='AA4831')
    for ll in [l1,l2]:
        setp(ll[1:],label='')
    ylabel('pressure [db]')
    xlabel('O2 concentration [uM]')
    xlim(260,325)
    ylim(3000,0)
    legend()
    if secnum=='1':
        title('West Greenland, LS line')
        savefig(figdir+'LSline_o2comp_prof.png',bbox_inches='tight')
    elif ('4' in secnum) | ('5' in secnum):
        title('East Greenland, CF line')
        savefig(figdir+'CFline_o2comp_prof.png',bbox_inches='tight')

comp_o2_profs('1')
comp_o2_profs('4 offshore')

comp_o2_profs('5')

def o2diff_per_sec(secnum):
    plot(grdat.where(grdat.sec=='section '+secnum)['odiff'],grdat.prs,color='k')
    ylabel('pressure [db]')
    xlabel('AA4831 - SBE43\nO2 concentration difference [uM]')
    xlim(0,30)
    ylim(3000,0)
    if secnum=='1':
        title('West Greenland, LS line')
        savefig(figdir+'LSline_o2diff_prof.png',bbox_inches='tight')
    elif ('4' in secnum) | ('5' in secnum):
        title('East Greenland, CF line')
        savefig(figdir+'CFline_o2diff_prof.png',bbox_inches='tight')

o2diff_per_sec('1')
o2diff_per_sec('4 offshore')
o2diff_per_sec('5')

def o2diff_alldeep(stas,tit,lab):
    import matplotlib as mpl
    fig,ax=subplots(1,1)
    minsta=stas.min()
    maxsta=stas.max()
    elcmap=cm.rainbow
    for dsta in stas:
        ax.plot(grdat.where(grdat.sta==dsta)['odiff'],grdat.prs,color=elcmap((dsta-minsta)/(maxsta-minsta)))
    ylabel('pressure [db]')
    xlabel('AA4831 - SBE43\nO2 concentration difference [uM]')
    xlim(0,30)
    ylim(3000,0)
    title(tit)
    norm = mpl.colors.Normalize(vmin=minsta, vmax=maxsta)
    cax=fig.add_axes([0.925,0.15,0.02,0.7])
    cb1 = mpl.colorbar.ColorbarBase(cax,cmap=elcmap,norm=norm, label='Station number')
    savefig(figdir+lab+'_o2diff_prof.png',bbox_inches='tight')

o2diff_alldeep(grdat['sta'].where(grdat.max_prs>550),'All stations deeper than 550m','Alldeep')

o2diff_alldeep(grdat['sta'].where(grdat.sec=='section 1'),'West Greenland, LS line','LSline')
o2diff_alldeep(grdat['sta'].where(grdat.sec=='section 4 offshore'),'East Greenland, CF line','CFline')
o2diff_alldeep(grdat['sta'].where(grdat.sec=='section 5'),'East Greenland, just upstream of CF line','CFline')
