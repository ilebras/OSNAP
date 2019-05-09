from aux_funcs import *

versname='1810JHcal'
daily=pickle.load(open('../pickles/xarray/CF_xarray_gridplot_notid_'+versname+'.pickle','rb'))
daily['across track velocity']=-1*daily['across track velocity']

#################################################################################
###################### (Freshwater) Transport #################################
#################################################################################

mid_dist_plus=hstack((1.25,(diff(daily.distance)[:-1]+diff(daily.distance)[1:])/2,2.25))
middistmat_plus=transpose((tile(mid_dist_plus,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
mid_dist=mid_dist_plus.copy()
mid_dist[daily.distance<0]=0
mid_dist[daily.distance==0]=mid_dist_plus[daily.distance==0]/2
middistmat=transpose((tile(mid_dist,[len(daily.depth)-1,len(daily.date),1])),(2,0,1))
depthdiffmat=transpose((tile(diff(daily.depth),[len(daily.distance),len(daily.date),1])),(0,2,1))

daily['xport plus']=daily['across track velocity'][:,:-1,:]*depthdiffmat*middistmat_plus/1e3

onesxr=daily.salinity/daily.salinity

sep=9
cc={}
cc['trans']=daily['xport plus'][:sep,:,:].sum('depth').sum('distance')

cc['area']=(onesxr[:sep,:-1,:]*depthdiffmat[:sep,:,:]*middistmat[:sep,:,:]/1e3).sum('depth').sum('distance')
cc['sal_all']=(daily['xport plus'][:sep,:-1,:]*daily['salinity'][:sep,:-1,:]).sum('distance').sum('depth')/cc['trans']
cc['sal'][abs(cc['sal_all']-35)>3]=nanmean(cc['sal_all'])
cc['sal'].mean()

cc['sal_all'].plot()
cc['sal'].plot()

srefvec=arange(33,35.1,0.1)#[33,34,35,36]#29,30,31,32,
srefvec[-3]

cc_fresh=zeros((len(srefvec),len(daily.date)))
all_fresher=zeros((len(srefvec),len(daily.date)))
for ii,sref in enumerate(srefvec):
    cc_fresh[ii,:]=-(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-sref)/sref).sum('depth').sum('distance')
    all_fresher[ii,:]=-(daily['xport plus'][:,:-1,:].where(daily.salinity<=sref)*1e3*(daily.salinity[:,:-1,:]-sref)/sref).sum('depth').sum('distance')

def plot_refdep():
    plot(srefvec,mean(cc_fresh,axis=1),'-',label='Coastal current',linewidth=3)
    plot(srefvec,mean(all_fresher,axis=1),'-',label='Fresher than reference salinity',linewidth=3)
    legend(markerscale=2,fontsize=13)
    # plot(34.9,cc_fresh[-3,:].mean(),'ko')
    axvline(34.9,color='k')
    axhline(9,color='r')
    text(34,10,'Southeast Greenland',color='r',fontsize=13)
    text(34,5,'(Bamber et al. 2018)',color='r',fontsize=13)
    # text(33.15,43,'Referenced to 34.9',color='k',fontsize=13)
    # text(33.15,39,'(Le Bras et al. 2018)',color='k',fontsize=13)
    # grid('on')
    xlim(33.1,35)
    ylim(0,80)
    xlabel('Reference salinity')
    ylabel('Freshwater transport [mSv]')
    savefig('FreshwaterThoughts/refsaldep.png',bbox_inches='tight')

plot_refdep()


daily['salinity'][sep,:,:].mean(dim='depth').mean()

daily['salinity'][sep,:,:].mean(dim='depth').plot()

cc['fresh']=-(daily['across track velocity'][:sep,:-1,:]*depthdiffmat[:sep,:]*middistmat_plus[:sep,:]*(daily.salinity[:sep,:-1,:]-sref)/sref).sum('depth').sum('distance')
