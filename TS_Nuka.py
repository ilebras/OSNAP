from aux_funcs import *


### Testing the loading of Nuka data from Jamie



dat=io.loadmat(datadir+'aux_data/nuka/jamie/osnapsssfromnuka.mat')


plot(dat['OSNAP_East'][:,0],dat['OSNAP_East'][:,1],'o')
plot(CFlon,CFlat)

# Load all large TS dictionaries

[date,month,prs,sal,tmp]=pickle.load(open('../pickles/TSdailydic/TS_daily_dic_wcorr.pickle','rb'))


# set up basic ts plotting framework

def spruce_TS():
    dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k',zorder=20)
    clabel(dens)
    ylabel('temperature')
    xlabel('salinity')
    legend(loc=(1.05,0.2),numpoints=1,markerscale=2)


def TSplot_seasonal(moor,prs):
    scatter(sal[moor][prs],tmp[moor][prs],c=month[moor][prs],
            cmap=cm.brg_r,lw=0,alpha=0.5,vmin=1,vmax=12)


def spruce_seasonal():
    grid('on')
    cbar=colorbar(ticks=range(1,13,2))
    cbar.set_ticklabels(['January','March','May','July','September','November'])
    ylabel('temperature')
    xlabel('salinity')
    dens=contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k')
    clabel(dens)


def TSall(col='grey',alf=0.5):
    for key in sal:
        for key2 in sal[key]:
            plot(sal[key][key2],tmp[key][key2],
                 '.',color=col,alpha=alf,
                mec=None)

for key in range(1,8):
    figure(figsize=(12,8))
    TSall()
    spruce_TS()
    for key2 in sal[key]:
        plot(sal[key][key2],tmp[key][key2],'r.')
    plot(dat['osnapsi'][key,:],dat['osnapti'][key,:],'bo')
    title('CF'+str(key))
    savefig('../figures/TS/CF'+str(key)+'_TS_all.png')
    savefig('../figures/TS/CF'+str(key)+'_TS_all.pdf')

for key in month:
    for key2 in month[key]:
        month[key][key2]=array(month[key][key2])


# def TSmonthly():
for key in range(1,8):
    print(key)
    fig, axarr = subplots(3,4, sharex=True, sharey=True, figsize=(20,15))
    for ii in range(12):
        for keyall in sal:
            for key2all in sal[keyall]:
                axarr[int(ii/4),mod(ii,4)].plot(sal[keyall][key2all],
                                             tmp[keyall][key2all],
                     '.',color='grey',alpha=0.3,
                    mec=None)
        for key2 in sal[key]:
            yrvec=array([dd.year for dd in date[key][key2][month[key][key2]==ii+1]])
            im=axarr[int(ii/4),mod(ii,4)].scatter(sal[key][key2][month[key][key2]==ii+1],
                    tmp[key][key2][month[key][key2]==ii+1],
                    c=yrvec,cmap=cm.rainbow,vmin=2014,vmax=2016,lw=0,zorder=20)
            axarr[int(ii/4),mod(ii,4)].plot(dat['osnapsi'][key,ii],
                                        dat['osnapti'][key,ii],'bo',markersize=20,zorder=30)

        axarr[int(ii/4),mod(ii,4)].set_title(ii+1)
        dens=axarr[int(ii/4),mod(ii,4)].contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k',zorder=20)

        if (ii+1)==8:
            cax = fig.add_axes([1, 0.15, 0.03, 0.7])
            cbar=fig.colorbar(im, cax=cax,ticks=[2014,2015,2016])
            cbar.ax.set_yticklabels(['2014','2015','2016'])
            axarr[int(ii/4),mod(ii,4)].set_ylim(-2,10)
            axarr[int(ii/4),mod(ii,4)].set_xlim(31,36)
    fig.text(0.5, 1.04,'CF'+str(key),
             fontsize=20,ha='center')
    fig.text(0.5, -0.04, 'salinity', ha='center')
    fig.text(-0.04, 0.5, 'temperature', va='center', rotation='vertical')
    fig.tight_layout()
    savefig('../figures/TS/monthly/CF'+str(key)+'_TS_monthly.png',
                bbox_inches='tight')

# def TSmonthly_prs():
#     for key in range(1,8):
#         print(key)
#         fig, axarr = subplots(3,4, sharex=True, sharey=True, figsize=(20,15))
#         for ii in range(12):
#             for keyall in sal:
#                 for key2all in sal[keyall]:
#                     axarr[int(ii/4),mod(ii,4)].plot(sal[keyall][key2all],
#                                                  tmp[keyall][key2all],
#                          '.',color='grey',alpha=0.3,
#                         mec=None)
#             for key2 in sal[key]:
#                 im=axarr[int(ii/4),mod(ii,4)].scatter(sal[key][key2][month[key][key2]==ii+1],
#                         tmp[key][key2][month[key][key2]==ii+1],
#                         c=prs[key][key2][month[key][key2]==ii+1],cmap=cm.rainbow_r,
#                         vmin=0,vmax=1000,lw=0,zorder=20)
#                 axarr[int(ii/4),mod(ii,4)].plot(dat['osnapsi'][key,ii],
#                                             dat['osnapti'][key,ii],'bo',markersize=20,zorder=30)
#
#             axarr[int(ii/4),mod(ii,4)].set_title(ii+1)
#             dens=axarr[int(ii/4),mod(ii,4)].contour(salvec,tmpvec,pdenmat,arange(24,29,0.5),colors='k',zorder=20)
#             axarr[int(ii/4),mod(ii,4)].grid('on')
#
#             if (ii+1)==8:
#                 cax = fig.add_axes([1, 0.15, 0.03, 0.7])
#                 cbar=fig.colorbar(im, cax=cax,ticks=[0,500,1e3])
#                 axarr[int(ii/4),mod(ii,4)].set_ylim(-2,10)
#                 axarr[int(ii/4),mod(ii,4)].set_xlim(31,36)
#         fig.text(0.5, 1.04,'CF'+str(key),
#                  fontsize=20,ha='center')
#         fig.text(0.5, -0.04, 'salinity', ha='center')
#         fig.text(-0.04, 0.5, 'temperature', va='center', rotation='vertical')
#         fig.tight_layout()
#         savefig('../figures/TS/monthly/CF'+str(key)+'_TS_monthly_prs.png',
#                     bbox_inches='tight')



#########################################################################
################## Load all Nuka data #######################################
########################################################################

nukall=io.loadmat(datadir+'aux_data/nuka/jamie/nukaosnapdata.mat')

nukall.keys()

nukall['osnapx']

plot(nukall['x'],nukall['y'],'.');
plot(nukall['osnapx'],nukall['osnapy'],'k*',markersize=10);



shape(nukall['x'])
shape(nukall['saltracks'])

TSall()
scatter(nukall['saltracks'],nukall['temptracks'],'.');

barrier

#########################################################################
################## Sermilik data #######################################
########################################################################

#
# sermlist=glob.glob(datadir+'/aux_data/sermilik2015ctds/*')
#
#
# sermT=[]
# sermS=[]
# sermP=[]
# for ss in sermlist:
#     sermdat=io.loadmat(ss)
#     sermS=hstack((sermS,sermdat['sal'].flatten()))
#     sermT=hstack((sermT,sermdat['temp'].flatten()))
#     sermP=hstack((sermP,sermdat['pres'].flatten()))
#
#
# TSall()
# plot(sermS,sermT,'r.',alpha=0.5)
# xlim([31,36])
# spruce_TS()
# title('Sermilik and OSNAP TS space')
# savefig('../figures/TS/Sermilik_comp.pdf')
# savefig('../figures/TS/Sermilik_comp.pdf')


#########################################################################
################## Feili's data #######################################
########################################################################

#
#
# feili=io.loadmat(datadir+'Feili/OSNAP_CF_regridded_170613_Feili_resaved.mat')
#
#
# plot(feili['sal'].flatten(),feili['tmp'].flatten(),'r.')
# TSall(col='k',alf=1)
# # spruce_TS()
# grid('on')
# xlabel('salinity')
# ylabel('temperature')
# savefig('../figures/TS/Feili_comp.pdf')
# savefig('../figures/TS/Feili_comp.png')
#
# #doesn't actually make sense to compare here until on same grid (see Inspect_saltmpinterp.py)
#
# hexbin(feili['sal'].flatten(),feili['tmp'].flatten(),bins='log',cmap=cm.hot_r)
# colorbar()
# spruce_TS()
# xlim([31,36])
# ylim([-2,10])
# savefig('../figures/TS/TS_total_log_feili.png')
# savefig('../figures/TS/TS_total_log_feili.pdf')
