from aux_funcs import *

adcp_dist=array(pd.DataFrame.from_csv(datadir+'Shipboard/adcp_distances.dat').index)

with open(datadir+'Shipboard/adcp_distances.dat', 'r') as f:
    reader = csv.reader(f)
    adcp_dist = list(reader)

adcp_dist=[float(dd[0]) for dd in adcp_dist]

adcp_14=io.loadmat(datadir+'Shipboard/kn221_2014/kn221_2014_vm_adcp_ilebras.mat')
adcp_16=io.loadmat(datadir+'Shipboard/ar07_2016/ar07_2016_vm_adcp_ilebras.mat')

adcp_14.keys()

plot(adcp_14['lon'],adcp_14['lat'],'ko');
plot(adcp_16['lon'],adcp_16['lat'],'ro');
plot(360+CFlon,CFlat,'yo');
ylim([60,60.2])
xlim([316.75,318])

def acrosstrackname(adic):
    adic['across track velocity']=adic['v_dt']*sin(theta)+adic['u_dt']*cos(theta)

    return adic

adcp_14=acrosstrackname(adcp_14)
adcp_16=acrosstrackname(adcp_16)


def contadcp(distvec,adic,tit):
    figure()
    contourf(distvec,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.8,0.05),cmap=cm.RdBu_r);
    colorbar(label='m/s')
    contour(distvec,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2],colors='k');
    fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    gca().invert_yaxis()
    title('Vessel mounted ADCP across-track velocity, '+tit)
    ylim([800,0])
    xlabel('distance (km)')
    ylabel('depth (m)')
    savefig('../figures/shipboard/ADCPsection_full_'+tit+'.pdf')
    xlim([-20,90])
    savefig('../figures/shipboard/ADCPsection_zoom_'+tit+'.pdf')





contadcp(adcp_dist,adcp_14,'2014')


contadcp(adcp_dist[:-3],adcp_16,'2016')

#################################################################################
# Make four-panel comparison plot between shipboard and mooring
#################################################################################

dat=pickle.load(open('../pickles/CF_xarray_gridplot_notid.pickle','rb'))

dat.date[0]

dat.date[-1]

bathdist=hstack((-20,bathdist))
bathbath=hstack((bathbath[0],bathbath))



#new function just has minimal framework
def contadcp_min(axspec,distvec,adic):
    axspec.contourf(distvec,adic['depth'][:,0],adic['across track velocity'],arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(distvec,adic['depth'][:,0],adic['across track velocity'],[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)

#same for mooring version
def contmoor_min(axspec,datechoice):
    im=axspec.contourf(dat.distance,dat.depth,dat['across track velocity'][:,:,datechoice].T,arange(-0.8,0.805,0.05),cmap=cm.RdBu_r);
    axspec.contour(dat.distance,dat.depth,dat['across track velocity'][:,:,datechoice].T,[-0.4,-0.2,0],colors='k');
    axspec.fill_between(bathdist,bathbath,1800*ones(len(bathbath)),color='k',zorder=22)
    [axspec.axvline(mm,color='w',linewidth=2) for mm in distvec]
    [axspec.axvline(mm,color='k',linewidth=0.8) for mm in distvec]

    return im


def plt4pan():
    fs=11
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8,6))
    contadcp_min(ax1,adcp_dist,adcp_14)
    ax1.set_xlim([-15,110])
    ax1.set_ylim([800,0])
    im=contmoor_min(ax2,0)
    contadcp_min(ax3,adcp_dist[:-3],adcp_16)
    im=contmoor_min(ax4,-1)
    ax1.set_title('Vessel mounted ADCP',fontsize=fs+2)
    ax2.set_title('Moored array measurements',fontsize=fs+2)
    ax1.set_ylabel('August, 2014 \n \n depth (m)',fontsize=fs)
    ax3.set_ylabel('July, 2016 \n \n depth (m)',fontsize=fs)
    # ax1.set_ylabel('depth (m)')
    # ax3.set_ylabel('depth (m)')
    ax3.set_xlabel('distance (km)',fontsize=fs)
    ax4.set_xlabel('distance (km)',fontsize=fs)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax,label='across-track velocity [m/s]')

    savefig('../figures/shipboard/ADCPmoorcomp_4panel.pdf')
    savefig('../figures/shipboard/ADCPmoorcomp_4panel.png')


plt4pan()

#################################################################################
# Compare vertically integrated velocity profiles as a function of distance from the coast
#################################################################################

def compdistvel(afunc,tit,ylab):
    intdpth=300
    plot(adcp_dist,afunc(adcp_14['across track velocity'][adcp_14['depth'][:,0]<intdpth,:],axis=0),label='August 2014')

    plot(adcp_dist[:-3],afunc(adcp_16['across track velocity'][adcp_16['depth'][:,0]<intdpth,:],axis=0),label='August 2016')
    plot(0,0,'k',label='Vessel mounted ADCP')
    plot(0,0,'k--',label='Moored array measurements')
    plot(dat.distance,afunc(dat['across track velocity'][:,dat.depth<intdpth,0],axis=1),'--',color='blue')
    plotinstpos('v')
    plot(dat.distance,afunc(dat['across track velocity'][:,dat.depth<intdpth,-1],axis=1),'--',color='orange')
    xlim([-15,110])
    axhline(0,color='k')
    xlabel('distance (km)')
    ylabel(ylab)
    legend(loc=(1.05,0.2))
    title(tit)

dat.salinity.min()

compdistvel(nanmean,'Average currents (to 300m)','velocity [m/s]')
savefig('../figures/shipboard/ADCPmoorcomp_avevel.png',bbox_inches='tight')

compdistvel(nansum,'Integrated currents (to 300m)','integrated velocity [m/s]')
savefig('../figures/shipboard/ADCPmoorcomp_intvel.png',bbox_inches='tight')
