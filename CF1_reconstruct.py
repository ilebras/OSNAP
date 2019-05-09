#############################################################################
##############################################################################
#######  Reconstruct CF1 time series post icey blowdown ###############
#############################################################################
#############################################################################

from aux_funcs import *

moorname='CF1'

# Load corrected "daily dic" version

[date_all,month_all,prs_all,sal_all,tmp_all]=pickle.load(open('../pickles/TSdailydic/TS_daily_dic_wcorr.pickle','rb'))

date=date_all[1]
tmp=tmp_all[1]
sal=sal_all[1]
prs=prs_all[1]
month=month_all[1]

#############################################################################
################## First look at data #################
#############################################################################


figure(figsize=(13,10))
for ii in sal:
    subplot(311)
    plot(date[ii],sal[ii])
    subplot(312)
    plot(date[ii],tmp[ii])
    subplot(313)
    plot(date[ii],prs[ii])

for key in date:
    print(date[key][-1])
#############################################################################
################## Look at linear regressions between instruments #################
#############################################################################

stind=0
for ii in sal:
  if ii!=140:
    if ii!=150:
        timeind=(date[150]<=max(date[ii])) & (date[150]>=min(date[ii]))
        sal2=sal[150][timeind]
        tmp2=tmp[150][timeind]
        sal_line=poly1d(polyfit(sal2,sal[ii][stind:], 1))
        tmp_line=poly1d(polyfit(tmp2,tmp[ii][stind:], 1))

        figure(1,figsize=(10,4))
        subplot(121)
        plot(sal[150][timeind],sal[ii][stind:],'o',alpha=0.5)
        plot(sort(sal2),sal_line(sort(sal2)),linewidth=2,zorder=8)
        xlim([32.5,35])
        ylim([32.5,35])
        subplot(122)
        plot(tmp[150][timeind],tmp[ii][stind:],'o',alpha=0.5,label=str(int(mean(prs[ii]))))
        plot(sort(tmp2),tmp_line(sort(tmp2)),linewidth=2,zorder=8)
        xlim([-2,7])
        ylim([-2,7])

        figure(figsize=(10,4))
        subplot(121)
        scatter(sal[150][timeind],sal[ii][stind:],c=month[ii][stind:],vmin=1,vmax=12)
        plot(sort(sal2),sal_line(sort(sal2)),'k',linewidth=2)
        xlim([32.5,35])
        ylim([32.5,35])
        subplot(122)
        scatter(tmp[150][timeind],tmp[ii][stind:],c=month[ii][stind:],vmin=1,vmax=12)
        plot(sort(tmp2),tmp_line(sort(tmp2)),'k',linewidth=2)
        xlim([-2,7])
        ylim([-2,7])
        suptitle('Mean pressure ='+str(int(mean(prs[ii]))))
        colorbar()

figure(1)
legend(loc=(1.05,0.5),title='Mean pressure \n of other sensor',numpoints=1)
suptitle('Linear regression between instruments')
subplot(121)
xlabel('Salinity at bottom (170db)')
ylabel('Salinity measured by other sensors')
subplot(122)
xlabel('Temperature at bottom (170db)')
ylabel('Temperature measured by other sensors')


#############################################################################
################## Plot method up and make reconstructions #################
#############################################################################

time={}
for ii in sal:
    time[ii]=[datetime.datetime.toordinal(ddd) for ddd in date[ii]]



salcorrvec={}
fig, axarr=subplots(2,1, sharex=True,figsize=(10,6))
for ii in sal:
  if ii!=150:
    timeind=(date[150]<=max(date[ii])) & (date[150]>=min(date[ii]))
    lpass=30

    axarr[0].plot(date[ii][stind:],sal[ii][stind:]-sal[150][timeind])
    axarr[0].plot(date[ii][stind:][::lpass],hrly_ave(sal[ii][stind:]-sal[150][timeind],lpass),'o-')
    mnthcorr=hrly_ave(sal[ii][stind:]-sal[150][timeind],lpass)
    f1=interp1d(month[ii][stind:][::lpass],mnthcorr)
    eachmonthcorr=f1([8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7])
    twoyearcorr=hstack((eachmonthcorr,eachmonthcorr,eachmonthcorr[:2]))
    twoyeartime=hstack((time[150][::lpass],time[150][::lpass][-1]+30))
    f2=interp1d(twoyeartime,twoyearcorr)
    intermonthcorr=f2(time[150][stind:])
    salcorrvec[ii]=intermonthcorr
    axarr[0].plot(date[150][stind:],intermonthcorr,linewidth=2)
    axarr[0].set_ylabel('salinity difference')

    axarr[1].plot(date[ii][stind:],tmp[ii][stind:]-tmp[150][timeind],
                  label=str(ii)+'db')
    axarr[1].plot(date[ii][stind:][::lpass],hrly_ave(tmp[ii][stind:]-tmp[150][timeind],lpass),'o-')
    axarr[1].plot(date[150][stind:],nanmean(tmp[ii][stind:]-tmp[150][timeind])*ones(len(date[150][stind:])),linewidth=2)
    xlabel('date')
    axarr[1].set_ylabel('temperature difference')
    legend(loc=(1.05,0.8),title='Nominal instrument pressure')
savefig('../figures/CF1_recon/CF1_recon_method.png',bbox_inches='tight')



salcorrvec[150]=zeros(len(date[150]))


datevec=date[150]


def reconvec(whichone):
    newvec=zeros(len(datevec))
    newvec[timeind]=whichone[ii][stind:]
    newvec[~timeind]=whichone[150][~timeind]+mean(whichone[ii][stind:]-whichone[150][timeind])
    return newvec

def salrecon(whichone):
    newvec=zeros(len(datevec))
    newvec[timeind]=whichone[ii][stind:]
    newvec[~timeind]=whichone[150][~timeind]+salcorrvec[ii][~timeind]
    return newvec



#############################################################################
# ### For salinity, use a monthly mean historical shift for each instrument and, as before, place them at the mean past position.
#############################################################################

newsal={}
newtmp={}
newprs={}
newmnprs={}
newdate={}
newtime={}
salreg={}
tmpreg={}
prsreg={}
for ii in sal:
    timeind=(date[150]<=max(date[ii])) & (date[150]>=min(date[ii]))

    stind=0
    sal2=sal[150][timeind]
    tmp2=tmp[150][timeind]
    sal_line=poly1d(polyfit(sal2,sal[ii][stind:], 1))
    tmp_line=poly1d(polyfit(tmp2,tmp[ii][stind:], 1))

    salreg[ii]=sal_line(sal[150])
    tmpreg[ii]=tmp_line(tmp[150])

    newsal[ii]=salrecon(sal)
    newtmp[ii]=reconvec(tmp)
    newprs[ii]=reconvec(prs)
    newmnprs[ii]=mean(newprs[ii])*ones(len(datevec))
    newdate[ii]=datevec
    # newtime[ii]=timevec

    figure(1,figsize=(13,10))
    subplot(311)
    plot(newdate[ii],newprs[ii])
    subplot(312)
    plot(newdate[ii],newtmp[ii],label=str(int(round(nanmean(newprs[ii])/50.0)*50.0)))
    subplot(313)
    plot(newdate[ii],newsal[ii])

    figure(figsize=(13,10))
    subplot(311)
    plot(newdate[ii],newprs[ii])
    ylabel('Pressure')
    title('Reconstruction of '+str(int(round(nanmean(newprs[ii])/50.0)*50.0))+'db instrument')
    subplot(312)
    plot(newdate[ii],newtmp[ii],label='using shift method')
    plot(datevec,tmpreg[ii],label='using linear regression')
    legend(loc=(1.05,0.5))
    ylabel('Temperature')
    subplot(313)
    plot(newdate[ii],newsal[ii])
    plot(datevec,salreg[ii])
    ylabel('Salinity')
    xlabel('Date')
    savefig('../figures/CF1_recon/CF1_recon_'+str(ii)+'.png',bbox_inches='tight')

figure(1)
subplot(311)
ylabel('Pressure')
subplot(312)
legend(loc=(1.05,0.5),title='Nominal pressure')
ylabel('Temperature')
subplot(313)
ylabel('Salinity')
xlabel('Date')
# savefig('../notes/figures/CF1_recon.pdf',bbox_inches='tight')
savefig('../figures/CF1_recon/CF1_recon.png',bbox_inches='tight')

#############################################################################
# Save reconstructed CF1 fields to reload in Interp_CTD
#############################################################################

pickle.dump([newdate,newtime,newprs,newmnprs,newsal,newtmp],
            open('../pickles/CF1recon/CF1_recon.pickle','wb'))
