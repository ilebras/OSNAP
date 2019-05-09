#################################################################################
#################################################################################
#################################################################################
#################################################################################
## Testing out atom while inspecting density inversions
## and effects on TS space of interpolated TS
#################################################################################
#################################################################################
#################################################################################
#################################################################################

from aux_funcs import *

moormin=3
moormax=4

#################################################################################
#################################################################################
##       Load
#################################################################################
#################################################################################


sal={}
tmp={}
pden={}


for ii in range(moormin,moormax):
    moorname='CF'+str(ii)
    sal[ii]={}
    tmp[ii]={}
    pden[ii]={}

    sal[ii][0],tmp[ii][0],pden[ii][0]=pd.read_pickle(open('../pickles/TSinterp/'+moorname+'_saltmpinterp_atom_notid.pickle','rb'))
    sal[ii][1],tmp[ii][1],pden[ii][1]=pd.read_pickle(open('../pickles/TSinterp/'+moorname+'_saltmpinterp_atom.pickle','rb'))

#################################################################################
#################################################################################
##       Inversion histograms
#################################################################################
#################################################################################


for ii in range(moormin,moormax):
    moorname='CF'+str(ii)

    figure()

    tid=pd.to_numeric(diff(pden[ii][1].values,axis=0).flatten())
    notid=pd.to_numeric(diff(pden[ii][0].values,axis=0).flatten())

    hist(tid[~isnan(tid)],bins=arange(-0.02,0.06,0.002),color='red')
    hist(notid[~isnan(notid)],bins=arange(-0.02,0.06,0.002),color='blue',alpha=0.5)

    title(moorname)
    xlim([-0.03,0.07])
    ylabel('Number in each bin')
    xlabel('Vertical difference of Density profile')
    axvline(0,color='k')

    savefig('../figures/invercheck/'+moorname+'_inv_hist.png')

#################################################################################
#################################################################################
##       TS density plots
#################################################################################
#################################################################################


##     First need to group all together

for ii in range(moormin,moormax):
    if ii>1:
        saltot=hstack((saltot,vstack((sal[ii][0].values.flatten(),sal[ii][1].values.flatten()))))
        tmptot=hstack((tmptot,vstack((tmp[ii][0].values.flatten(),tmp[ii][1].values.flatten()))))
    else:
        saltot=vstack((sal[ii][0].values.flatten(),sal[ii][1].values.flatten()))
        tmptot=vstack((tmp[ii][0].values.flatten(),tmp[ii][1].values.flatten()))

##     Monthly version

for ii in range(moormin,moormax):
    if ii>1:
        salmonth=hstack((saltot,vstack((sal[ii][0].T.resample('M').mean().values.flatten(),sal[ii][1].T.resample('M').mean().values.flatten()))))
        tmpmonth=hstack((tmptot,vstack((tmp[ii][0].T.resample('M').mean().values.flatten(),tmp[ii][1].T.resample('M').mean().values.flatten()))))
    else:
        salmonth=vstack((sal[ii][0].T.resample('M').mean().values.flatten(),sal[ii][1].T.resample('M').mean().values.flatten()))
        tmpmonth=vstack((tmp[ii][0].T.resample('M').mean().values.flatten(),tmp[ii][1].T.resample('M').mean().values.flatten()))


def hexTS(salvers,tmpvers,namevers):
    figure()
    hexbin(salvers, tmpvers, bins='log', cmap=cm.hot_r)
    colorbar()
    xlim([31,36])
    ylim([-2,10])
    pcont=contour(salvec,tmpvec,pdenmat,colors='k')
    clabel(pcont)
    xlabel('salinity')
    ylabel('potential temperature ($^\circ C$)')
    savefig('../figures/TS/TS_total_log_interp_'+namevers+'.png')
    savefig('../figures/TS/TS_total_log_interp_'+namevers+'.pdf')

hexTS(saltot[0,:], tmptot[0,:],'notid')
hexTS(saltot[1,:], tmptot[1,:],'wtid')

hexTS(salmonth[0,:], tmpmonth[0,:],'monthly_notid')
hexTS(salmonth[1,:], tmpmonth[1,:],'monthly_wtid')

#################################################################################
##        INDIVIDUAL MOORING TS PLOTS W+/OUT TIDBITS
#################################################################################

for ii in range(moormin,moormax):
        figure()
        scatter(sal[ii][1].values.flatten(),tmp[ii][1].values.flatten(),
                c=pden[ii][1].values.flatten())
        colorbar()
        contour(salvec,tmpvec,pdenmat,colors='k')
        ylim([-2,10])
        xlabel('salinity')
        ylabel('temperature')
        title('CF'+str(ii)+' with tidbits')
        savefig('../figures/invercheck/CF'+str(ii)+'_TS_wtid.png')



for ii in range(moormin,moormax):
    figure()
    scatter(sal[ii][0].values.flatten(),tmp[ii][0].values.flatten(),
            c=pden[ii][0].values.flatten())
    colorbar()
    contour(salvec,tmpvec,pdenmat,colors='k')
    ylim([-2,10])
    xlabel('salinity')
    ylabel('temperature')
    title('CF'+str(ii)+' no tidbits')

    savefig('../figures/invercheck/CF'+str(ii)+'_TS_notid.png')


#################################################################################
#################################################################################
##        PLOT PROFILES
#################################################################################
#################################################################################

def plotprofs(field,fieldname,savename):
        figure()
        for ii in range(moormin,moormax):
            moorname='CF'+str(ii)

            figure(figsize=(10,4))
            for jj in range(2):
                subplot(1,2,jj+1)
                plot(field[ii][jj].values,field[ii][jj].index)
                gca().invert_yaxis()
                suptitle(moorname+' '+fieldname)

            subplot(121)
            title('no tidbits')
            ylabel('pressure')
            subplot(122)
            title('with tidbits')
            savefig('../figures/profiles/prof_'+savename+'_CF'+str(ii)+'_both.png')
            savefig('../figures/profiles/pdfprof_'+savename+'_CF'+str(ii)+'_both.pdf')

# plotprofs(sal,'Salinity','sal')
#
# plotprofs(tmp,'Temperature','tmp')
#
# plotprofs(pden,'Potential Density','pden')


#################################################################################
#################################################################################
##       PLOT ONLY THE PROFILES THAT INCLUDE INVERSIONS
#################################################################################
#################################################################################


def plotprofs_inv(field,fieldname,savename):
    figure()
    for ii in range(moormin,moormax):
        moorname='CF'+str(ii)


        figure(figsize=(10,4))
        for jj in range(2):
            pdendiff=diff(pden[ii][jj],axis=0)
            subplot(1,2,jj+1)
            cc=0
            for tt in range(shape(field[ii][jj])[1]):
                if sum(pdendiff[:,tt]<0)>0:
                    plot(field[ii][jj].iloc[:,tt],field[ii][jj].index)
                    cc+=1
            if jj==0:
                title('no tidbits:'+' '+str(cc)+' profiles')
            else:
                title('with tidbits:'+' '+str(cc)+' profiles')
            gca().invert_yaxis()
        suptitle(moorname+' '+fieldname)
        subplot(121)
        ylabel('pressure')

        savefig('../figures/profiles/prof_'+savename+'_CF'+str(ii)+'.png')
        savefig('../figures/profiles/pdfprof_'+savename+'_CF'+str(ii)+'.pdf')


plotprofs_inv(pden,'Potential Density Inversions','pden_inv_both')


plotprofs_inv(sal,'Salinity of Inversions','sal_inv_both')

plotprofs_inv(tmp,'Temperature of Inversions','tmp_inv_both')
