
# coding: utf-8

# ## Make sections of T, S, potential density and velocity

# In[34]:


get_ipython().magic('run aux_funcs.ipynb')


# # Want to reload as an xarray!! --> Transport notebook does this already

# ## Load sal, tmp, potential density

# In[18]:


sal={}
tmp={}
pden={}

for moor in range(1,9):
    [sal[moor],
    tmp[moor],
    pden[moor]]=pickle.load(open('../pickles/TSinterp/CF'+str(moor)+'_saltmpinterp_dailyspline_mixedlayer.pickle','rb'))


# ## Load u,v

# In[19]:


u={}
v={}

u_along={}
u_across={}

for moor in range(1,9):
    [u[moor],v[moor]]=pickle.load(open('../pickles/VELinterp/CF'+str(moor)+'_uvinterp.pickle','rb'))
    u_along[moor]=-u[moor]*sin(theta)+v[moor]*cos(theta)
    u_across[moor]=u[moor]*cos(theta)+v[moor]*sin(theta)


# ## Take means, resample
#

# In[20]:


def meanpd(apanda):
    themeanestpanda=pd.DataFrame()
    for moor in range(8,0,-1):
        themeanestpanda[distvec[moor-1]]=mean(apanda[moor],axis=1)
    return themeanestpanda


# In[23]:


meansal=meanpd(sal)
meantmp=meanpd(tmp)
meanpden=meanpd(pden)
mean_ualong=meanpd(u_along)
mean_across=meanpd(u_across)


# In[26]:


datesubvec=['2014-09-15', '2014-12-15', '2015-03-15',
            '2015-06-15', '2015-09-15', '2015-12-15',
            '2016-03-15', '2016-06-15', '2016-09-15']


# In[29]:


def mkresampanda(apanda):
    resampanda={}
    for panel in range(9):
        resampanda[datesubvec[panel]]=pd.DataFrame()
        for moor in range(8,0,-1):
            resampanda[datesubvec[panel]][distvec[moor-1]]=apanda[moor].T.resample('3M',closed='left').iloc[panel]
    return resampanda


# In[31]:


sal.keys()


# In[33]:


len(distvec)


# In[30]:


salresamp=mkresampanda(sal)
tmpresamp=mkresampanda(tmp)
pdenresamp=mkresampanda(pden)
uac_resamp=mkresampanda(u_across)


# ## Plotting

# In[22]:


def onecont(apanda,tit,vrange,coloor,saldp,hlevs):
    ax1=contourf(apanda.columns,apanda.index,apanda.values,vrange,cmap=coloor)
    ax2=contour(apanda.columns,apanda.index,apanda.values,levels=hlevs,colors='k')

    fill_between(bathdist,bathbath,max(saldp)*ones(len(bathbath)),color='k',zorder=22)
    gca().invert_yaxis()
#     colorbar()
    xlabel('distance (km)')
    ylabel('pressure (db)')
    xlim([0,newdistvec[-1]])
    title(tit)

    return ax1,vrange


# In[23]:


def gridnbathy(apanda):
    # this function should extract a "bathymetry", i.e. maximum depth for each mooring in this panda
    # also, extend all fields up and down with the top or bottommost value respectively
    # and add an equivalent profile 5km on either side

    newpanda=apanda.copy()
    newpanda.columns=apanda.columns+5

    newpanda[min(newpanda.columns)-5]=newpanda[min(newpanda.columns)]
    newpanda[max(newpanda.columns)+5]=newpanda[max(newpanda.columns)]

    dpthvec=zeros(len(newpanda.columns))
    for cc in range(len(newpanda.columns)):
        nans_greaterthan100=newpanda.index[(isnan(newpanda.T.iloc[cc])) & (newpanda.index>100)]
        if sum(nans_greaterthan100)==0:
            dpthvec[cc]=max(newpanda.index)
        else:
            dpthvec[cc]=nans_greaterthan100[0]

    for cc in newpanda.columns:
        tester=newpanda[cc]
        tester[isnan(tester) &(tester.index>100)]=tester[~isnan(tester)].iloc[-1]
        tester[isnan(tester) &(tester.index<=100)]=tester[~isnan(tester)].iloc[0]
        newpanda[cc]=tester

    newpanda=newpanda.reindex_axis(sorted(newpanda.columns), axis=1)

    dpthvec=sort(dpthvec)

    return dpthvec,newpanda


# Get and save maximum instrument depth for re-use

# In[24]:


def regridandsmooth(elfieldo,zinty=40,dinty=1):
    saldp,fieldinterp=gridnbathy(elfieldo)
#         print fieldinterp.head()
    #re-grid to smooth a little more
    dold,zold=meshgrid(fieldinterp.columns,fieldinterp.index)
    #new distance grid
    dnew=fieldinterp.columns
    dfunc=interp1d(range(len(fieldinterp.columns)),fieldinterp.columns)
    dnew=dfunc(arange(0,len(fieldinterp.columns)-1+1./dinty,1./dinty))
#     dnew=arange(fieldinterp.columns[0],fieldinterp.columns[-1]+dinty,dinty)
#     print range(len(fieldinterp.columns)),arange(0,len(fieldinterp.columns),1./dinty)
    #new pressure grid
    znew=arange(fieldinterp.index[0],fieldinterp.index[-1],zinty)
    newfield=griddata((dold.flatten(),zold.flatten()),fieldinterp.values.flatten(),(dnew[None,:],znew[:,None]),method='cubic')
    newpanda=pd.DataFrame(index=znew,columns=dnew)
    newpanda[:]=newfield

    return saldp,newpanda


# In[25]:


mean_u_across_smooth=regridandsmooth(mean_u_across)
meanpden_smooth=regridandsmooth(meanpden)
meansal_smooth=regridandsmooth(meansal)
meantmp_smooth=regridandsmooth(meantmp)


# In[29]:

titvec=['August - October','November - January','February - April','May - July']

def fourplot(fieldo,tit,vrange,coloor,savename,hlevs,secondyear=0):
    
    fig=figure(figsize=(12,10))
    for num in range(4):

        if secondyear==1:
            datenum=nm+4
            datename='second'
        else:
            datenum=num
            datename='first'

        subplot(2,2,num+1)
        ax1,vrange=onecont(newpanda,titvec[num],vrange,coloor,saldp,hlevs)

        ylim([1000,0])
        xlim([0,newdistvec[-1]])
        plotinstpos(savename,saldp)


    subplot(221)
    gca().set_xticklabels('')
    xlabel('')
    subplot(222)
    gca().set_xticklabels('')
    gca().set_yticklabels('')
    xlabel('')
    ylabel('')
    subplot(224)
    gca().set_yticklabels('')
    ylabel('')

    cbaxes = fig.add_axes([0.95, 0.25, 0.02, 0.5])
    cbar=colorbar(ax1, cax = cbaxes,ticks=vrange[::2])
    suptitle(tit,fontsize=20)
    savefig('../figures/sections/'+datename+'year_seasonal_'+savename+'_fulltime.png')
    savefig('../figures/sections/'+datename+'year_seasonal_'+savename+'_fulltime.pdf')


# In[27]:


version='dropbot'


# In[ ]:


firstfour(pdenresamp,'Seasonal cycle of potential density at OSNAP Cape Farewell array, 2014 - 2015',linspace(26,27.75,36),cm.BuPu,'pden',[27.25,27.35])


# In[ ]:


firstfour(salresamp,'Seasonal cycle of salinity at OSNAP Cape Farewell array, 2014 - 2015',linspace(32.7,35.1,13),cm.YlGnBu_r,'sal',[34.0,34.3,34.5,34.7])


# In[ ]:


secondfour(salresamp,'Seasonal cycle of salinity at OSNAP Cape Farewell array, 2014 - 2015',linspace(32.7,35.1,13),cm.YlGnBu_r,'sal',[34.0,34.5])


# In[ ]:


firstfour(tmpresamp,'Seasonal cycle of temperature $(^\circ C)$ at OSNAP Cape Farewell array, 2014 - 2015',linspace(-1,8,31),cm.RdYlBu_r,'tmp',[4])


# In[30]:


firstfour(vresamp,'Seasonal cycle of meridional velocity (m/s) at OSNAP Cape Farewell array, 2014 - 2015',arange(-0.6,0.1,0.05),parula_map,'v',[0,-0.2,-0.4])


# In[31]:


fourplot(uac_resamp,'Seasonal cycle of across-track velocity (m/s) at OSNAP Cape Farewell array, 2014 - 2015',arange(-0.6,0.1,0.05),parula_map,'v'+version,[0,-0.2,-0.4])


# In[36]:


def doublecont(apanda1,apanda2,tit,vrange,coloor,saldp,hlevs):
    ax1=contourf(apanda1.columns,apanda1.index,apanda1.values,vrange,cmap=coloor)
    ax2=contour(apanda2.columns,apanda2.index,apanda2.values,levels=hlevs,colors='r',linewidths=2)
    fill_between(bathdist,bathbath,max(saldp)*ones(len(bathbath)),color='k',zorder=22)
    gca().invert_yaxis()
#     colorbar()
    xlabel('distance (km)')
#     ylabel('pressure (db)')
    xlim([0,newdistvec[-1]])
    title(tit)

    ylim([1000,0])
    xlim([0,newdistvec[-1]])
    plotinstpos(savename,saldp)

    return ax1,ax2,vrange


# In[37]:


def firstfour_double(field1,field2,tit,vrange,coloor,collab,savename,hlevs,setem):

    fig=figure(figsize=(7,6))
    for num in range(4):
        subplot(2,2,num+1)
        saldp,newpanda1=regridandsmooth(field1[datesubvec[num]])
        saldp,newpanda2=regridandsmooth(field2[datesubvec[num]])

        ax1,ax2,vrange=doublecont(newpanda1,newpanda2,titvec[num],vrange,coloor,saldp,hlevs)

        ylim([1000,0])
        xlim([0,newdistvec[-1]])
        plotinstpos(savename,saldp)


    subplot(221)
    gca().set_xticklabels('')
    xlabel('')
    subplot(222)
    gca().set_xticklabels('')
    gca().set_yticklabels('')
    xlabel('')
    ylabel('')
    subplot(224)
    gca().set_yticklabels('')
    ylabel('')


    cbaxes = fig.add_axes([0.95, 0.25, 0.02, 0.5])
#     cbaxes2 = fig.add_axes([0.25, -0.005, 0.5, 0.02])
    cbar=colorbar(ax1, cax = cbaxes,ticks=setem,label=collab)
#     cbar=colorbar(ax2, cax = cbaxes2,ticks=hlevs[::2],label='Salinity',orientation='horizontal')
    suptitle(tit,fontsize=12)
    savefig('../figures/sections/firstyear_seasonal_'+savename+'_fulltime.png',bbox_inches='tight')
    savefig('../figures/sections/firstyear_seasonal_'+savename+'_fulltime.pdf',bbox_inches='tight')


# In[39]:


from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


# In[40]:


tmpmap=shiftedColorMap(cm.RdYlBu_r,midpoint=0.7)


# In[ ]:


firstfour_double(tmpresamp,salresamp,'',
                 linspace(-1,8,31),tmpmap,'Temperature $(^\circ C)$','salandtmp_proposal'+version,[34.5],range(-4,10,2))


# In[ ]:


firstfour_double(uac_resamp,salresamp,'',
                arange(-0.6,0.05,0.05),cm.YlGnBu_r,'Across track velocity (m/s)','vacrosswhaline_proposal'+version,[34.5],arange(-0.6,0.1,0.1))


# In[39]:


get_ipython().magic('run aux_funcs.ipynb')


# In[ ]:


newdistvec


# In[34]:


def ADCPtri():
    plot(newdistvec[1],160,'r^',markersize=10,zorder=40)
    plot(newdistvec[2],160,'r^',markersize=10,zorder=40)
    plot(newdistvec[3],160,'r^',markersize=10,zorder=40)
    plot(newdistvec[4],340,'r^',markersize=10,zorder=40)
    plot(newdistvec[5],90,'r^',markersize=10,zorder=40)
    plot(newdistvec[6],90,'r^',markersize=10,zorder=40)
    plot(newdistvec[7],90,'r^',markersize=10,zorder=40)


# In[41]:


fig=figure(figsize=(7,6))

hlevs=[34.5]
savename='tmp'
labx=5

subplot(2,2,1)

saldp,newpanda1=regridandsmooth(tmpresamp[datesubvec[0]])
saldp,sal1=regridandsmooth(salresamp[datesubvec[0]])
tmpax,ax2,vrange=doublecont(newpanda1,sal1,titvec[0],linspace(-1,8,31),tmpmap,saldp,hlevs)

text(labx,900,'Temperature',color='white',zorder=100)

subplot(2,2,2)
saldp,newpanda1=regridandsmooth(tmpresamp[datesubvec[2]])
saldp,sal2=regridandsmooth(salresamp[datesubvec[2]])
tmpax2,ax2,vrange=doublecont(newpanda1,sal2,titvec[2],linspace(-1,8,31),tmpmap,saldp,hlevs)

text(labx,900,'Temperature',color='white',zorder=100)



savename='v'

subplot(2,2,3)
saldp,newpanda3=regridandsmooth(uac_resamp[datesubvec[0]])
vax,ax2,vrange=doublecont(newpanda3,sal1,'',arange(-0.6,0.05,0.05),cm.YlGnBu_r,saldp,hlevs)

text(labx,900,'Velocity',color='white',zorder=100)
ADCPtri()

subplot(2,2,4)
saldp,newpanda3=regridandsmooth(uac_resamp[datesubvec[2]])
vax,ax2,vrange=doublecont(newpanda3,sal2,'',arange(-0.6,0.05,0.05),cm.YlGnBu_r,saldp,hlevs)

text(labx,900,'Velocity',color='white',zorder=100)
ADCPtri()





subplot(221)
gca().set_xticklabels('')
xlabel('')
ylabel('pressure (db)')
subplot(222)
gca().set_xticklabels('')
gca().set_yticklabels('')
xlabel('')
ylabel('')
subplot(223)
ylabel('pressure (db)')
subplot(224)
gca().set_yticklabels('')
ylabel('')


cbaxes = fig.add_axes([0.95, 0.57, 0.02, 0.3])
cbar=colorbar(tmpax2, cax = cbaxes,ticks=range(-4,10,2),label='$^\circ C$')

cbaxes = fig.add_axes([0.95, 0.15, 0.02, 0.3])
cbar=colorbar(vax, cax = cbaxes,ticks=arange(-0.6,0.2,0.2),label='m/s')

savefig('../figures/sections/firstyear_seasonal_tmpandvelwisohal_reverse.pdf',bbox_inches='tight')


# In[ ]:


def fullmean(fieldo,tit,vrange,coloor,savename,hlevs,ylimd,collab):

    fig=figure(figsize=(6,5))

    saldp,newpanda=regridandsmooth(fieldo)

    ax1,vrange=onecont(newpanda,'',vrange,coloor,saldp,hlevs)

    ylim(ylimd,0)


    plotinstpos(savename,saldp)

    cbaxes = fig.add_axes([0.95, 0.25, 0.02, 0.5])
    cbar=colorbar(ax1, cax = cbaxes,ticks=vrange[::2],label=collab)
    suptitle(tit,fontsize=14)

    savefig('../figures/sections/fullmean_'+savename+'_fulltime.png',bbox_inches='tight')
    savefig('../figures/sections/fullmean_'+savename+'_fulltime.pdf',bbox_inches='tight')


# In[ ]:


fullmean(meanv,'Mean meridional velocity at OSNAP Cape Farewell array, 2014 - 2016',arange(-0.8,0.8,0.05),cm.RdBu_r,'vshipcomp',[-0.2,-0.4],800,'m/s')


# In[35]:


def fullmean_proposal(fieldo,tit,vrange,coloor,savename,hlevs,ylimd,collab):

    fig=figure(figsize=(4,3))

    saldp,newpanda=regridandsmooth(fieldo,zinty=50,dinty=1)

    ax1,vrange=onecont(newpanda,'',vrange,coloor,saldp,hlevs)

    ylim(ylimd,0)


    plotinstpos(savename,saldp)

    cbaxes = fig.add_axes([0.95, 0.25, 0.02, 0.5])
    cbar=colorbar(ax1, cax = cbaxes,ticks=vrange[::5],label=collab)
    suptitle(tit,fontsize=12)

    savefig('../figures/sections/fullmean_'+savename+'_fulltime_proposal.png',bbox_inches='tight')
    savefig('../figures/sections/fullmean_'+savename+'_fulltime_proposal.pdf',bbox_inches='tight')


# In[ ]:


fullmean_proposal(mean_u_across,'Mean across-track velocity at OSNAP Cape Farewell array, 2014 - 2016',arange(-0.5,0,0.02),
                  cm.YlGnBu_r,'vacross',arange(0,-0.5,-0.2),1900,'m/s')


# In[ ]:


fullmean_proposal(mean_u_along,'Mean along-track velocity at OSNAP Cape Farewell array, 2014 - 2016',arange(-0.2,0.2,0.01),cm.RdBu_r,'valong',[-0.1,0,0.1],1900,'m/s')


# In[ ]:


fullmean(meantmp,'Mean temperature at OSNAP Cape Farewell array, 2014 - 2016',linspace(-1,8,31),cm.RdYlBu_r,'tmp',[4],1900,'($^\circ$ C)')


# In[ ]:


fullmean(meansal,'Mean salinity at OSNAP Cape Farewell array, 2014 - 2016',linspace(32.7,35.1,13),cm.YlGnBu_r,'sal',[34.0,34.3,34.5,34.7],1900,'')


# In[ ]:


def differentmeans(apanda,aletterfortime,howmany):
    resampanda={}
    for num in range(howmany):
        resampanda[str(num)]=pd.DataFrame()
        for moor in range(7,0,-1):
            resampanda[str(num)][distvec[moor-1]]=apanda[str(moor)].T.resample(aletterfortime,closed='left').iloc[num]
    return resampanda


# In[ ]:


vdaily=differentmeans(u_across,'D',700)


# In[ ]:


def plotmany(fieldo,tit,vrange,coloor,savename,hlevs,howmany,const=0):
    for num in range(howmany):
        figure(figsize=(6,4))
        saldp,newpanda=regridandsmooth(fieldo[str(num+const)])

        ax1,vrange=onecont(newpanda,'',vrange,coloor,saldp,hlevs)


        ylim([1000,0])
        xlim([0,newdistvec[-1]])
        plotinstpos(savename,saldp)
        title(str(num))
        savefig('../figures/vmovie/'+savename+'_'+str(num).zfill(3)+'.png')

#         colorbar(ticks=vrange[::2])


# In[ ]:


plotmany(vdaily,'',arange(-1,0.3,0.05),cm.YlGnBu_r,'vdaily',[0,-0.2,-0.4],650,const=0)
