from aux_funcs import *

datadir

sermlist=glob.glob(datadir+'/sermilik2015ctds/*')

sermlist

sermT=[]
sermS=[]
sermP=[]
for ss in sermlist:
    sermdat=io.loadmat(ss)
    sermS=hstack((sermS,sermdat['sal'].flatten()))
    sermT=hstack((sermT,sermdat['temp'].flatten()))
    sermP=hstack((sermP,sermdat['pres'].flatten()))

sermS





dat=pickle.load(open('../pickles/xarray/CF_xarray_notid_1803extrap.pickle','rb'))

CF8=dat.mean(dim='date').sel(distance=95)
CF1=dat.mean(dim='date').sel(distance=0)

salvec=linspace(25,35.1,100)
tmpvec=linspace(-2,10,100)
salmat,tmpmat=meshgrid(salvec,tmpvec)

SA_vec=gsw.SA_from_SP(salvec,zeros(len(salvec)),CFlon[3],CFlat[4])
CT_vec=gsw.CT_from_t(SA_vec,tmpvec,zeros(len(salvec)))
pdenmat=zeros((shape(salmat)))
for ii in range(len(salvec)):
    for jj in range(len(tmpvec)):
        pdenmat[jj,ii]=gsw.sigma0(salvec[ii],tmpvec[jj])


figure(figsize=(5,5))
subplot(121)
plot(sermS,sermP,label='Fjord Water');
# plot(CF8.salinity,CF8.depth,label='Atlantic Water')
# legend()
title('Salinity')
xlim([25,35.1])
ylim([600,0])
axhline(300,color='k')
text(28,150,'cold/fresh')
ylabel('depth [m]')
subplot(122)
gca().set_yticklabels('')
plot(sermT,sermP);
# plot(CF8.temperature,CF8.depth)
title('Potential temperature ($\Theta$) [$^\circ$C]')
text(-1,500,'warm/salty')
ylim([600,0])
tight_layout()
axhline(300,color='k')
savefig('../figures/Sermilik/TSprofiles_Sermilik.pdf')

figure(figsize=(5,5))
subplot(121)
# plot(sermS,sermP,label='Fjord Water');
plot(CF1.salinity,CF1.depth,color='red',linewidth=3)
# legend()
title('Salinity')
# xlim([25,35.1])
ylim([175,0])
ylabel('depth [m]')
subplot(122)
gca().set_yticklabels('')
# plot(sermT,sermP);
plot(CF1.temperature,CF1.depth,color='red',linewidth=3)
title('Potential temperature ($\Theta$) [$^\circ$C]')
ylim([175,0])
tight_layout()
savefig('../figures/Sermilik/TSprofiles_shelf.pdf')


figure(figsize=(5,5))
subplot(121)
# plot(sermS,sermP,label='Fjord Water');
plot(CF8.salinity,CF8.depth,label='Atlantic Water',color='orange',linewidth=3)
text(34.902,250,'mid-depth sal max')
# legend()
title('Salinity')
# xlim([25,35.1])
ylim([2000,0])
ylabel('depth [m]')
subplot(122)
gca().set_yticklabels('')
# plot(sermT,sermP);
plot(CF8.temperature,CF8.depth,color='orange',linewidth=3)
title('Potential temperature ($\Theta$) [$^\circ$C]')
ylim([2000,0])
tight_layout()
savefig('../figures/Sermilik/TSprofiles_AW.pdf')

figure(figsize=(5,5))
subplot(121)
plot(sermS,sermP,label='Fjord Water');
plot(CF8.salinity,CF8.depth,label='Atlantic Water',color='orange',linewidth=3)
plot(CF1.salinity,CF1.depth,label='Shelf Water',color='red',linewidth=3)
legend()
title('Salinity')
xlim([25,35.1])
ylim([600,0])
axhline(300,color='k')
ylabel('depth [m]')
subplot(122)
gca().set_yticklabels('')
plot(sermT,sermP);
plot(CF8.temperature,CF8.depth,label='Atlantic Water',color='orange',linewidth=3)
plot(CF1.temperature,CF1.depth,label='Shelf Water',color='red',linewidth=3)
title('Potential temperature ($\Theta$) [$^\circ$C]')
ylim([600,0])
tight_layout()
axhline(300,color='k')
savefig('../figures/Sermilik/TSprofiles_all.pdf')

figure()
plot(sermS,sermT,'o',label='Fjord Water');
plot(CF8.salinity,CF8.temperature,'o',label='Atlantic Water')
plot(CF1.salinity,CF1.temperature,'ro',label='Shelf Water')
dens=contour(salvec,tmpvec,pdenmat,arange(20.5,29,1),colors='k',zorder=20)
legend(loc=(1.05,0.2))
xlim([25,35.1])
xlabel('salinity')
ylabel('potential temperature ($\Theta$) [$^\circ$C]')
savefig('../figures/Sermilik/TSdiagram.pdf',bbox_inches='tight')
