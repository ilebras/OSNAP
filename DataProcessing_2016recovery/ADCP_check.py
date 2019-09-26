from aux_funcs import *

for mm in range(1,2):

    moorname='cf'+str(mm)

    adcplist=glob.glob(datadir+'ADCP_Data_CF/*'+moorname+'*Final1_ilebras.mat')

    adcplist

    dat=io.loadmat(adcplist[0])

    figure(figsize=(15,5))
    subplot(121)
    plot(dat['u'],dat['z'],'o');
    ylabel('prs [db]')
    title('u')
    gca().invert_yaxis()
    subplot(122)
    plot(dat['v'],dat['z'],'o');
    title('v')
    gca().invert_yaxis()
    savefig('../../figures/ADCP_raw/uvprof_'+moorname+'.png')

    for ii in range(shape(dat['u'])[0]):
        figure()
        subplot(311)
        plot(dat['z'][ii,:]);
        ylabel('prs')
        subplot(312)
        ylabel('u')
        plot(dat['u'][ii,:]);
        subplot(313)
        plot(dat['v'][ii,:]);
        ylabel('v')
        savefig('../../figures/ADCP_raw/'+moorname+'_uvtseries_'+str(ii)+'.png')
