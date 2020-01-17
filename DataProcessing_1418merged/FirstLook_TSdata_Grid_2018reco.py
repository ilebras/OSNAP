from firstfuncs_1618 import *

dat={}
for moornum in range(6,8):
    dlist=glob.glob(datadir+'OSNAP2018recovery/mcat_nc/CF'+str(moornum)+'_MCTD*')

    dat[moornum]=xr.open_dataset(dlist[0]).resample(TIME='1H').mean(dim='TIME')
    for dd in dlist[1:]:
        tmp=xr.open_dataset(dd).resample(TIME='1H').mean(dim='TIME')
        dat[moornum]=xr.concat([dat[moornum],tmp],dim='DEPTH')
    dat[moornum]=dat[moornum].transpose('TIME','DEPTH')
    add_SA_CT_PT(dat[moornum])
    dat[moornum]=dat[moornum].drop('PSAL_QC').drop('TEMP_QC').drop('CNDC').drop('CNDC_QC').drop('PRES_QC')

for ii in range(5,8):
    dat[ii].to_netcdf(datadir+'OSNAP2018recovery/Hourly_netcdf/CF'+str(ii)+'_mcat_2018recovery_hourly.nc','w',format='netCDF4')


# for ii in range(1,8):
#     dat[ii].to_netcdf(datadir+'OSNAP2018recovery/Hourly_netcdf/CF'+str(ii)+'_mcat_2018recovery_hourlymerged.nc','w',format='netCDF4')

####These exist for all moorings. Theres a parallel set without "merged" which I made just for moorings 5-7 for Astrid's analysis
dat={}
for ii in range(1,8):
    dat[ii]=xr.open_dataset(datadir+'OSNAP2018recovery/Hourly_netcdf/CF'+str(ii)+'_2018recovery_hourlymerged.nc')


def plot_overview(dat,moornum,xlab):
    f,[ax1,ax2,ax3,ax4,ax5]=subplots(5,1,figsize=(12,15),sharex=True)
    ax1.plot(dat.TIME,dat.PRES.T)
    ax1.invert_yaxis()
    ax1.set_ylabel('pressure [db]')
    ax2.plot(dat.TIME,dat.TEMP.T)
    ax2.set_ylabel('temperature [$^\circ$C]')
    ax3.plot(dat.TIME,dat.PSAL.T)
    ax3.set_ylabel('practical salinity []')
    ax4.plot(dat.TIME,dat.TEMP_QC.T)
    ax4.set_ylabel('temperature QC')
    ax5.plot(dat.TIME,dat.PSAL_QC.T)
    ax5.set_ylabel('salinity QC')
    ax1.set_title('CF'+str(moornum))
    savefig(fig18+'overview_CF'+str(moornum)+xlab+'.png')

def plot_TS(dat,moornum,xlab):
    figure(figsize=(5,4))
    for ii in range(len(dat.DEPTH)):
        plot(dat.PSAL[ii,:],dat.TEMP[ii,:],'.')
    xlabel('salinity')
    ylabel('temperature')
    title('CF'+str(moornum))
    xlim(32,35.5)
    ylim(-2,10)
    savefig(fig18+'TS_CF'+str(moornum)+xlab+'.png')

for mm in range(1,8):
    plot_overview(dat[mm],mm,'')

for mm in range(1,5):
    plot_TS(dat[mm],mm,'')

for mm in range(5,8):
    plot_TS(dat[mm],mm,'')
    xlim(34.5,35.2)
    ylim(1,8)
    savefig(fig18+'TSzoom_CF'+str(mm)+'.png')

for dd in dat:
    print(dd)

def plot_TSall():
    figure(figsize=(5,4))
    for dd in dat:
        plot(dat[dd].PSAL,dat[dd].TEMP,'.',color='C'+str(dd),label=str(dd))
    legend()
    xlabel('salinity')
    ylabel('temperature')
    xlim(32,35.5)
    ylim(-2,10)
    savefig(fig18+'TSall.png')

plot_TSall()



dat_daily={}
for ii in range(1,8):
    dat_daily[ii]=dat[ii].resample(TIME='1D').mean(dim='TIME')





for ii in range(1,8):
    add_SA_CT_PT(dat_daily[ii])


for ii in range(1,8):
    dat_daily[ii].to_netcdf(datadir+'OSNAP2018recovery/Daily_netcdf/CF'+str(ii)+'_mcat_2018recovery_daily.nc','w',format='netCDF4')

def plot_overview(dat,moornum,xlab):
    f,[ax1,ax2,ax3,ax4,ax5]=subplots(5,1,figsize=(12,15),sharex=True)
    ax1.plot(dat.TIME,dat.PRES)
    ax1.invert_yaxis()
    ax1.set_ylabel('pressure [db]')
    ax2.plot(dat.TIME,dat.TEMP)
    ax2.set_ylabel('temperature [$^\circ$C]')
    ax3.plot(dat.TIME,dat.PSAL)
    ax3.set_ylabel('practical salinity []')
    ax4.plot(dat.TIME,dat.TEMP_QC)
    ax4.set_ylabel('temperature QC')
    ax5.plot(dat.TIME,dat.PSAL_QC)
    ax5.set_ylabel('salinity QC')
    ax1.set_title('CF'+str(moornum))
    savefig(fig18+'overview_CF'+str(moornum)+'_'+xlab+'.png')

for mm in range(1,8):
    plot_overview(dat_daily[mm],mm,'daily')


def plot_TS(dat,moornum,xlab):
    figure(figsize=(5,4))
    for ii in range(len(dat.DEPTH)):
        plot(dat.PSAL[:,ii],dat.TEMP[:,ii],'.')
    xlabel('salinity')
    ylabel('temperature')
    title('CF'+str(moornum))
    xlim(32,35.5)
    ylim(-2,10)
    savefig(fig18+'TS_CF'+str(moornum)+'_'+xlab+'.png')


for mm in range(1,5):
    plot_TS(dat_daily[mm],mm,'daily')

for mm in range(5,8):
    plot_TS(dat_daily[mm],mm,'daily')
    xlim(34.5,35.2)
    ylim(1,8)
    savefig(fig18+'TSzoom_CF'+str(mm)+'_daily.png')

def plot_TSall_daily():
    figure(figsize=(5,4))
    for dd in dat:
        plot(dat_daily[dd].PSAL,dat_daily[dd].TEMP,'.',color='C'+str(dd),label=str(dd))
    xlabel('salinity')
    ylabel('temperature')
    xlim(32,35.5)
    ylim(-2,10)
    savefig(fig18+'TSall_daily.png')
    xlim(34.5,35.2)
    ylim(1,8)
    savefig(fig18+'TSzoom_all_daily.png')
    ylim(3.2,5)
    xlim(34.865,34.98)
    savefig(fig18+'TSzoommore_all_daily.png')

plot_TSall_daily()
