from aux_funcs import *

dendat=xr.open_dataset(datadir+'OSNAP2016recovery/gridded_CF-OOI/density_gridded_props_cf5-oom_from10m.nc')


samp_str=['1D','3D','1W','1M']
# theta_var={}
# for ss in samp_str:
#     theta_var[ss]=sqrt((dendat.ptmp.resample(date=ss).mean(dim='date')**2).mean(dim='date'))
#
# theta_mgrad={}
# for ss in samp_str:
#     bgtmp=dendat.ptmp.resample(date=ss).mean(dim='date')
#     theta_mgrad[ss]=NaN*bgtmp.mean(dim='date')
#     for dd in range(1,3):
#         theta_mgrad[ss][dd,:]=((bgtmp[dd+1,:,:].mean(dim='date')-bgtmp[dd-1,:,:].mean(dim='date'))/(dendat.distance[dd+1]-dendat.distance[dd-1]))/1e3
#

sal_var={}
for ss in samp_str:
    sal_var[ss]=sqrt((dendat.psal.resample(date=ss).mean(dim='date')**2).mean(dim='date'))

sal_mgrad={}
for ss in samp_str:
    bgsal=dendat.psal.resample(date=ss).mean(dim='date')
    sal_mgrad[ss]=NaN*bgsal.mean(dim='date')
    for dd in range(1,3):
        sal_mgrad[ss][dd,:]=((bgsal[dd+1,:,:].mean(dim='date')-bgsal[dd-1,:,:].mean(dim='date'))/(dendat.distance[dd+1]-dendat.distance[dd-1]))/1e3


sal_var['1D'].plot()
sal_mgrad['1D'].plot()

mixing_length={}
for ss in samp_str:
    mixing_length[ss]=theta_var[ss][1:-1,:]/theta_mgrad[ss]/1e3

for ii in range(2):
    figure()
    for ss in ['1D']:
        plot(mixing_length[ss][ii,:],dendat.den,label=ss)
    legend()

WATER MASS RELATIONSHIPS ARE TOO COMPLEX!!

for ii in range(4):
    figure()
    dendat.ptmp[ii,:,:].plot()
    ylim(27.85,27.5)
    [axhline(dd,color='r') for dd in [d1,d2,d3]]
