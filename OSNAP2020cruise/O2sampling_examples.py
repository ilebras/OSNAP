from aux_funcs_2020 import *

datadir='/home/isabela/Documents/projects/OSNAP/data/Shipboard/netcdf/'
dat=xr.open_dataset(datadir+'AllCTD20_2m_Leahproc.nc')

dat['o2_sbe43']

figdir='/home/isabela/Documents/projects/OSNAP/figures/AR45/'
def O2prof_samp(station,depths):
    figure(figsize=(4,7))
    plot(dat['o2_sbe43'].sel(sta=station),dat['prs'])
    ylabel('pressure [db]')
    gca().invert_yaxis()
    xlabel('O2 concentration [uM]')
    for dd in depths:
        axhline(dd,color='r')
    savefig(figdir+'O2profile_sampling_sta'+str(station)+'.png',bbox_inches='tight')



O2prof_samp(1,[145,1252,1950,2905])
O2prof_samp(3,[1100,2100,2483])
O2prof_samp(86,[192,975,1500,1725])
O2prof_samp(89,[1200,1500,1837])
O2prof_samp(92,[583,950,1310])
