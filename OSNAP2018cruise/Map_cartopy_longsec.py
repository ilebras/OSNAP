from map_funcs import *
import cartopy
import cartopy.crs as ccrs
from AR30_funcs import *
ctd=xr.open_dataset('OSNAP2018cruise/data/CTD_2m_Leahproc.nc')



central_lon, central_lat = -30, 60
extent = [-55, -35, 58, 62]

# Plot the long sections only: 1,2,3,5:
seclist=['section 1','section 2','section 3','section 5']



def map_setup():
    figure(figsize=(14,10))
    ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))
    ax.set_extent(extent)
    ax.coastlines(resolution='50m')
    for ss in seclist:
        ax.plot(ctd.lon[seclab[ss]],ctd.lat[seclab[ss]],'o',label=ss, transform=ccrs.PlateCarree())


map_setup()

#Need that etopo1!!
