from pylab import *
import gsw
cd /home/isabela/Documents/projects/OSNAP/OSNAP_hydrogen/Oxygen_Processing
from AAoptode_ctdconvert_funcs import *

from seabird import *
import glob


#######Load the data you want to run through:

bin_list=sort(glob.glob('/home/isabela/Documents/cruises/OSNAP2020_oxygen_AR45/data/ctd/*'))
bin_list
profile = fCNV(bin_list[1])
profile
profile.keys()

####### Run it through:


O2corr=aaoptode_ctdconvert(profile['PSAL'],profile['TEMP'],profile['PRES'],profile['v7'])

profile.keys()

def plot_prof(var):
    plot(var,profile['PRES'])
    gca().invert_yaxis()
    ylabel('pressure [db]')
    xlabel('Oxygen')


plot(var,profile['PRES'])

plot(profile['v7'])

plot_prof(O2corr)
savefig('/home/isabela/Documents/cruises/OSNAP2020_oxygen_AR45/data/ctd/O2profile_sta000.png',bbox_inches='tight')
