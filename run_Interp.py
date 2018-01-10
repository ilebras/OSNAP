from Interp_CTD import *
from Interp_AQD_ADCP import *

# Loop through all moorings
# Options such as surface reconstruction option and CF1 reconstruction option will be specified here when I have them

for mm in range(1,9):
    # Interp_CTD(mm)
    Interp_vel(mm)
