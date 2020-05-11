from firstfuncs_1618 import *

#convert from km^3/month to Sv

def convert_to_Sv(num):
    num_Sv=num*1e3**3/30/24/60**2/1e6*1000
    return num_Sv

convert_to_Sv(21)
convert_to_Sv(540)



def from_km3yr_to_mSv(var):
    ans=var/1000*31.7
    return ans

from_km3yr_to_mSv(1000)
