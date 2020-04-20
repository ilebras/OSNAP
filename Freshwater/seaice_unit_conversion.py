from firstfuncs_1618 import *

#convert from km^3/month to Sv

def convert_to_Sv(num):
    num_Sv=num*1e3**3/30/24/60**2/1e6*1000
    return num_Sv

convert_to_Sv(21)
convert_to_Sv(540)
