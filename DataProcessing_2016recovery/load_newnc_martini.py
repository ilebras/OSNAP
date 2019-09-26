from aux_funcs import *

fnames=glob.glob('../data/MCTD_Data_CF/KMartini/CF4*')

fnames

dat=xr.open_dataset(fnames[1])

dat.FileName

dat.s0.plot()
dat
