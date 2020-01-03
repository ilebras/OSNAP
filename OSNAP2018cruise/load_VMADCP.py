from AR30_funcs import *

datadir='/home/isabela/Documents/projects/OSNAP/data/Shipboard/AR30_2018/'

secvec=[1,2,3,5]

u={}
v={}
lon={}
lat={}
depth={}
sta={}

for nn in secvec:
    matdat=glob.glob(datadir+'ADCP/vmadcp/vmpro/sect'+str(nn)+'os/*.mat')
    # ('cid', 'O'), ('sectid', 'O'), ('u', 'O'), ('v', 'O'), ('depth', 'O'), ('lon', 'O'), ('lat', 'O'), ('stID', 'O'), ('sampN', 'O'), ('u_dt', 'O'), ('v_dt', 'O'), ('DTmodel', 'O'), ('DTmodelZ', 'O'), ('DTmodelu', 'O'), ('DTmodelv', 'O'), ('trueBT', 'O')
    dat=io.loadmat(matdat[0])['vm_data']
    u[nn]=dat[0][0][2]
    v[nn]=dat[0][0][3]
    lon[nn]=dat[0][0][5].flatten()-360
    lat[nn]=dat[0][0][6].flatten()
    sta[nn]=dat[0][0][7].flatten()
    depth[nn]=dat[0][0][4]

for nn in secvec:
    plot(lon[nn],lat[nn],'o')
def makevec(sta):
    return hstack((sta[1],sta[2],sta[3],sta[5]))

udat=xr.Dataset({'u':(['prs_u','sta'],makevec(u)),'v':(['prs_u','sta'],makevec(v)),'lon':(['sta'],makevec(lon)),'lat':(['sta'],makevec(lat))},coords={'prs_u':gsw.p_from_z(-depth[1][:,0],59.5),'sta':makevec(sta)},)

udat.to_netcdf('OSNAP2018cruise/data/VMADCP_os.nc','w')

for nn in secvec:
    figure()
    pcolor(lon[nn],depth[nn],v[nn])
