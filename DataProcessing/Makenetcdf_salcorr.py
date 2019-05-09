from aux_funcs import *

nclist=sort(glob.glob(datadir+'MCTD_Data_CF/NetCDF/ilebras_salcorr/*'))

nclist

[date,month,prs,sal,tmp]=pd.read_pickle(open('../pickles/TSdailydic/TS_daily_dic_wcorr_15min.pickle','rb'))

corrlist=[[1,50],[2,100],[2,200],[2,50],[4,100],[4,350],[4,400],[4,50],[5,1000],[5,100],[5,1300],[5,250],
            [5,500],[5,750],[6,150],[6,100],[7,100]]


for ii,nn in enumerate(nclist):

    dset = Dataset(nn,'r+')
    print(nn)
    print(corrlist[ii][0],corrlist[ii][1])

    # print(len(dset['PSAL']))
    # print(len(sal[corrlist[ii][0]][corrlist[ii][1]]))

    figure(figsize=(12,3))
    plot(dset['PSAL'])
    plot(sal[corrlist[ii][0]][corrlist[ii][1]])

    dset['PSAL'][:]=NaN*ones(len(dset['PSAL']))
    dset['PSAL'][:]=sal[corrlist[ii][0]][corrlist[ii][1]]

    dset.close()
