from aux_funcs import *

[cc,eg]=pickle.load(open('../pickles/geoturndic.pickle','rb'))

def getcorrs(cur):
    print(corrcoef(cur['abs trans'],cur['geo trans'])[1,0])

from scipy.optimize import leastsq

def fitseas(t,data,guess_mean,guess_phase,guess_std,guess_period):

    first_guess=guess_std*np.sin(2*pi*(t+guess_phase)/guess_period) +guess_mean

    optimize_func = lambda x: x[0]*np.sin(2*pi*(t+x[1])/guess_period) + x[2] - data

    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

    # tgrid=linspace(t[0],t[-1],100)

    data_fit = est_std*np.sin(2*pi*(t+est_phase)/guess_period) + est_mean

    return data_fit,est_std

tlen=len(cc['abs trans'])

eg_sin={}
eg_sin['abs trans'],eg_sin['abs amp']=fitseas(arange(tlen),eg['abs trans'],-20,0,10,365.25)
eg_sin['geo trans'],eg_sin['geo amp']=fitseas(arange(tlen),eg['geo trans'],-15,0,10,365.25)
eg_sin['baro trans'],eg_sin['baro amp']=fitseas(arange(tlen),eg['baro trans'],-5,0,10,365.25)

cc_sin={}
cc_sin['abs trans'],cc_sin['abs amp']=fitseas(arange(tlen),cc['abs trans'],-1,0,10,365.25)
cc_sin['geo trans'],cc_sin['geo amp']=fitseas(arange(tlen),cc['geo trans'],-0.5,0,10,365.25)
cc_sin['baro trans'],cc_sin['baro amp']=fitseas(arange(tlen),cc['baro trans'],-0.5,0,10,365.25)

datevec=eg['abs trans'].date

def velcompall(dic,dic2,colc,tit,ylimi):
    figure(figsize=(12,3))
    dic['abs trans'].plot(color=colc,linestyle='-',label='',alpha=0.8)
    dic['geo trans'].plot(color='purple',linestyle='-',label='',alpha=0.8)
    dic['baro trans'].plot(color='grey',linestyle='-',label='',alpha=0.8)
    plot(datevec,dic2['abs trans'],color=colc,linestyle='-',label='absolute transport')
    plot(datevec,dic2['geo trans'],color='purple',linestyle='-',label='geostrophic transport')
    plot(datevec,dic2['baro trans'],color='grey',linestyle='-',label='"barotropic" transport')
    legend()
    title(tit)
    ylim([ylimi,0])
    ylabel('Transport [Sv]')
    savefig('../figures/geocomp/transeas_'+tit+'.png',bbox_inches='tight')


def geoseas(dic,dic2,colc,tit,ylimi):
    figure(figsize=(12,3))
    dic['geo trans'].plot(color=colc,linestyle='-',label='',alpha=0.8)
    plot(datevec,dic2['geo trans'],color=colc,linestyle='-',label='geostrophic transport')
    title(tit)
    legend()
    ylim([ylimi,0])
    ylabel('Transport [Sv]')
    savefig('../figures/geocomp/geoseas_'+tit+'.png',bbox_inches='tight')

def baroseas(dic,dic2,colc,tit,ylimi):
    figure(figsize=(12,3))
    dic['baro trans'].plot(color=colc,linestyle='-',label='',alpha=0.8)
    plot(datevec,dic2['baro trans'],color=colc,linestyle='-',label='barotropic transport')
    title(tit)
    legend()
    ylim([ylimi,0])
    ylabel('Transport [Sv]')
    savefig('../figures/geocomp/baroseas_'+tit+'.png',bbox_inches='tight')

baroseas(eg,eg_sin,egicol,'Slope Current-0.05',-15)
geoseas(eg,eg_sin,egicol,'Slope Current-0.05',-25)
velcompall(eg,eg_sin,egicol,'Slope Current-0.05',-40)

eg_sin['baro amp']
eg_sin['geo amp']
eg_sin['abs amp']


baroseas(cc,cc_sin,ccol,'Coastal Current-0.05',-3)

geoseas(cc,cc_sin,ccol,'Coastal Current-0.1',-3)

velcompall(cc,cc_sin,ccol,'Coastal Current-0.1',-3)

cc_sin['baro amp']
cc_sin['geo amp']
cc_sin['abs amp']

cc.keys()



XXXXXXXXXXXXXXXXXXXXXXXX

# from Spectrum_funcs import *
#
# freqs={}
# ps={}
# psd={}
#
# freqs['cc']={}
# ps['cc']={}
# psd['cc']={}
#
# freqs['eg']={}
# ps['eg']={}
# psd['eg']={}
#
# def calcspec(cur,curstr,sty):
#     freqs[curstr][sty], ps[curstr][sty], psd[curstr][sty] = spectrum4(cur[sty],nsmooth=2)
#     return freqs,ps,psd
#
# def runallspec():
#     freqs,ps,psd = calcspec(cc,'cc','abs trans')
#     freqs,ps,psd = calcspec(cc,'cc','geo trans')
#     freqs,ps,psd = calcspec(cc,'cc','baro trans')
#
#     freqs,ps,psd = calcspec(eg,'eg','abs trans')
#     freqs,ps,psd = calcspec(eg,'eg','geo trans')
#     freqs,ps,psd = calcspec(eg,'eg','baro trans')
#     return freqs,ps,psd
#
# freqs,ps,psd=runallspec()
#
# def compspectra(cur,tit):
#     loglog(freqs[cur]['abs trans'],psd[cur]['abs trans'],color='purple',label='absolute')
#     loglog(freqs[cur]['geo trans'],psd[cur]['geo trans'],color='red',label='geostrophic')
#     loglog(freqs[cur]['baro trans'],psd[cur]['baro trans'],color='blue',label='barotropic')
#     # loglog(freqs[cur]['baro trans'],psdii[cur]['baro trans']+psdii[cur]['geo trans'],color='k',label='abs recon')
#     #summing spectra does not recreate the spectrum of their sum
#     [axvline(1/(365.25/ii),color='k',alpha=0.5) for ii in [1,2,4,8,16,32,64,128]]
#     xticks([1/(365.25/ii) for ii in [1,2,4,8,16,32,64,128]])
#     minorticks_off()
#     gca().set_xticklabels([1,2,4,8,16,32,64,128])
#     xlabel('cycles/year')
#     legend()
#     ylabel('Power spectral density [Sv$^2$ x day]')
#     title(tit)
#     savefig('../figures/geocomp/spec'+cur+'.png')
#
# compspectra('eg','Slope Current')
#
# compspectra('cc','Coastal Current')
#
# loglog(freqs['cc']['abs trans'],psd['cc']['abs trans'],color=ccol,label='coastal current')
# loglog(freqs['eg']['abs trans'],psd['eg']['abs trans'],color=egicol,label='slope current')
# legend()
#
# loglog(freqs['cc']['abs trans'],psd['cc']['abs trans']/psd['cc']['abs trans'][0],color=ccol,label='coastal current')
# loglog(freqs['eg']['abs trans'],psd['eg']['abs trans']/psd['eg']['abs trans'][0],color=egicol,label='slope current')
# legend()
