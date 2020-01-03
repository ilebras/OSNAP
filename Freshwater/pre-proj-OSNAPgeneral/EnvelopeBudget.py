from pylab import *

## v1s1+v2s2=v3s3 --> v1s1=v3s3 (since s2=0)
## v1+v2=v3
## v1s1=(v1+v2)s3
## s3/s1=v1/(v1+v2)

v1vec=arange(0.2,3,0.1)

v2=0.009

sratio=v1vec/(v1vec+v2)

def plotfunc():
    plot(v1vec,sratio)
    xlabel('Transport [Sv]')
    ylabel('Ratio of downstream/upstream salinity')


plotfunc()

#Change from mixing 1 Sv of water with salinity 32:
sratio[8]*32

def plotsal():
    # figure(figsize=(4,3))
    plot(v1vec,32*sratio)
    ylabel('Downstream salinity \n if upstream salinity = 32')
    xlabel('Coastal current transport [Sv]')
    title('Mixing 9 mSv in a (portion of a) coastal current')
    savefig('FreshwaterThoughts/VolumeDep.png',bbox_inches='tight')

plotsal()
