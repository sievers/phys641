import numpy as np
from matplotlib import pyplot as plt


n=5000
sigs=[5.0, 50.0, 500.0]
a=0.9

x=np.arange(2*n)
x[x>n]=x[x>n]-2*n

#make a gaussian
template=np.exp(-0.5*(x**2)/(50.0**2))

plt.ion()
plt.clf();
vecs=[None]*len(sigs)
g=np.random.randn(len(x))
gft=np.fft.rfft(g)
for i in range(len(sigs)):
    sig=sigs[i]
    mycorr=a*np.exp(-0.5*(x**2/sig**2))
    mycorr[0]+=(1-a)
    corrft=np.fft.rfft(mycorr)
    corrft=np.real(corrft)
    print 'for sig ',sig,' total expected noise power is ',corrft.sum()
    dat_sim=np.fft.irfft(gft*np.sqrt(corrft))
    vecs[i]=dat_sim
    plt.loglog(np.abs(corrft))

    
template_ft=np.fft.rfft(template)
template_ps=np.abs(template_ft)**2
plt.plot(template_ps/template_ps[0]*101,':')
ax=plt.axis()
ax=np.asarray(ax)
ax[2]=1e-4
plt.axis(ax)
