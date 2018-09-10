import numpy
from matplotlib import pyplot as plt
import time

dx=0.01
x=numpy.arange(-10,10,dx)
n=len(x)
noise=1+numpy.random.rand(n)
Ninv=1.0/noise**2



x0=0
amp_true=2.0

sig=0.3
template=numpy.exp(-0.5*(x-x0)**2/sig**2)

dat=template*amp_true+numpy.random.randn(n)*noise

snr=numpy.zeros(n)
amp=numpy.zeros(n)
dat_filt=Ninv*dat
denom=(numpy.dot(template,Ninv*template))
rt_denom=numpy.sqrt(denom)
t1=time.time()
big_lhs=numpy.zeros(n)
for i in range(n):
    template_shifted=numpy.exp(-0.5*(x-x[i])**2/sig**2)
    rhs=numpy.dot(template_shifted,dat_filt)
    lhs=numpy.dot(template_shifted**2,Ninv)
    big_lhs[i]=lhs
    snr[i]=rhs/numpy.sqrt(lhs)
    amp[i]=rhs/lhs
t2=time.time()
print 'took ',t2-t1,' seconds to brute force'



t1=time.time()
ft1=numpy.fft.rfft(dat*Ninv)
ft2=numpy.fft.rfft(template)
ft3=numpy.fft.rfft(template**2)
noiseft=numpy.fft.rfft(Ninv)
rhs=numpy.fft.irfft(ft1*ft2)
lhs=numpy.fft.irfft(ft3*noiseft)
amp2=rhs/lhs
snr2=rhs/numpy.sqrt(lhs)
print 'amp maxes are ',amp.max(),amp2.max()
print 'snr_maxes are',snr.max(),snr2.max()

#plt.clf();
#plt.plot(x,snr);
#plt.title('Signal-to-Noise Plot')
#plt.savefig('snr_plot.png')
#
#plt.clf();
#plt.plot(x,amp)
#plt.title('Source Amplitude Plot')
#plt.savefig('amp_plot.png')
#
#plt.clf();
#plt.plot(x,dat)
#plt.title('Raw Data')
#plt.savefig('dat_raw.png')

