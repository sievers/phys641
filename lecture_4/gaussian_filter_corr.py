import numpy
from matplotlib import pyplot as plt
import time

dx=0.01
noise=0.5
Ninv=1.0/noise**2
x=numpy.arange(-10,10,dx)
n=len(x)


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
for i in range(n):
    template_shifted=numpy.exp(-0.5*(x-x[i])**2/sig**2)
    rhs=numpy.dot(template_shifted,dat_filt)
    snr[i]=rhs/rt_denom
    amp[i]=rhs/denom
t2=time.time()
dt_slow=t2-t1
print 'took ',dt_slow,' seconds to brute force'
t1=time.time()
ft1=numpy.fft.rfft(dat)
ft2=numpy.fft.rfft(template)
rhs=numpy.fft.irfft(numpy.conj(ft1)*ft2*Ninv)


snr_conv=rhs/rt_denom
amp_conv=rhs/denom
t2=time.time()
dt_fast=t2-t1
print 'took ',dt_fast,' seconds to do FFT convolutions, ', dt_slow/dt_fast,' times faster than brute force'
print 'max snrs are ',snr.max(),snr_conv.max()


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

