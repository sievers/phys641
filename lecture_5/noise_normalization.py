import numpy as np

N=1500
sig=1.0
dat=np.random.randn(N)*sig

var_real=np.mean(dat**2)

datft=np.fft.fft(dat)
var_fourier=np.real(np.mean(datft*np.conj(datft)))
print 'in real space, my variance is ',var_real,' while in Fourier space it is ',var_fourier


chi2_real=np.sum(dat**2/var_real)

dat_filt=np.fft.ifft(np.fft.fft(dat)/var_fourier)
chi2_fourier=np.real(np.dot(dat,dat_filt))
print 'chi2 values calcualted via real and (wrong) Fourier transform way are ',chi2_real,chi2_fourier

#now make sure the variance is what we'd get if the Fourier transform were normalized to be a strict rotation
var_fourier=np.real(np.mean(datft*np.conj(datft)))/N
dat_filt=np.fft.ifft(np.fft.fft(dat)/var_fourier)
chi2_fourier=np.real(np.dot(dat,dat_filt))
print 'chi2 values calcualted via real and (right) Fourier transform way are ',chi2_real,chi2_fourier
