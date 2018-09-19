import numpy as np
from matplotlib import pyplot as plt

#specify what sort of noise/data we want
n=10000
alpha=-2.0
knee=100.0

#one easy thing to do is generate white noise, then scale its Fourier transform
dat=np.random.randn(n)
datft=np.fft.rfft(dat)
nuvec=np.arange(len(datft))+1 #the plus one is so the bottom frequency isn't zero


filtvec=np.sqrt(1+(nuvec/knee)**alpha)
datft=datft*filtvec
dat_pink=np.fft.irfft(datft)

plt.ion()
plt.clf()
plt.plot(dat_pink)
