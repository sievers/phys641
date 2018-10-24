import numpy as np
from matplotlib import pyplot as plt

n=1000
x=np.arange(n)
x[n/2:]=x[n/2:]-n

dat=np.random.randn(n)
datft=np.fft.rfft(dat)
k=np.arange(len(datft))+0.0
k[0]=0.5
var_k=1e-4+1.0/k**2
datft=datft*np.sqrt(var_k)
dat_back=np.fft.irfft(datft)
