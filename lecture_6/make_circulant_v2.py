import numpy as np
from matplotlib import pyplot as plt

n=500

V=2500
vec=np.arange(n)
mycorr=V-0.5*np.abs(0-vec)

mycorr2=np.append(mycorr,np.flipud(mycorr[1:-1]))

g=np.random.randn(len(mycorr2))
gft=np.fft.rfft(g)
mycorr_ft=np.fft.rfft(mycorr2)
print np.mean(np.abs(np.real(mycorr_ft))), np.mean(np.abs(np.imag(mycorr_ft)))
mycorr_ft=np.real(mycorr_ft)
mycorr_ft[mycorr_ft<0]=0
myps=np.sqrt(mycorr_ft)
mydat=np.fft.irfft(myps*gft)
mydat_use=mydat[:n]
