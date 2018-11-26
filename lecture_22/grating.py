import numpy as np
from matplotlib import pyplot as plt
plt.ion()

spacing=20
width=15
nn=3000
nline=nn/spacing/2
xx=np.zeros(nn)
prof=np.zeros(nn)
prof[:width]=1.0
prof=prof/prof.sum()
for i in range(nline):
    xx[i*spacing]=1.0

beamft=np.fft.fft(xx)*np.fft.fft(prof)
beam=np.abs(np.fft.ifft(beamft))**2

