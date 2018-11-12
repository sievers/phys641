import numpy as np
from matplotlib import pyplot as plt
plt.ion()
npt=2048
x=np.linspace(0,1-1.0/npt,npt)


nu_range=np.linspace(305,306,11)
win=0.5-0.5*np.cos(2*np.pi*x)
for nu in nu_range:
    y=np.cos(2*np.pi*x*nu)
    yft=np.abs(np.fft.rfft(y))
    yft=yft/np.sqrt(np.sum(y**2))/np.sqrt(npt)
    y2=y*win
    y2ft=np.abs(np.fft.rfft(y2))
    y2ft=y2ft/np.sqrt((np.sum(y2**2)))/np.sqrt(npt)
    print 'max vals are for unwindowed/windowed are',yft.max(),y2ft.max(), ' with freq ',nu
    #print 'norms are ',np.sum(yft**2),np.sum(y2ft**2),
    #plt.clf();plt.semilogy((yft)**2)
    #plt.plot(y2ft**2)
