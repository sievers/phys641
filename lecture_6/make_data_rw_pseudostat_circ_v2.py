import numpy as np

n=5000
V0=5*n #something a lot larger than we care about 
#only calculate first row
N=V0-0.5*np.arange(n)
#now paste the covariance back on itself
NN=np.append(N,np.flipud(N[1:-1]))
Nft=np.fft.rfft(NN)
print 'real average is ',np.mean(np.abs(np.real(Nft))),' with imaginay average ',np.mean(np.abs(np.imag(Nft)))
Nft=np.real(Nft)
Nft[Nft<0]=0 
g=np.random.randn(len(NN))
gft=np.fft.rfft(g)
dat_big=np.fft.irfft(gft*np.sqrt(Nft))
dat=dat_big[:n]

ii=Nft>1e-6
x=np.arange(len(Nft))+1
xx=x[ii]
yy=Nft[ii]

i=np.int(np.floor(np.sqrt(len(xx))))
yy_pred=yy[i]*(1.0*xx/xx[i])**-2


