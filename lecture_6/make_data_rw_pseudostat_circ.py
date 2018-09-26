import numpy as np


n=500

V0=5*n #something a lot larger than we care about
N=np.zeros([n,n])
jvec=np.arange(n)
for i in range(n):
    N[i,:]=V0-0.5*np.abs(i-jvec)

Nplus=np.zeros([n,n])
for i in range(n):
    Nplus[i,:]=V0-0.5*np.abs(i-jvec+n)

Nminus=np.zeros([n,n])
for i in range(n):
    Nminus[i,:]=V0-0.5*np.abs(i-jvec-n)

ii=Nplus>N
N[ii]=Nplus[ii]
ii=Nminus>N
N[ii]=Nminus[ii]

N=0.5*(N+N.transpose())
Nft=np.fft.rfft(N[0,:])
print 'real average is ',np.mean(np.abs(np.real(Nft))),' with imaginay average ',np.mean(np.abs(np.imag(Nft)))
Nft=np.real(Nft)
Nft[Nft<0]=0  #roundoff error has made some of these negative
sig=np.sqrt(Nft)
g=np.random.randn(n)
gft=np.fft.rfft(g)
dat=np.fft.irfft(gft*sig)
assert(1==0)
r=np.linalg.cholesky(N)
g=np.random.randn(n)
dat=np.dot(r,g)
