import numpy as np


n=500

V0=5*n #something a lot larger than we care about
N=np.zeros([n,n])
jvec=np.arange(n)
for i in range(n):
    N[i,:]=V0-0.5*np.abs(i-jvec)

r=np.linalg.cholesky(N)
g=np.random.randn(n)
dat=np.dot(r,g)
