import numpy as np
import time

n=500

NN=np.zeros([n,n])
t1=time.time()
for i in range(n):
    for j in range(n):
        NN[i,j]=(i+j-np.abs(i-j))/2
t2=time.time()
print 'first way took ',t2-t1

N=np.zeros([n,n])
jvec=np.arange(n)
t1=time.time()
for i in range(n):
        N[i,:]=(i+jvec-np.abs(i-jvec))/2
t2=time.time()
print 'second way took ',t2-t1
print 'average difference is ',np.mean(np.abs(N-NN))
print 'average deviation from symmetry is ',np.mean(np.abs(N-N.transpose()))

#r=np.linalg.cholesky(N)
N=N+1
r=np.linalg.cholesky(N)
g=np.random.randn(n)
dat=np.dot(r,g)
