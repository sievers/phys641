import numpy as np

n=50 #how many random variables to create
niter=100000   #how many trials to average over to confirm our correlation matrix is the same

#first, create a positive-definite matrix
true_cov=np.random.randn(n,n)
true_cov=np.dot(true_cov.transpose(),true_cov)  #this is guaranteed to be positive-definite
true_cov=true_cov+np.eye(n)  #lets just push the eigenvalues a bit away from zero, just to be safe.

L=np.linalg.cholesky(true_cov)

g=np.random.randn(n,niter)
d=np.dot(L,g)
mycov=np.dot(d,d.transpose())/niter

print 'average error is ',np.mean(np.abs(mycov-true_cov)), ' out of an average value of ',np.mean(np.abs(true_cov))
