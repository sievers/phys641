import numpy as np


def simulate_corrnoise(mat,nset=1):
    #if you work through the math, you need
    #to scale gaussian noise by the square root of the eigenvalues
    #then multiply by the eigenvectors.  Note that eigh assumes
    #input matrix is symmetric, and is more stable than eig for our purposes.
    #also, if we want to simulate many sets of data, there's no point finding
    #the eigenvalues/eigenvectors lots of times.  Let nset be the number of simulated datasets you want
    e,v=np.linalg.eigh(mat)
    e[e<0]=0 #make sure we don't have any negative eigenvalues due to roundoff
    n=len(e)
    #make gaussian random variables
    g=np.random.randn(n,nset)
    #now scale them by the square root of the eigenvalues
    rte=np.sqrt(e)
    for i in range(nset):
        g[:,i]=g[:,i]*rte
    #and rotate back into the original space
    dat=np.dot(v,g)
    return dat
    

n=100
nset=100000
mat=np.ones([n,n])+np.eye(n) #make the noise matrix that is one everywhere but 2 along diagonal

dat=simulate_corrnoise(mat,nset)

mat_sim=np.dot(dat,dat.transpose())/nset
print 'RMS error is ',np.std(mat-mat_sim)

