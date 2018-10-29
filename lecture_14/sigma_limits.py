import numpy as np
from scipy.special import erfc

npoints=400
niter=1000

plus1=np.zeros(niter)
minus1=np.zeros(niter)
plus2=np.zeros(niter)
minus2=np.zeros(niter)

m1_ind=1-0.5*erfc(-1/np.sqrt(2))
p1_ind=1-0.5*erfc(1/np.sqrt(2))

m2_ind=1-0.5*erfc(-2/np.sqrt(2))
p2_ind=1-0.5*erfc(2/np.sqrt(2))


m1_ind=np.int(m1_ind*npoints)
m2_ind=np.int(m2_ind*npoints)
p1_ind=np.int(p1_ind*npoints)
p2_ind=np.int(p2_ind*npoints)

print 'probs are ',[m2_ind,m1_ind,p1_ind,p2_ind], 'with ',npoints,' independent samples'

for ii in range(niter):
    dat=np.random.randn(npoints)
    dat.sort()
    plus1[ii]=dat[p1_ind]
    plus2[ii]=dat[p2_ind]

    minus1[ii]=dat[m1_ind]
    minus2[ii]=dat[m2_ind]
    
print 'mean/std of -2 sigma is ',np.mean(minus2),np.std(minus2)
print 'mean/std of -1 sigma is ',np.mean(minus1),np.std(minus1)
print 'mean/std of +1 sigma is ',np.mean(plus1),np.std(plus1)
print 'mean/std of +2 sigma is ',np.mean(plus2),np.std(plus2)
