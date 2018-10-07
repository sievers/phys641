#We want to show how many events we need to see before a Gaussian and Poisson agree within a factor 
#of 2 at both 3- and 5-sigma.  We'll do this by taking the 3/5-sigma levels of a Gaussian, and evaluating the
#difference in the log probability.

import numpy
from matplotlib import pyplot as plt


#For convenience, lets generate log(n!) for a bunch of values of n.
nmax=10000
ivec=numpy.arange(nmax)
ivec[0]=1 #special case since 0! is 1

log_ivec=numpy.log(ivec)
log_fact=numpy.cumsum(log_ivec)
#print a check value
icheck=10
print 'comparison of ',icheck,' factorial is ',numpy.math.factorial(icheck),' which we think is ',numpy.exp(log_fact[icheck])


lamda=numpy.arange(3,nmax-10*numpy.sqrt(nmax))

#sigs=[-5,5,-3,3]
sigs=[5,3,-3]
for nsig in sigs:
    nn=numpy.asarray(numpy.round(lamda+nsig*numpy.sqrt(lamda)),dtype='int64')
    logg=-0.5*(nn-lamda)**2/lamda-0.5*numpy.log(2*numpy.pi*lamda)
    logp=0*logg
    for i in range(len(logp)):
        logp[i]=-lamda[i]+nn[i]*numpy.log(lamda[i])-log_fact[nn[i]]
        
    if logp[0]>logg[0]:
        myind=numpy.min(numpy.where(logp-logg<numpy.log(2.0)))
    else:
        myind=numpy.min(numpy.where(logp-logg>numpy.log(0.5)))
    print 'at ',nsig,' sigma, we need about ',lamda[myind],nn[myind],' events to be good within a factor of 2.'
