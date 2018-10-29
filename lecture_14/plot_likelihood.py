import numpy as np
from matplotlib import pyplot as plt

plt.ion()

n=10
x=np.random.randn(n)
sig=np.linspace(0.2,3,1000)
chisq=np.sum(x**2)/sig**2

logdet=n*np.log(sig**2)

loglike=-0.5*chisq-0.5*logdet
plt.clf()
plt.plot(sig,loglike)
plt.xlabel('Standard Deviation')
plt.ylabel('Data Likelihood')
plt.savefig('variance_likelihood.png')

