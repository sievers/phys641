import numpy as np

npt=100 #number of points to use
nset=100000  #number of independent sets to use

sig=10
amp=5

x=np.arange(npt)
x=x-np.mean(x)

y_true=amp*np.exp(-0.5*(x**2/sig**2))
A=y_true

fit_amp=np.zeros(nset)
fit_var=np.zeros(nset)
for i in range(nset):
    d=y_true+np.random.randn(npt)
    mynoise=np.std(d)
    Ninv=np.eye(npt)/mynoise**2
    lhs=np.dot(A,np.dot(Ninv,A))
    rhs=np.dot(A,np.dot(Ninv,d))
    fit_amp[i]=rhs/lhs
    fit_var[i]=1.0/lhs

print 'mean recovered amplitude is ',fit_amp.mean()

wts=1/fit_var
tot_mean=np.sum(wts*fit_amp)/np.sum(wts)
wt_tot=np.sum(wts)
err_tot=1/np.sqrt(wt_tot)
nsig=(1-tot_mean)/err_tot
print 'mean recovered weighted amplitude is ',tot_mean,' plus/minus ',err_tot, 'for an error of ',nsig,' standard deviations'
