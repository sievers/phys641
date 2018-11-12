import numpy as np
from matplotlib import pyplot as plt

x=np.arange(-20,20,0.1)
s_true=1
true_amp=20
y_true=true_amp*np.exp(-0.5*x**2/s_true**2)
y=y_true+np.random.randn(len(x))


npt=5000
s_linear=np.linspace(0.1,4,npt)
chisq_linear=0*s_linear
for i in range(npt):
    pred=np.exp(-0.5*(x**2)/s_linear[i]**2)
    lhs=np.dot(pred,pred)
    rhs=np.dot(pred,y)
    amp=rhs/lhs
    delt=y-amp*pred
    chisq_linear[i]=np.sum(delt**2)

var=np.linspace(s_linear[0]**2,s_linear[-1]**2,npt)
chisq_var=0*var
for i in range(npt):
    pred=np.exp(-0.5*(x**2)/var[i])
    lhs=np.dot(pred,pred)
    rhs=np.dot(pred,y)
    amp=rhs/lhs
    delt=y-amp*pred
    chisq_var[i]=np.sum(delt**2)

sigma_marg=np.sum(s_linear*np.exp(-0.5*chisq_linear))/np.sum(np.exp(-0.5*chisq_linear))
var_marg=np.sum(var*np.exp(-0.5*chisq_var))/np.sum(np.exp(-0.5*chisq_var))

print 'fitting for sigma gives me ',sigma_marg,' with data amplitude ',true_amp
print 'fitting for variance gives me ',var_marg,' which works out to a sigma of ',np.sqrt(var_marg)
