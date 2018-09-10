import numpy as np
from matplotlib import pyplot as plt
npt=1000
x=np.linspace(-1,1,npt)
y_true=x**3-0.5*x+0.2
noise=0.1
y=y_true+np.random.randn(npt)*noise


for ord in range(1,25):
    a=np.zeros([npt,ord+1])
    a[:,0]=1.0
    for i in range(ord):
        a[:,i+1]=a[:,i]*x
        
    u,s,v=np.linalg.svd(a,0) #0 says to compute compact form of SVD
    utd=np.dot(u.transpose(),y)
    sutd=utd/s
    #always check what's being returned.  help(np.linalg.svd) says the matrix
    # is USV, not USV^T, as is often standard
    fitp=np.dot(v.transpose(),sutd) 

    pred=np.dot(a,fitp)
    rms=np.std(pred-y_true)
    print 'for order=',ord,' rms err is ',rms
plt.ion()
plt.clf();
plt.plot(x,y,'*');
plt.plot(x,pred,'r')
