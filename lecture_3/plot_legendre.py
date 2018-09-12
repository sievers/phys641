import numpy as np
from matplotlib import pyplot as plt

x=np.linspace(-1,1,10001)
ord=5 #max order of legendre

p=np.zeros([len(x),ord+1])
p[:,0]=1.0
p[:,1]=x
for n in range(1,ord):
    tmp=(2*n+1)*x*p[:,n]-n*p[:,n-1]
    p[:,n+1]=tmp/(n+1)

plt.ion()
plt.clf()
plt.plot(x,p)
plt.savefig('legendre_polys.png')

ata=np.dot(p.transpose(),p)
ee,vv=np.linalg.eig(ata)

cond_pred=(2.0/(2*0+1))/(2.0/(2*ord+1)) #which simplifies to 2*ord+1
print 'condition number is',ee.max()/ee.min(),' compared to expected ',cond_pred
