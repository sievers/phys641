import numpy as np
from matplotlib import pyplot as plt
x=np.linspace(-5,5,1000)
y_true=4.0*x**3+10*x**2-3*x+8
y=y_true+50*np.random.randn(len(x))

ord=25
a=np.zeros([len(x),ord+1])
a[:,0]=1.0
for i in range(ord):
    a[:,i+1]=a[:,i]*x

ata=np.dot(a.transpose(),a)
rhs=np.dot(a.transpose(),y)
ata_inv=np.linalg.inv(ata)
fitp=np.dot(ata_inv,rhs)
u,s,v=np.linalg.svd(a,0)
utd=np.dot(u.transpose(),y)
sutd=utd/s
fitp2=np.dot(v.transpose(),sutd)
y_pred=np.dot(a,fitp)
y_pred2=np.dot(a,fitp2)
