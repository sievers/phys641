import numpy as np
from matplotlib import pyplot as plt

def planck_g(x):
    return x**2*np.exp(x)/(np.exp(x)-1)**2
def planck(x):
    return x/(np.exp(x)-1)


dx=0.001
xmax=5
x=np.arange(dx,xmax,dx)
print len(x)

plt.ion()
plt.clf();
plt.plot(x,planck_g(x))
plt.plot(x,planck(x))
plt.xlabel('x')
plt.ylabel('Planck Functions')
plt.legend(['Planck g','B_nu'])
plt.savefig('planck_g.png')

T=2.725
h=6.63e-27
k=1.38e-16

#let's look at some values of Planck functions for typical ground-based CMB windows
nus=[30, 90, 150, 220,270, 350]
for nu in nus:
    x=nu*1e9*h/k/T
    f=planck(x)
    g=planck_g(x)
    print 'at ',nu,' GHz, x is ',x,' Planck function is ',f,' an Planck g is ',g
