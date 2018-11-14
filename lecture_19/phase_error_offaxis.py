import numpy as np
from matplotlib import pyplot as plt
#figure out how far away an off-axis source focuses, and what phase errors look like
D_dish=0.5 
f_ratio=10
f=D_dish*f_ratio;
th_deg=0.5;th=th_deg*np.pi/180 #are we in focus at a distance of th_deg from center?

rdish=D_dish/2.0
a=0.25/f #equation for parabola
lamda=500e-9

x=np.linspace(-rdish,rdish,20)
d=np.linspace(0.95,1.05,100001)*f

dtot=np.zeros([len(d),len(x)])


xf=-d*np.sin(th)   #x/y coordinates of possible focal points
yf=d*np.cos(th)

for  i in range(len(x)):
    x0=x[i]
    y0=a*x0**2
    d1=-(np.cos(th)*y0+np.sin(th)*x0) #distance from infinity to dish
    d2=np.sqrt((yf-y0)**2+(xf-x0)**2) #distance from dish to focus

    dtot[:,i]=d1+d2


plt.ion();
plt.clf();
mystd=np.std(dtot,axis=1)
ii=np.argmin(mystd)
print 'min scatter at angle ',th_deg, ' degrees is ',mystd[ii]
print 'in wavelengths that is ',mystd[ii]/lamda
print 'distance from dish center of focus is ',d[ii]
plt.clf();plt.semilogy(d,mystd)

#for i in range(1,len(x)):
#    plt.plot(d,dtot[:,i]-dtot[:,0])

