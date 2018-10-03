import numpy as np
from numpy import pi


msun=2e33
G=6.67e-8
c=3e10

m1=1.4*msun
m2=m1
M=m1+m2


#first, find max frequency we're going to hit.
isco=6*G*M/c**2
omega_sqr=G*M/isco**3
fmax=np.sqrt(omega_sqr)/2/pi


#now let's integrate df/dt=96/f pi^(8/3) (Gm_chirp/c^3)^(5/3) f^(11/3)
#solution = df*f^(-11/3)=a*dt
#-3/8 f^(-8/3) = a t + k
#f=(-3(at+k)/8)^(-3/8)
#if f_crit happens at t=0, then 
#k=-3/8f_crit^(-8/3)

k=-3.0/8.0*fmax**(-8.0/3.0)

mc=(m1*m2)**(3.0/5.0)/(m1+m2)**(1.0/5.0)
mfac=G*mc/c**3
a=96.0/5.0*pi**(8.0/3.0)*mfac**(5.0/3.0)

dt=1e-4
t0=-20
t=np.arange(t0,0,dt)

f=(-3.0*(a*t+k)/8.0)**(-3.0/8.0)







