import numpy as np
x=np.arange(1000)
n=len(x)
#make the signal vector we want
y=np.exp(-0.5*(x-500.0)**2/50.0**2)

#next, make a matrix that has delta x in it.
#lots of ways to do this, but this is one which is fast
xmat=np.repeat([x],len(x),axis=0)
dx=xmat-xmat.transpose()

for a in [0.1,0.5,0.9]:
    for sig in [5.0,50.0,500.0]:
        N=a*np.exp(-0.5*dx**2/sig**2)+(1-a)*np.eye(n)
        Ninv=np.linalg.inv(N)
        lhs=np.dot(y,np.dot(Ninv,y))
        err=1/np.sqrt(lhs)
        print 'for a ', a,' and sigma ',sig,' error is ',err



#part b:  the worst error bars are for a=0.9 and sigma=50.  The reason is that
#the correlated noise in this case looks just like the signal we're searching for,
#which also has 50 sample width.  So, it's going to be very hard to tell a true source
#from a noise fluctuation.  The best case is a=0.9 and sigma=500.  For large values of a,
#most of the noise is in the correlations, not the white part.  When the correlation length is very long, 
#then the noise has a hard time making short bumps, so our source is going to stand out.  

